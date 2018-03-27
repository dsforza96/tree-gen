#include "yocto/yocto_gl.h"
#include "yocto/yocto_gl.cpp"
#include "voro++-0.4.6/src/voro++.hh"

using namespace ygl;

const auto e = 2.71828184f;     // esponente per sommatoria dei raggi
const auto r0 = 0.01f;          // raggio iniziale
const auto eps = 0.02f;

/* Genera in modo randomico i punti di attrazione e crea il
   diagramma di Voronoi */
voro::container throw_darts(int N, const vec2f& p0, const vec2f& p1, const vec2f& t0, const vec2f& t1)
{
    auto points = vector<vec3f>();
    auto bbox = bbox3f();
    auto rng = init_rng(time(nullptr));

    for (auto i = 0; i < N; i++)
    {
        auto y = next_rand1f(rng, p0.x, p1.x);
        auto r = eval_bezier_cubic(p1, t1, t0, p0, y).y;
        auto p = sqrtf(next_rand1f(rng, 0, r * r));
        auto t = next_rand1f(rng, 0, 2 * pif);

        points.push_back({cosf(t) * p, y, sinf(t) * p});

        bbox += points.back();
    }

    auto vorodiag = voro::container(bbox.min.x - eps, bbox.max.x + eps, bbox.min.y - eps, bbox.max.y + eps,
                                    bbox.min.z - eps, bbox.max.z + eps, N / 5, N / 5, N / 5, false, false, false, 8);

    // aggiungo tutti i punti di attrazione
    for (auto i = 0; i < N; i++)
        vorodiag.put(i, points[i].x, points[i].y, points[i].z);

    return vorodiag;
}

void add_branch(vec3f node, int node_id, float di, float D, const voro::container& voro_attr,
                const voro::container& voro_nodes, unordered_set<int>& dead_attr,
                unordered_set<int>& computed_attr, unordered_set<int>& dead_nodes,
                vector<pair<vec3f, int>>& new_nodes)
{
    auto attr_id = 0;
    auto x = 0.0, y = 0.0, z = 0.0, reder = 0.0;

    auto node_loop = voro::c_loop_subset(voro_nodes);
    auto attr_loop = voro::c_loop_subset(voro_attr);
    attr_loop.setup_sphere(node.x, node.y, node.z, di, true);

    auto sum = vec3f{0, 0, 0};
    auto summed = false;

    if (attr_loop.start())
        do
        {
            attr_loop.pos(attr_id, x, y, z, reder);

            if (dead_attr.count(attr_id) || computed_attr.count(attr_id))
                continue;

            auto attr = vec3f{(float) x, (float) y, (float) z};

            node_loop.setup_sphere(x, y, z, length(attr - node) - eps, true);
            if (node_loop.start())
                continue;

            computed_attr.insert(attr_id);

            sum += normalize(attr - node);
            summed = true;
        } while (attr_loop.inc());

    if (summed)
    {
        auto new_node = node + D * normalize(sum);
        new_nodes.push_back({new_node, node_id});
    }
    else
        dead_nodes.insert(node_id);
}

void kill_points(vec3f node, float  dk, const voro::container& voro_attr, unordered_set<int>& dead_attr)
{
    auto attr_loop = voro::c_loop_subset(voro_attr);
    attr_loop.setup_sphere(node.x, node.y, node.z, dk, true);

    auto attr_id = 0;
    auto x = 0.0, y = 0.0, z = 0.0, reder = 0.0;

    if (attr_loop.start())
        do
        {
            attr_loop.pos(attr_id, x, y, z, reder);
            dead_attr.insert(attr_id);
        } while (attr_loop.inc());
}

// Crescita dell'albero
vector<vec3f> grow(int iter_num, int N, float D, float di, float dk,
                   const voro::container& voro_attr, vector<int>& parents)
{
    auto positions = vector<vec3f>();
    parents = vector<int>();

    positions.push_back({0, 0, 0});
    parents.push_back(0);

    auto attr_loop = voro::c_loop_subset(voro_attr);

    auto voro_nodes = voro::container(voro_attr.ax, voro_attr.bx, min(voro_attr.ay, -eps), voro_attr.by,
                                      voro_attr.az, voro_attr.bz, N / 5, N / 5, N / 5, false, false, false, 8);
    voro_nodes.put(0, 0.0f, 0.0f, 0.0f);
    auto nodes_loop = voro::c_loop_all(voro_nodes);

    auto dead_attr = unordered_set<int>();
    auto computed_attr = unordered_set<int>();
    auto dead_nodes = unordered_set<int>();
    auto new_nodes = vector<pair<vec3f, int>>();

    auto threads = vector<std::thread>();

    auto node_id = 0, attr_id = 0;
    auto x = 0.0, y = 0.0, z = 0.0, reder = 0.0;

    // uccidi gli attraction points vicini alla radice
    attr_loop.setup_sphere(x, y, z, dk, true);

    if (attr_loop.start())
        do
        {
            attr_loop.pos(attr_id, x, y, z, reder);
            dead_attr.insert(attr_id);
        } while (attr_loop.inc());

    // loop di crescita
    for (auto i = 0; i < iter_num; i++)
    {
        if (dead_attr.size() == N)
            break;

        computed_attr.clear();

        if (nodes_loop.start())
            do
            {
                nodes_loop.pos(node_id, x, y, z, reder);

                if (dead_nodes.count(node_id))
                    continue;

                auto node = vec3f {(float) x, (float) y, (float) z};

                threads.push_back(std::thread(add_branch, node, node_id, di, D,
                                              std::ref(voro_attr), std::ref(voro_nodes), ref(dead_attr),
                                              ref(computed_attr), ref(dead_nodes), ref(new_nodes)));
            } while (nodes_loop.inc());

        for (auto& t : threads)
            t.join();
        threads.clear();

        if (new_nodes.empty())
            break;

        while (!new_nodes.empty())
        {
            auto new_node = new_nodes.back().first;
            voro_nodes.put(positions.size(), new_node.x, new_node.y, new_node.z);
            positions.push_back(new_node);
            parents.push_back(new_nodes.back().second);

            threads.push_back(std::thread(kill_points, new_node, dk, std::ref(voro_attr), ref(dead_attr)));

            new_nodes.pop_back();
        }

        for (auto& t : threads)
            t.join();
        threads.clear();

        if (i % 10 == 9 || !i)
            log_info("{} nodes added after {} iterations...", positions.size(), i + 1);
    }

    return positions;
}

void make_cylinder(shape* tree, const vec3f& node, const vec3f& p_node, float r)
{
    auto points = (int) tree->pos.size();

    auto axis = node - p_node;
    auto h = length(axis);
    auto f = make_frame_fromz(p_node, axis);

    for (auto i = 0; i <= 16; i++)
        for (auto j = 0; j <= 16; j++)
    {
            auto u = (float) i / 16;
            auto v = (float) j / 16;

            auto c = transform_point(f, {0, 0, v * h});

            tree->pos.push_back(transform_point(f, {cosf(u * 2 * pif) * r, sinf(u * 2 * pif) * r, v * h}));
            tree->norm.push_back(normalize(tree->pos.back() - c));
            tree->texcoord.push_back({u, v});

            if (i != 16 && j != 16)
                tree->quads.push_back(
                        {points + i * (16 + 1) + j, points + (i + 1) * (16 + 1) + j,
                         points + (i + 1) * (16 + 1) + j + 1, points + i * (16 + 1) + j + 1});
        }
}

// Crea la shape dell'albero con i quad dei cilindri
shape* draw_tree(const vector<vec3f> positions, const vector<int>& parents)
{
    auto shp = new shape{"tree"};

    auto rad = vector<float>(positions.size(), 0.0f);

    for (auto i = (int) positions.size() - 1; i > 0; i--)
    {
        auto pos = positions[i];

        auto par = parents[i];
        auto ppos = positions[par];

        auto r = rad[i] == 0.0f ? r0 : pow(rad[i], 1/e);
        rad[par] += pow(r, e);

        make_cylinder(shp, pos, ppos, r);
    }

    return shp;
}

int main(int argc, char** argv)
{
    auto parser = make_parser(argc, argv, "tree-gen", "Generate a stochastic tree");
    auto N = parse_opt(parser, "--attr-points", "-N", "Number of attraction points", 1200);
    auto D = parse_opt(parser, "--distance", "-D", "Distance between two nodes", 1.0f);
    auto di = parse_opt(parser, "--influence-radius", "-di", "Radius of influence, equals <val> * D", 17);
    auto dk = parse_opt(parser, "--kill-distance", "-dk", "Kill distance, equals <val> * D", 2);
    auto iter_num = parse_opt(parser, "--iter-num", "-i", "Number of iterations", 100);
    auto path = parse_opt(parser, "--output", "-o", "Output file", "out/out.obj"s);

    if (parser._usage || should_exit(parser))
    {
        printf("%s\n", get_usage(parser).c_str());
        return 1;
    }

    auto scn = new scene();

    auto p0 = vec2f{2, 6};
    auto p1 = vec2f{10, 6};
    auto t0 = vec2f{11, 6};
    auto t1 = vec2f{11, 6};

    log_info("Generating attraction points...");

    auto vorodiag = throw_darts(N, p0, p1, t0, t1);

    log_info("Beginning tree generation...");

    vector<int> par;
    auto pos = grow(iter_num, N, D, di, dk, vorodiag, par);

    log_info("Drawing tree...");

    auto tree = draw_tree(pos, par);

    log_info("Saving model...");

    //Test scena
    tree->mat = add_test_material(scn, test_material_type::matte_colored);
    scn->shapes.push_back(tree);
    auto inst = new instance{"tree", identity_frame3f, tree};
    scn->instances.push_back(inst);
    add_test_camera(scn, test_camera_type::cam3);
    add_test_lights(scn, test_light_type::envlight);
    add_test_instance(scn, test_shape_type::floor, test_material_type::matte_green, identity_frame3f);
    add_test_environment(scn, test_environment_type::sky1, identity_frame3f);

    save_scene(path, scn, save_options());

    return 0;
}