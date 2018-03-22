#include "yocto/yocto_gl.h"
#include "yocto/yocto_gl.cpp"
#include "voro++-0.4.6/src/voro++.hh"

using namespace ygl;

/* Genera in modo randomico i punti di attrazione e crea il
   diagramma di Voronoi */
voro::container throw_darts(int N, vec2f p0, vec2f p1, vec2f t0, vec2f t1)
{
    auto points = std::vector<vec3f>();
    auto bbox = bbox3f();
    const auto eps = 0.02f;
    auto rng = init_rng(time(nullptr));

    for (auto i = 0; i < N; i++)
    {
        auto y = next_rand1f(rng, p0.x, p1.x);
        auto r = next_rand1f(rng, 0, eval_bezier_cubic(p1, t1, t0, p0, y).y);
        auto t = next_rand1f(rng, 0, 2 * pif);

        points.push_back({cosf(t) * r, y, sinf(t) * r});

        bbox += points.back();
    }

    auto vorodiag = voro::container(bbox.min.x - eps, bbox.max.x + eps, bbox.min.y - eps, bbox.max.y + eps,
                                    bbox.min.z - eps, bbox.max.z + eps, N / 5, N / 5, N / 5, false, false, false, 8);

    // aggiungo tutti i punti di attrazione
    for (auto i = 0; i < N; i++)
        vorodiag.put(i, points[i].x, points[i].y, points[i].z);

    return vorodiag;
}

// Crescita dell'albero
voro::container grow(int iter_num, int N, float D, float dk, float di, const voro::container& voro_attr,
                     vector<vec3f>& positions, vector<int>& parents)
{
    auto nodes_id = 0;

    positions.push_back({0, 0, 0});
    parents.push_back(0);

    auto attr_loop = voro::c_loop_subset(voro_attr);

    auto voro_nodes = voro::container(voro_attr.ax, voro_attr.bx, min(voro_attr.ay, -0.02f), voro_attr.by,
                                      voro_attr.az, voro_attr.bz, N / 5, N / 5, N / 5, false, false, false, 8);
    voro_nodes.put(nodes_id++, 0.0f, 0.0f, 0.0f);
    auto nodes_loop = voro::c_loop_all(voro_nodes);

    auto dead_attr = std::unordered_set<int>();
    auto computed_attr = std::unordered_set<int>();
    auto dead_nodes = std::unordered_set<int>();
    auto new_nodes = std::vector<vec3f>();

    auto node_id = 0, attr_id = 0, search_id = 0;
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

                        voro_nodes.find_voronoi_cell(attr.x, attr.y, attr.z, x, y, z, search_id);

                        if (search_id != node_id)
                            continue;

                        computed_attr.insert(attr_id);

                        sum += normalize(attr - node);
                        summed = true;
                    } while (attr_loop.inc());

                if (summed)
                {
                    new_nodes.insert(new_nodes.begin(), node + D * normalize(sum));
                    parents.push_back(node_id);
                }
                else
                    dead_nodes.insert(node_id);
            } while (nodes_loop.inc());

        while (!new_nodes.empty())
        {
            auto new_node = new_nodes.back();
            voro_nodes.put(nodes_id++, new_node.x, new_node.y, new_node.z);
            positions.push_back(new_node);

            attr_loop.setup_sphere(new_node.x, new_node.y, new_node.z, dk, true);

            if (attr_loop.start())
                do
                {
                    attr_loop.pos(attr_id, x, y, z, reder);
                    dead_attr.insert(attr_id);
                } while (attr_loop.inc());

            new_nodes.pop_back();
        }
    }

    return voro_nodes;
}

void make_cylinder(shape* tree, vec3f node, vec3f p_node, float r)
{
    auto points = (int) tree->pos.size();

    for (auto i = 0; i <= 16; i++)
        for (auto j = 0; j <= 16; j++)
        {
            auto u = 2 * pif * i / 16;
            auto v = (float) j / 16;

            auto p = v * p_node + (1 - v) * node;

            tree->pos.push_back({cosf(u) * r + p.x, p.y, sinf(u) * r + p.z});
            tree->norm.push_back(normalize(p - tree->pos.back()));
            tree->texcoord.push_back({u / (2 * pif), v});

            if (i != 16 && j != 16)
                tree->quads.push_back(
                        {points + i * (16 + 1) + j, points + (i + 1) * (16 + 1) + j, points + (i + 1) * (16 + 1) + j + 1,
                         points + i * (16 + 1) + j + 1});
        }
}

// Crea la shape dell'albero con i quad dei cilindri
shape* draw_tree(const voro::container& voro_nodes, const std::vector<vec3f> positions, const std::vector<int>& parents)
{
    auto shp = new shape{"tree"};

    for (auto i = (int) positions.size() - 1; i > 0; i--)
    {
        auto pos = positions[i];

        auto par = parents[i];
        auto ppos = positions[par];

        make_cylinder(shp, pos, ppos, 0.01f);
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

    auto p0 = vec2f{0, 0};
    auto p1 = vec2f{10, 0};
    auto t0 = vec2f{6, 2};
    auto t1 = vec2f{1, 4};

    auto voro_attr = throw_darts(N, p0, p1, t0, t1);

    auto pos = vector<vec3f>();
    auto par = vector<int>();

    auto voro_nodes = grow(iter_num, N, D, dk, di, voro_attr, pos, par);

    auto tree = draw_tree(voro_nodes, pos, par);

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