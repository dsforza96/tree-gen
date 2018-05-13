#include "yocto/yocto_gl.h"
#include "yocto/yocto_gl.cpp"
#include "voro++/src/voro++.hh"

using namespace ygl;

const auto e = 2.71828184f;     // esponente per sommatoria dei raggi
const auto r0 = 0.01f;          // raggio iniziale
const auto leaf_threshold = 0.02f;
const auto leaf_prob = 0.2f;
const auto eps = 1e-16f;

// Semafori per la mutua esclusione
std::mutex mtx_dattr;
std::mutex mtx_comp;
std::mutex mtx_new;
std::mutex mtx_dnode;

/* Genera in modo randomico i punti di attrazione e crea il
   diagramma di Voronoi */
voro::container throw_darts(int N, const vec2f& p0, const vec2f& p1, const vec2f& t0, const vec2f& t1)
{
    auto points = std::vector<vec3f>();
    auto bbox = bbox3f();
    auto rng = init_rng(time(nullptr));

    for (auto i = 0; i < N; i++)
    {
        auto y = next_rand1f(rng, p0.x, p1.x);
        auto r = interpolate_bezier(p1, t1, t0, p0, y).y;
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
                const voro::container& voro_nodes, std::unordered_set<int>& dead_attr,
                std::unordered_set<int>& computed_attr, std::unordered_set<int>& dead_nodes,
                std::vector<std::pair<vec3f, int>>& new_nodes)
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

            mtx_comp.lock();
            computed_attr.insert(attr_id);
            mtx_comp.unlock();

            sum += normalize(attr - node);
            summed = true;
        } while (attr_loop.inc());

    if (summed)
    {
        auto new_node = node + D * normalize(sum);
        mtx_new.lock();
        new_nodes.push_back({new_node, node_id});
        mtx_new.unlock();
    }
    else
    {
        mtx_dnode.lock();
        dead_nodes.insert(node_id);
        mtx_dnode.unlock();
    }
}

void kill_points(vec3f node, float  dk, const voro::container& voro_attr, std::unordered_set<int>& dead_attr)
{
    auto attr_loop = voro::c_loop_subset(voro_attr);
    attr_loop.setup_sphere(node.x, node.y, node.z, dk, true);

    auto attr_id = 0;
    auto x = 0.0, y = 0.0, z = 0.0, reder = 0.0;

    if (attr_loop.start())
        do
        {
            attr_loop.pos(attr_id, x, y, z, reder);

            if (dead_attr.count(attr_id))
                continue;

            mtx_dattr.lock();
            dead_attr.insert(attr_id);
            mtx_dattr.unlock();
        } while (attr_loop.inc());
}

// Crescita dell'albero
std::vector<vec3f> grow(int iter_num, int N, float D, float di, float dk,
                   const voro::container& voro_attr, std::vector<int>& parents)
{
    auto positions = std::vector<vec3f>();
    parents = std::vector<int>();

    positions.push_back({0, 0, 0});
    parents.push_back(0);

    auto attr_loop = voro::c_loop_subset(voro_attr);

    auto voro_nodes = voro::container(voro_attr.ax, voro_attr.bx, min(voro_attr.ay, (double) -eps), voro_attr.by,
                                      voro_attr.az, voro_attr.bz, N / 5, N / 5, N / 5, false, false, false, 8);
    voro_nodes.put(0, 0.0f, 0.0f, 0.0f);
    auto nodes_loop = voro::c_loop_all(voro_nodes);

    auto dead_attr = std::unordered_set<int>();
    auto computed_attr = std::unordered_set<int>();
    auto dead_nodes = std::unordered_set<int>();
    auto new_nodes = std::vector<std::pair<vec3f, int>>();

    auto threads = std::vector<std::thread>();

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
                                              std::ref(voro_attr), std::ref(voro_nodes), std::ref(dead_attr),
                                              std::ref(computed_attr), std::ref(dead_nodes), std::ref(new_nodes)));
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

            threads.push_back(std::thread(kill_points, new_node, dk, std::ref(voro_attr), std::ref(dead_attr)));

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

void load_leaf(scene* scn, const std::string& name)
{
    auto shp = new shape{"leaf"};

    shp->pos = std::vector<vec3f>{{-0.1, 0, 0}, {-0.1, 0, 0.6}, {0.1, 0, 0.6}, {0.1, 0, 0},
                                  {-0.1, -1e-4, 0}, {-0.1, -1e-4, 0.6}, {0.1, -1e-4, 0.6}, {0.1, -1e-4, 0}};
    shp->quads = std::vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7}};
    shp->norm = std::vector<vec3f>{{0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0},
                                   {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
    shp->texcoord = std::vector<vec2f>{{0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}};

    auto txt = new texture{"leaf", "leaf.png"};
    txt->ldr = load_image4b(name);
    scn->textures.push_back(txt);
    auto mat = new material{"leaf", true};
    mat->kd = {1, 1, 1};
    mat->kd_txt = txt;
    shp->mat = mat;
    scn->materials.push_back(mat);

    auto group = new shape_group{"leaf", "", std::vector<shape*>()};
    group->shapes.push_back(shp);
    scn->shapes.push_back(group);
}

inline frame3f compute_frame(const vec3f& pos, const vec3f& tangent, const frame3f& pframe)
{
    auto b = cross(tangent, pframe.z);

    if (!length(b))
        return make_frame_fromzx(pos, tangent, pframe.x);

    b = normalize(b);
    auto t = acosf(dot(tangent, pframe.z));

    return make_frame_fromzx(pos, tangent, transform_point(rotation_frame(b, t), pframe.x));
}

inline void add_leaf(rng_pcg32& rng, scene* scn, const vec3f& pos, const vec3f& ppos, const vec3f& norm)
{
    auto alpha = next_rand1f(rng);
    auto o = alpha * pos + (1 - alpha) * ppos;
    auto inst = new instance{"leaf", make_frame_fromz(o, norm), scn->shapes.front()};
    scn->instances.push_back(inst);
}

// Crea la shape dell'albero con i quad dei cilindri
void draw_tree(scene* scn, float D, const std::vector<vec3f> positions, const std::vector<int>& parents)
{
    auto tree = new shape{"tree"};

    auto txt = new texture{"bark", "bark.png"};
    txt->ldr = load_image4b("resources/bark.png");
    scn->textures.push_back(txt);
    auto mat = new material{"bark", true};
    mat->kd = {1, 1, 1};
    mat->kd_txt = txt;
    tree->mat = mat;
    scn->materials.push_back(mat);

    auto rad = std::vector<float>(positions.size(), 0.0f);
    auto children = std::vector<int>(positions.size(), 0);

    for (auto i = (int) positions.size() - 1; i >= 0; i--)
    {
        rad[i] = rad[i] ? pow(rad[i], 1 / e) : r0;
        rad[parents[i]] += pow(rad[i], e);

        children[parents[i]]++;
    }

    auto rng = init_rng(time(nullptr));
    auto ii = 0;

    for (auto i = (int) positions.size() - 1; i > 0; i--)
    {
        if (rad[i] != r0)
            continue;

        auto child = i;
        auto pframe = frame3f();
        auto j = i;

        while (child)
        {
            auto pos = positions[j];
            auto tangent = normalize(positions[child] - positions[parents[j]]);
            auto f = rad[j] != r0 ? compute_frame(pos, tangent, pframe)
                                  : make_frame_fromz(pos, tangent);

            for (auto jj = 0; jj <= 16; jj++)
            {
                auto u = (float) jj / 16;

                tree->pos.push_back(transform_point(f, {cosf(u * 2 * pif) * rad[j], sinf(u * 2 * pif) * rad[j], 0}));
                tree->norm.push_back(normalize(tree->pos.back() - pos));
                tree->texcoord.push_back({u, ii * D});

                if (rad[j] != r0 && jj)
                {
                    tree->quads.push_back({ii * (16 + 1) + jj, (ii - 1) * (16 + 1) + jj,
                                           (ii - 1) * (16 + 1) + jj - 1, ii * (16 + 1) + jj - 1});

                    if (rad[j] < leaf_threshold && next_rand1f(rng) < leaf_prob)
                        add_leaf(rng, scn, tree->pos.back(), tree->pos[(ii - 1) * (16 + 1) + jj], tree->norm.back());
                }
            }
            ii++;
            if (children[j])
                break;

            children[parents[j]]--;
            child = j;
            pframe = f;
            j = parents[j];
        }
    }

    auto group = new shape_group{"tree", "", std::vector<shape*>()};
    group->shapes.push_back(tree);
    scn->shapes.push_back(group);
}

int main(int argc, char** argv)
{
    auto parser = make_parser(argc, argv, "tree-gen", "Generate a stochastic tree");
    auto N = parse_opt(parser, "--attr-points", "-N", "Number of attraction points", 1200);
    auto D = parse_opt(parser, "--distance", "-D", "Distance between two nodes", 1.0f);
    auto di = parse_opt(parser, "--influence-radius", "-di", "Radius of influence, equals <val> * D", 17) * D;
    auto dk = parse_opt(parser, "--kill-distance", "-dk", "Kill distance, equals <val> * D", 2) * D;
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

    std::vector<int> par;
    auto pos = grow(iter_num, N, D, di, dk, vorodiag, par);

    log_info("Drawing tree...");

    load_leaf(scn, "resources/leaf.png");

    draw_tree(scn, D, pos, par);

    log_info("Saving model...");

    auto inst = new instance{"tree", identity_frame3f, scn->shapes.back()};
    scn->instances.push_back(inst);

    save_scene(path, scn, save_options());

    return 0;
}