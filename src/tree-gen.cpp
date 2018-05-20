#include "yocto/yocto_gl.h"
#include "voro++/src/voro++.hh"

using namespace ygl;

const auto e = 2.71828184f;     // esponente per sommatoria dei raggi
const auto r0 = 0.01f;          // raggio iniziale
const auto leaf_threshold = r0 * 2;
const auto tree_leaf_ratio = 100;
const auto eps = 1e-4f;

// Semafori per la mutua esclusione
std::mutex mtx_dattr;
std::mutex mtx_comp;
std::mutex mtx_new;
std::mutex mtx_dnode;

void mkdir(const std::string& dir)
{
    #ifndef _MSC_VER
        auto e = system(("mkdir -p " + dir).c_str());
    #else
        auto fdir = dir;
        for (auto& c : fdir)
            if (c == '/') c = '\\';
        auto e = system(("mkdir " + fdir).c_str());
    #endif
}

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
        auto r = interpolate_bezier(p0, t0, t1, p1, (y - p0.x) / (p1.x - p0.x)).y;
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
    auto attr_id = 0, near_id = 0;
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

            auto skip = false;
            auto r = length(attr - node) + 1e-12;
            node_loop.setup_sphere(x, y, z, r, true);

            if (node_loop.start())
                do
                {
                    node_loop.pos(near_id, x, y, z, reder);
                    if (near_id != node_id)
                        skip = length(attr - vec3f{(float) x, (float) y, (float) z}) < r;
                }while(node_loop.inc());

            if(skip)
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

inline bool contains(const voro::container& voro, const vec3f& p)
{
    return p.x >= voro.ax && p.x <= voro.bx && p.y >= voro.ay && p.y <= voro.by && p.z >= voro.az && p.z <= voro.bz;
}

// Crescita dell'albero
std::vector<vec3f> grow(int iter_num, int N, float D, float di, float dk,
                   const voro::container& voro_attr, std::vector<int>& parents)
{
    auto positions = std::vector<vec3f>();
    auto positions_set = std::unordered_set<vec3f>();
    parents = std::vector<int>();

    positions.push_back({0, 0, 0});
    positions_set.insert({0, 0, 0});
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
    auto new_positions_set = std::unordered_set<vec3f>();

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
            auto par_node = new_nodes.back().second;

            if (contains(voro_nodes, new_node) && !positions_set.count(new_node))
            {
                voro_nodes.put(positions.size(), new_node.x, new_node.y, new_node.z);
                positions.push_back(new_node);
                new_positions_set.insert(new_node);
                parents.push_back(par_node);

                threads.push_back(std::thread(kill_points, new_node, dk, std::ref(voro_attr), std::ref(dead_attr)));
            }
            else
                dead_nodes.insert(par_node);

            new_nodes.pop_back();
        }

        positions_set.clear();
        positions_set = new_positions_set;
        new_positions_set.clear();

        for (auto& t : threads)
            t.join();
        threads.clear();

        if (i % 10 == 9 || !i)
            log_info("{} nodes added after {} iterations...", positions.size(), i + 1);
    }

    return positions;
}

shape* load_leaf(scene* scn, float scale)
{
    auto shp = new shape{"leaf"};

    shp->pos = std::vector<vec3f>{{-0.5, 0, 0}, {-0.5, 0, 3}, {0.5, 0, 3}, {0.5, 0, 0},
                                  {-0.5, 0, 0}, {-0.5, 0, 3}, {0.5, 0, 3}, {0.5, 0, 0}};

    for (auto i = 0; i < 8; i++)
    {
        shp->pos[i] = shp->pos[i] * scale;
        if (i > 3)
            shp->pos[i] = shp->pos[i] - vec3f{0, eps, 0};
    }

    shp->quads = std::vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4}};
    shp->norm = std::vector<vec3f>{{0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0},
                                   {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
    shp->texcoord = std::vector<vec2f>{{0, 1}, {1, 1}, {1, 0}, {0, 0},{0, 1}, {1, 1}, {1, 0}, {0, 0}};

    auto txt = new texture{"leaf", "leaf.png"};
    txt->ldr = load_image4b("resources/leaf.png");
    scn->textures.push_back(txt);
    auto mat = new material{"leaf", true};
    mat->kd = {1, 1, 1};
    mat->kd_txt = txt;
    shp->mat = mat;
    scn->materials.push_back(mat);

    return shp;
}

inline frame3f parallel_trans_frame(const vec3f &pos, const vec3f &tangent, const frame3f &pframe)
{
    auto b = cross(tangent, pframe.z);

    if (!length(b))
        return make_frame_fromzx(pos, tangent, pframe.x);

    b = normalize(b);
    auto t = acosf(dot(tangent, pframe.z));

    return make_frame_fromzx(pos, tangent, transform_point(rotation_frame(b, t), pframe.x));
}

inline void add_leaf(rng_pcg32& rng, scene* scn, const vec3f& pos, const vec3f& ppos, const vec3f& norm, const shape* shp)
{
    auto alpha = next_rand1f(rng);
    auto o = alpha * pos + (1 - alpha) * ppos;
    auto f = make_frame_fromz(o, norm);
    auto leaf = new shape{"leaf_" + std::to_string(scn->shapes.front()->shapes.size())};

    for (auto i = 0; i < 8; i++)
    {
        leaf->pos.push_back(transform_point(f, shp->pos[i]));
        leaf->norm.push_back(i < 4 ? f.y : -f.y);
    }
    leaf->quads = shp->quads;
    leaf->texcoord = shp->texcoord;
    leaf->mat = shp->mat;

    scn->shapes.front()->shapes.push_back(leaf);
}

// Crea la shape dell'albero con i quad dei cilindri
void draw_tree(scene* scn, float D, const std::vector<vec3f>& positions, const std::vector<int>& parents, const shape* leaf)
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

    auto group = new shape_group{"leaf", "", std::vector<shape*>()};
    scn->shapes.push_back(group);

    auto rad = std::vector<float>(positions.size(), 0.0f);
    auto children = std::vector<int>(positions.size(), 0);

    for (auto i = (int) positions.size() - 1; i >= 0; i--)
    {
        if (!rad[i])
            children[i] = -1;

        rad[i] = rad[i] ? pow(rad[i], 1 / e) : r0;
        rad[parents[i]] += pow(rad[i], e);

        children[parents[i]]++;
    }

    auto rng = init_rng(time(nullptr));
    auto ii = 0;

    for (auto i = (int) positions.size() - 1; i > 0; i--)
    {
        if (children[i] >= 0)   // Se il nodo non e' una foglia saltalo
            continue;

        auto child = i;
        auto pframe = frame3f();
        auto j = i;

        while (child)
        {
            auto pos = positions[j];
            auto tangent = normalize(positions[child] - positions[parents[j]]);
            auto f = children[j] < 0 ? parallel_trans_frame(pos, tangent, pframe)
                                  : make_frame_fromz(pos, tangent);

            for (auto jj = 0; jj <= 16; jj++)
            {
                auto u = (float) jj / 16;

                tree->pos.push_back(transform_point(f, {cosf(u * 2 * pif) * rad[j], sinf(u * 2 * pif) * rad[j], 0}));
                tree->norm.push_back(normalize(tree->pos.back() - pos));
                tree->texcoord.push_back({u, ii * D});

                if (children[j] >= 0 && jj)     // Se il nodo e' una foglia non aggiungere i quad
                {
                    tree->quads.push_back({ii * (16 + 1) + jj, (ii - 1) * (16 + 1) + jj,
                                           (ii - 1) * (16 + 1) + jj - 1, ii * (16 + 1) + jj - 1});

                    if (rad[j] < leaf_threshold && next_rand1f(rng) < D)
                        add_leaf(rng, scn, tree->pos.back(), tree->pos[(ii - 1) * (16 + 1) + jj], tree->norm.back(), leaf);
                }
            }
            ii++;
            if (children[j] > 0)    // Prosegui solo se e' l'ultima volta che visiti il nodo
                break;

            children[parents[j]]--;
            child = j;
            pframe = f;
            j = parents[j];
        }
    }

    group = new shape_group{"tree", "", std::vector<shape*>()};
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
    auto path = parse_opt(parser, "--output", "-o", "Output directory", "out"s);
    auto make_scene = parse_flag(parser, "--make-scene", "-s", "Add camera and environment");
    auto crown = parse_arg(parser, "crown shape", "Crown's shape", ""s, true, {"CONICAL", "CYLINDRICAL", "BEZIER"});

    vec2f p0, p1, t0, t1;

    if (crown == "CYLINDRICAL"s || crown == "CONICAL"s)
    {
        auto height = parse_arg(parser, "crown height", "Crown's height", 0.0f, true);
        auto radius = parse_arg(parser, "crown radius", "Crown's radius", 0.0f, true);
        auto trunk = parse_arg(parser, "trunk height", "Trunk's height", 0.0f, true);

        p0 = {trunk, radius};
        p1 = crown == "CYLINDRICAL"s
             ? vec2f{trunk + height, radius}
             : vec2f{trunk + height, 1};
        t0 = p0;
        t1 = p1;
    }
    else if (crown == "BEZIER"s)
    {
        p0.x = parse_arg(parser, "p0.x", "Bezier's spline 1st point", 0.0f, true);
        p0.y = parse_arg(parser, "p0.y", "Bezier's spline 1st point", 0.0f, true);
        t0.x = parse_arg(parser, "p1.x", "Bezier's spline 2nd point", 0.0f, true);
        t0.y = parse_arg(parser, "p1.y", "Bezier's spline 2nd point", 0.0f, true);
        t1.x = parse_arg(parser, "p2.x", "Bezier's spline 3rd point", 0.0f, true);
        t1.y = parse_arg(parser, "p2.y", "Bezier's spline 3rd point", 0.0f, true);
        p1.x = parse_arg(parser, "p3.x", "Bezier's spline 4th point", 0.0f, true);
        p1.y = parse_arg(parser, "p3.y", "Bezier's spline 4th point", 0.0f, true);
    }

    if (parser._usage || should_exit(parser))
    {
        printf("%s\n", get_usage(parser).c_str());
        return 1;
    }

    auto scn = new scene();

    log_info("Generating attraction points...");

    auto vorodiag = throw_darts(N, p0, p1, t0, t1);

    log_info("Beginning tree generation...");

    std::vector<int> par;
    auto pos = grow(iter_num, N, D, di, dk, vorodiag, par);

    log_info("Drawing tree...");

    auto leaf = load_leaf(scn, p1.x / tree_leaf_ratio);

    draw_tree(scn, D, pos, par, leaf);

    if (make_scene)
    {
        log_info("Creating scene...");

        auto inst = new instance{"tree", identity_frame3f, scn->shapes.back()};
        scn->instances.push_back(inst);
        inst = new instance{"leaf", identity_frame3f, scn->shapes.front()};
        scn->instances.push_back(inst);

        auto sky = new texture{"sky", "sky.hdr"};
        sky->hdr = make_sunsky_image(720, pif / 6);
        scn->textures.push_back(sky);

        auto env = new environment();
        env->name = "sky";
        env->ke = {1, 1, 1};
        env->ke_txt = scn->textures.back();
        scn->environments.push_back(env);

        scn->cameras.push_back(make_view_camera(scn, 0));
    }

    log_info("Saving model...");

    mkdir(path);
    save_scene(path + "/tree.obj", scn, save_options());

    return 0;
}