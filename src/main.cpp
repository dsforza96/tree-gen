#include "yocto/yocto_gl.h"
#include "yocto/yocto_gl.cpp"
#include "voro++-0.4.6/src/voro++.hh"

using namespace ygl;

// Genera in modo randomico i punti di attrazione
voro::container throw_darts(int N, vec2f p0, vec2f p1, vec2f t0, vec2f t1)
{
    auto points = std::vector<vec3f>();
    auto bbox = bbox3f();
    const auto eps = 0.02f;

    auto rng = init_rng(time(nullptr));

    for (auto i = 0; i < N; i++)
    {
        auto y = next_rand1f(rng, p0.x, p1.x);
        auto p = eval_bezier_cubic(p1, t1, t0, p0, y);
        auto xz = next_rand2f(rng, -p, p);

        points.push_back({xz.x, y, xz.y});

        bbox += {xz.x, y, xz.y};
    }

    auto vorodiag = voro::container(bbox.min.x - eps, bbox.max.x + eps, bbox.min.y - eps, bbox.max.y + eps,
                                    bbox.min.z - eps, bbox.max.z + eps, N / 5, N / 5, N / 5, false, false, false, 8);

    for (auto i = 0; i < N; i++)
        vorodiag.put(i, points[i].x, points[i].y, points[i].z);

    return vorodiag;
}

void grow();

int main()
{
    auto N = 10;
    auto iter_num = 1;
    auto dist = 10.0f;
    auto p0 = vec2f{0, 0};
    auto p1 = vec2f{1, 0};
    auto t0 = vec2f{0.6, 0.2};
    auto t1 = vec2f{0.1, 0.4};
    auto root = vec3f{0, 0, 0};

    // vec2f x = eval_bezier_cubic(p1, t1, p0, t0, 0); py confronto

    auto vorodiag = throw_darts(N, p0, p1, t0, t1);
    auto voroloop = voro::c_loop_subset(vorodiag);

    for (auto i = 0; i < iter_num; i++)
    {
        //auto pid = 0;
        auto x = 0.0, y = 0.0, z = 0.0;
        //vorodiag->find_voronoi_cell(root.x, root.y, root.z, x, y, z, pid);

        voroloop.setup_sphere(root.x, root.y, root.z, dist, true);

        while (voroloop.inc())
        {
            voroloop.pos(x, y, z);
            printf("%d x: %f, y: %f, z: %f\n", ++i, x, y, z);
        }
    }

    return 0;
}