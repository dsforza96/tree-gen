#include "yocto/yocto_gl.h"
#include "yocto/yocto_gl.cpp"
#include "voro++-0.4.6/src/voro++.hh"

using namespace ygl;

// Genera in modo randomico i punti di attrazione
voro::container* throw_darts(int N, vec2f p0, vec2f p1, vec2f t0, vec2f t1)
{
    auto points = std::vector<vec3f>();
    auto xmax = 0.0f;
    auto xmin = 0.0f;
    auto zmax = 0.0f;
    auto zmin = 0.0f;

    auto rng = init_rng(time(nullptr));


    for (auto i = 0; i < N; i++)
    {
        auto y = next_rand1f(rng, p0.x, p1.x);
        auto p = eval_bezier_cubic(p1, t1, t0, p0, y);
        auto xz = next_rand2f(rng, -p, p);

        points.push_back({xz.x, y, xz.y});

        xmax = max(xmax, p.x);
        xmin = min(xmin, p.x);
        zmax = max(zmax, p.y);
        zmin = min(zmin, p.y);
    }

    auto vorodiag = new voro::container(xmin, xmax, p0.x, p1.x, zmin, zmax, N / 5, N / 5, N / 5, false, false, false, 8);

    for (auto i = 0; i < N; i++)
        vorodiag->put(i, points[i].x, points[i].y, points[i].z);

    return vorodiag;
}

void grow();

int main()
{
    auto N = 10;
    auto iter_num = 100;
    auto p0 = vec2f{0, 0};
    auto p1 = vec2f{1, 0};
    auto t0 = vec2f{0.6, 0.2};
    auto t1 = vec2f{0.1, 0.4};
    auto root = vec3f{0, 0, 0};

    // vec2f x = eval_bezier_cubic(p1, t1, p0, t0, 0); py confronto

    auto vorodiag = throw_darts(N, p0, p1, t0, t1);

    for (auto i = 0; i < iter_num; i++)
    {
        auto pid = 0;
        auto x = 0.0d, y = 0.0d, z = 0.0d;
        vorodiag->find_voronoi_cell(root.x, root.y, root.z, x, y, z, pid);
    }

    return 0;
}