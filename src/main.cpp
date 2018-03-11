#include "yocto/yocto_gl.h"
#include "yocto/yocto_gl.cpp"
#include "cppdelaunay/delaunay/Voronoi.h"
#include "cppdelaunay/delaunay/Voronoi.cpp"

using namespace ygl;

// Genera in modo randomico i punti di attrazione
shape* throw_darts(int N, vec2f p0, vec2f p1, vec2f t0, vec2f t1)
{
    auto shp = new shape{"points"};

    auto rng = init_rng(time(nullptr));

    for (auto i = 0; i < N; i++)
    {
        auto y = next_rand1f(rng, p0.x, p1.x);
        auto p = eval_bezier_cubic(p1, t1, t0, p0, y);
        auto xz = next_rand2f(rng, -p, p);

        shp->pos.push_back({xz.x, y, xz.y)});
    }

    return shp;
}

int main()
{
    auto N = 10;
    auto p0 = vec2f{0, 0};
    auto p1 = vec2f{1, 0};
    auto t0 = vec2f{0.6, 0.2};
    auto t1 = vec2f{0.1, 0.4};
    auto root = vec3f{0, 0, 0};

    // vec2f x = eval_bezier_cubic(p1, t1, p0, t0, 0); py confronto

    auto points = throw_darts(N, p0, p1, t0, t1);



    return 0;
}