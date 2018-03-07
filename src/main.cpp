#include "../src/yocto/yocto_gl.h"
#include "../src/yocto/yocto_gl.cpp"

using namespace ygl;



int main()
{
    auto p0 = vec2f{0, 0};
    auto p1 = vec2f{1, 0};
    auto t0 = vec2f{0.6, 0.2};
    auto t1 = vec2f{0.1, 0.4};

    // vec2f x = eval_bezier_cubic(p1, t1, p0, t0, 0); py confronto

    for (auto i = 0; i < 100; i++)
    {
        auto x = i / 100.0f;
        vec2f p = eval_bezier_cubic(p1, t1, t0, p0, x);
        for (auto j = 0; j < (int) (p.y * 100); j++)
        {
            printf(".");
        }
        printf("\n");
    }


    return 0;
}