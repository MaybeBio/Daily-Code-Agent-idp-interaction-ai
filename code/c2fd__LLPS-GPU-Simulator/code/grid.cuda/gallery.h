#pragma once

namespace gallery {
    // void constant(Vec &dest, myreal a, int i=0);
    // void flat(Vec &v, myreal eps, myreal sol, int i=0);
    // void bump(Vec &v, myreal eps, myreal sol, bool is_3d, int i=0);
    // void bump_tanh(Vec &v, myreal eps, myreal sol, bool is_3d, int i=0);
    // void bump8(Vec &v, myreal eps, myreal sol, bool is_3d, int i=0);
    // void wall8(Vec &v, myreal eps, myreal sol, bool is_3d, int i=0);
    // void rod (Vec &v, myreal eps, myreal sol, bool is_3d, int i=0);
    // void nail(Vec &v, myreal eps, myreal sol, bool is_3d, bool gradient=false, int i=0);
    // void column(Vec &v, myreal eps, myreal sol, bool is_3d, int i=0);
    void my_sine(Vec &v, int i=0);
    void spatial_test_pattern(Vec &v, myreal dx, myreal dy, myreal dz, int i=0);
    void test_helmholtz(Vec &v, myreal dx, myreal dy, myreal dz, int toggle, int i=0);

    void set(myreal *v, const std::string &shape, myreal phi_max, myreal width);
    void set(Vec &v, const std::string &shape, myreal phi_max, myreal width=1e-6, int i=0);
}
