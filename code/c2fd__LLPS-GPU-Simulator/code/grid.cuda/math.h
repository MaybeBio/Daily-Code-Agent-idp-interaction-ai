#pragma once
#include "grid.cuda/vec.h"
#include <cusp/blas.h>
namespace blas = cusp::blas;


//================
// math.inl

void laplacian(Vec &dest, const Vec &src, int i=0);
void divergent(Vec &dest, const Vec &velo);
void div_scalar_velo(Vec &dest, const Vec &a, const Vec3 &velo);
void div_scalar_grad(Vec &dest, const Vec &a, const Vec &b);
void v_dot_grad(Vec &dest, const Vec3 &velo, const Vec &x);


//================
// math.cu

// myreal sum(Vec &v);
myreal sum(Vec &v, int i=0);
myreal min(Vec &dest);
// myreal max(Vec &dest);
myreal max(Vec &dest, int i=0);
// myreal maxabs(Vec &dest);
myreal maxabs(Vec &dest, int i=0);

void clip(Vec &v, myreal min, myreal max);
void clip_and_print(Vec &x, myreal min, myreal max, const std::string &var_name, bool print_it);
void xpa(Vec &x, myreal a);

inline void project(const Vec &now, const Vec &old, myreal omega, Vec &dest) {
    // Need these myreal() casts.
    // Otherwise, cusp::axpby --> cusp::axpby<double>.
    // Since we don't use double in our CUDA card,
    // these double operations yield a bunch of zeros.
    blas::axpby(now, old, dest, omega, myreal(1) - omega);
}
