#pragma once
#include "grid.cuda/grid.h"

// #include <cusp/blas.h>
// namespace blas = cusp::blas;


// template <class T>
// class TwoOf {
// public:
//     T now, old;		// Can't use "new" since it is a reserved word.

//     TwoOf<T>(){}
//     template<typename ValueType>
//     TwoOf<T>(ValueType v) : now(v), old(v) {}
//     template<typename ValueType, typename ValueType2>
//     TwoOf<T>(ValueType v, ValueType2 w) : now(v,w), old(v,w) {}

//     // T& operator=(T a) {old = now; now = a; return now;}
//     operator T&() {return now;}
//     operator const T&() const {return now;}

//     // (now, old) = (step n,n-1).
//     // Project to step n+alpha.
//     void project(T &dest, myreal omega = 2.) const;
// };

// template<>
// inline void TwoOf<Vec>::project(Vec &dest, myreal omega) const {
//     // Need these myreal() casts.
//     // Otherwise, cusp::axpby --> cusp::axpby<double>.
//     // Since we don't use double in our CUDA card,
//     // these double operations yield a bunch of zeros.
//     blas::axpby(now, old, dest, omega, myreal(1) - omega);
// }
