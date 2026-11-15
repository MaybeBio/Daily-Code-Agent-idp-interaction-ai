#include "grid.cuda/math.h"

#include "grid.cuda/grid_params.h"
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <string>


// Important:
//    The functions pad_reflex, pad_velo, ...
//    only pad the sides, but edges and corners of the padding
//    still contain junk data.  Don't let it affect eta_ave.
//    We only sum/max/min over data[1:-1, : , : ].


// myreal sum(Vec &v) {
//     return thrust::reduce(v.begin(), v.end());
// }
myreal sum(Vec &v, int i) {
    return thrust::reduce(v.begin() + i*gp.Ndata + gp.sy, v.begin() + (i+1)*gp.Ndata - gp.sy);
}


myreal min(Vec &v) {
    return thrust::reduce( v.begin() + gp.sy, v.end() - gp.sy, myreal(0.), thrust::minimum<myreal>());
}


// myreal max(Vec &v) {
//     return thrust::reduce(
// 	v.begin(), v.end(),
// 	myreal(0.), thrust::maximum<myreal>());
// }
myreal max(Vec &v, int i) {
    return thrust::reduce( v.begin() + i*gp.Ndata + gp.sy, v.begin() + (i+1)*gp.Ndata - gp.sy, myreal(0.), thrust::maximum<myreal>());
}


// myreal maxabs(Vec &v) {
//     return thrust::transform_reduce(
// 	v.begin(), v.end(),
// 	thrust::absolute_value<myreal>(),
// 	myreal(0.), thrust::maximum<myreal>());
// }

template<typename T>
struct absolute_value : public std::unary_function<T,T> {
  __host__ __device__ T operator()(const T &x) const {return x < T(0) ? -x : x;}
};


myreal maxabs(Vec &v, int i) {
    // 2011-07-24 -- Thrust bug?  This returns 0. if there is any NaN entry.
    return thrust::transform_reduce(
	v.begin() + i*gp.Ndata + gp.sy,
	v.begin() + (i+1)*gp.Ndata - gp.sy,
	absolute_value<myreal>(),
	myreal(0.), thrust::maximum<myreal>());
}








struct Clip {
    myreal a, b;
    Clip(myreal a, myreal b) : a(a), b(b) {assert(a<=b);}
    __host__ __device__
    void operator()(myreal &t) {t = min(max(t,a),b);}
};

void clip(Vec &x, myreal min, myreal max) {
    thrust::for_each(x.begin(), x.end(), Clip(min, max));
}

void clip_and_print(Vec &x, myreal min, myreal max, const std::string &var_name, bool print_it) {
    if (print_it) {
        Vec old(x);
        clip(x,min,max);
        // old = old-x
        blas::axpy(x,old,myreal(-1));
        std::cout << var_name << " clipped.  sum() increases by approximately " << -sum(old)/gp.NypNzp << std::endl;
    } else {
        clip(x,min,max);
    }
}



struct XPA {
    myreal a;
    XPA(myreal a) : a(a) {}
    __host__ __device__
    void operator()(myreal &t) {t += a;}
};

void xpa(Vec &x, myreal a) {
    thrust::for_each(x.begin(), x.end(), XPA(a));
}
