#pragma once
#include "tools.h"
#include <cusp/array1d.h>

// typedef cusp::device_memory MemorySpace;

typedef cusp::array1d<myreal, cusp::device_memory> Vec;
typedef cusp::array1d<myreal, cusp::device_memory> Vec3;
typedef cusp::array1d<myreal, cusp::host_memory> H_Vec;
typedef cusp::array1d<myreal, cusp::host_memory> H_Vec3;


//================
// Unwrap pointer
//
// inline myreal *cast(Vec &v) {return thrust::raw_pointer_cast(&(v)[0]);}
// inline myreal *cast(H_Vec &v) {return thrust::raw_pointer_cast(&(v)[0]);}
// inline const myreal *cast(const Vec &v)   {return thrust::raw_pointer_cast(&(v)[0]);}
// inline const myreal *cast(const H_Vec &v) {return thrust::raw_pointer_cast(&(v)[0]);}

template <typename T>
inline myreal *cast(T &v) {
    return thrust::raw_pointer_cast(&(v)[0]);
}
template <typename T>
inline const myreal *cast(const T &v) {
    return thrust::raw_pointer_cast(&(v)[0]);
}


//================
// Wrap pointer
//
//   device pointer --> thrust::device_ptr --> cusp::array1d

// thrust::device_ptr<int> dev_ptr(raw_ptr);
// thrust::device_pointer_cast(raw_ptr)   // type inferred by template


template <typename T>
cusp::array1d_view<T> wrap(T *raw_ptr, size_t size) {
    thrust::device_ptr<T> dev_ptr(raw_ptr);    
    return cusp::array1d_view<T>(dev_ptr, dev_ptr + size);
    // same thing
    // return cusp::make_array1d_view(dev_ptr, dev_ptr + size);
}
