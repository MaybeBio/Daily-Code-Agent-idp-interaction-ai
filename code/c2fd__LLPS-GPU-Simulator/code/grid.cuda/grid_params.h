#pragma once
#include "tools.h"
#include <cassert>
#include <iostream>


class GridParams;
extern GridParams gp;


const int Ndim = 2;
const int NdimTau = Ndim*(Ndim+1)/2; // 6
const int Npad = 1;



inline int ceil_div(int a, int b) {return (a+b-1)/b;}

class GridParamsBase {
public:
    // The program assumes that these variables appears x/y/z consecutively in memory.
    int Ny , Nz;
    int Ny1, Nz1;
    int Nyp, Nzp;
    bool is_3d;

    myreal Ly, Lz;
    myreal dy, dz;
    myreal _dy, _dz ;	// 1/dx
    myreal _dy2, _dz2;	// 1/dx/dx



    // set by subclass
    // dimBlock ... in C ordering
    int dby;
    int dbz;
    int BLOCK_SIZE_3D;

    int Nypp, Nzpp;
    int sy, sz;		// stride, sz == 1
    size_t Ndata;
    size_t NypNzp;

    int wrap_y, wrap_z;
    int periodicity;		// 0=wall, 1=periodic, xyz at bit 012
    //int phi_floor_pattern_num;	// For gallery::pad_phi()

    GridParamsBase() {}		// For device side
    GridParamsBase(int Ny_, int Nz_, myreal Ly_, myreal Lz_) {
	Ny = Ny_; Nz = Nz_;
	Ly = Ly_; Lz = Lz_;

	Ny1 = Ny+Npad;
	Nz1 = Nz+Npad;
	Nyp = Ny+2*Npad;
	Nzp = Nz+2*Npad;

	is_3d = false; //(Nx > 2);

	dy = Ly/Ny;
	dz = Lz/Nz;
	_dy = 1/dy;
	_dz = 1/dz;
	_dy2 = _dy*_dy;
	_dz2 = _dz*_dz;
    }

    void set_periodicity(bool y, bool z) {
	// true:  periodic
	// false: wall
	assert(Ndata);

	periodicity = 0;
	if (y)  periodicity += 1<<0;
	if (z)  periodicity += 1<<1;

	wrap_y = (y ? Ny : 1) * sy;
	wrap_z = (z ? Nz : 1) * sz;
    }
    inline bool is_periodic_y() {return periodicity & 1<<0;}
    inline bool is_periodic_z() {return periodicity & 1<<1;}
    inline int ind(int j, int k)	const {return j*sy + k;}
};



class GridParams : public GridParamsBase {
public:
    dim3 dimBlock;
    dim3 dimGrid;
    int BLOCK_SIZE_1D;
    int BLOCK_SIZE_1D_MULT;

    GridParams() {}
    GridParams(int Ny, int Nz,
               int dby_, int dbz_,
	       myreal Ly, myreal Lz)
	: GridParamsBase(Ny, Nz, Ly, Lz) {

        dby = dby_;
        dbz = dbz_;
        BLOCK_SIZE_3D = dby*dbz;

	// dimGrid ... in C ordering
	int dgy = ceil_div(Nyp, dby);
	int dgz = ceil_div(Nzp, dbz);

	dimBlock = dim3(dbz,dby); // fortran
	dimGrid  = dim3(dgz,dgy); // fortran; sic *

	// GridParamsBase
	Nypp = dgy*dby;
	Nzpp = dgz*dbz;
	sz = 1;
	sy = Nzpp;
	Ndata = Nypp * sy;
	NypNzp = Nyp * Nzp;


	// Ensure that pad cell and the element that use it
	// get called in the same block.
	// One way to satisfy this is:
	//    gcd(Nxp,dbx)>1 && gcd(Nyp,dby)>1 && gcd(Nyz,dbz)>1

	assert(dby>=2*Npad && dbz>=2*Npad);
	assert((Nyp-1)/dby == (Ny1-Npad)/dby);
	assert((Nzp-1)/dbz == (Nz1-Npad)/dbz);

	// For device_interface.inl: laplacian, etc.
	BLOCK_SIZE_1D = 128;
	BLOCK_SIZE_1D_MULT = 32;
	assert(BLOCK_SIZE_1D <= 512);  // for compute 1.3
	assert(BLOCK_SIZE_1D <= 1024); // for compute 2.0
	assert(dimGrid.x<=65535 && dimGrid.y<=65535); // CUDA's spec
	assert(dimBlock.z<=64);
	if ((Ny+3)*(Nz+3) >= (1<<24)) {
	    // The +3 above are just for safety.
	    std::cout << "\n\n\n"
		      << "\n    Nx*Ny*Nz exceeds CUFFT's spec version 3.2"
		      << "\n    Want Nx*Ny*Nz < 256^3 == 2^24."
		      << "\n\n\n" << std::endl;
	    // http://developer.download.nvidia.com/compute/cuda/3_2_prod/toolkit/docs/CUDA_Toolkit_Release_Notes_Linux.txt
	    //
	    // CUFFT supports batched transforms > 512 MB. The previous
	    // version of CUFFT failed when (batch size * transform size *
	    // datatype size) for a 1D single-precision transform exceeded
	    // 512MB. This has been fixed so that now the total size can
	    // be as large as the device memory capacity allows. The exact
	    // size varies depending on whether operating in-place or
	    // out-of-place and depending on how much internal
	    // intermediate memory is required by API, which can vary
	    // depending on the actual size of the transform. Note
	    // however, that while the total size can be much larger, the
	    // size of each individual transform is still limited to 128 M
	    // elements (1 GB for single-precision.)
	}


	// Need this for __shared__ a[3*BLOCK_SIZE_3D]
	//assert(BLOCK_SIZE_3D == dimBlock.x * dimBlock.y*dimBlock.z);
	assert(BLOCK_SIZE_3D == dimBlock.x * dimBlock.y);


	// // From device.inl::let_idx_(ijk)
	// assert(dbx*dgx < (1<<10));
	// assert(dby*dgy < (1<<10));
	// assert(dbz*dgz < (1<<11));
    }

    void copy_to_device() const;
};

std::ostream& operator<<(std::ostream& os, const GridParams &gp);
