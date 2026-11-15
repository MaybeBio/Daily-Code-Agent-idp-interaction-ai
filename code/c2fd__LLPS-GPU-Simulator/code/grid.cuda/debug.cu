#include "grid.cuda/debug.h"
#include "grid.cuda/grid_params.h"
#include <iostream>
#include <iomanip>
#include <cuda_runtime.h>
#include <cstdio>
using namespace std;

GridException::GridException(const string text) {
	cout << "GridException: " << text << endl;
	print_backtrace();
}

void print_memory_usage() {
	size_t free, total;
	cudaMemGetInfo(&free, &total);
	cout << "Momory used+free (MB): "
		<< fixed << setprecision(2)
		   << setw(7) << (total-free)/1e6 
		   << " + " << setw(7) << free/1e6
		   << " = " << setw(7) << total/1e6 << endl;
}


void check(string label) {
	bool need_newline = false;
	if ((DEBUG & DEBUG_VERBOSE_STEPS && label!="") ||
			(DEBUG & DEBUG_VERBOSE_MEMORY)) {
		cout << left << setw(26) << label;
		need_newline = true;
	}

	if (DEBUG & DEBUG_CHECK) {
		cudaThreadSynchronize();
		cudaError_t e = cudaGetLastError();
		if (e!=cudaSuccess) {
			cout<<"CudaErrorString = " << cudaGetErrorString(e)<<endl;
			throw GridException("check");
		}
	}

	if (DEBUG & DEBUG_VERBOSE_MEMORY) {
		print_memory_usage();
	} else if (need_newline) {
		cout << endl;
	}
}



void check(const string &label, const Vec &v, vector<int> pads) {
	check(label);
	list_nan_inf(v, pads, true);
}

bool list_nan_inf_(const myreal *v, 
		int d,       // only for printing
		int px, int py, int pz,
		int qx, int qy, int qz) {
	// px == pad size at x = 0
	// qx == pad size at x = Lx
	bool found = false;

	for (int j=Npad-py; j<gp.Ny1+qy; j++)
		for (int k=Npad-pz; k<gp.Nz1+qz; k++) {
			int idx = j*gp.sy + k*gp.sz;
			myreal val = v[idx];
			if (!isfinite(val)) {
				// printf automatically flushes at newline.
				printf("(%d,%3d,%3d) = %g\n", d,j,k,val);
				found = true;
			}
		}

	return found;
}

void list_nan_inf(const Vec &v, vector<int> pads, bool quit) {
	if (!(DEBUG & DEBUG_LIST_NAN_INF))  return;

	bool found = false;
	H_Vec hv(v);

	if (DEBUG & DEBUG_LIST_NAN_INF_VERBOSE)
		cout<<"list_nan_inf()"<<endl;

	assert(v.size()%gp.Ndata==0);
	int n = v.size() / gp.Ndata;
	int ps = pads.size();
	assert(ps%6==0);
	assert(ps==6 || ps==n*6);
	int md = ps==6 ? 0 : 1;

	for (int d = 0; d<n; d++) {
		found = found || list_nan_inf_(&hv[d*gp.Ndata], d,
				pads[d*md  ], pads[d*md+1], pads[d*md+2],
				pads[d*md+3], pads[d*md+4], pads[d*md+5]);
	}
	if (quit && found)
		throw GridException("list_nan_inf");
}
