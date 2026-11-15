#include "grid.cuda/bicgstab.h"
using namespace std;

string bicg_save_prefix;        // used in cusp/krylov/bicgstab.inl


int bi_cgstab(const LinearOperator &A, Vec &x, const Vec &b,
		const int iteration_limit,
		const myreal relative_tolerance,
		bool verbose,
		string save_prefix) {
	const int &lim = iteration_limit;
	const myreal &tol = relative_tolerance;

	if ((OUTPUT & OUTPUT_BICG_STEPS) && globals::h5file)
		bicg_save_prefix = save_prefix;
	else
		bicg_save_prefix = "";


	if (verbose) {
		cusp::verbose_monitor<myreal> monitor(b, lim, tol);
		cusp::krylov::bicgstab(A, x, *const_cast<Vec *>(&b), monitor);
		return monitor.iteration_count();
	} else {
		cusp::default_monitor<myreal> monitor(b, lim, tol);
		cusp::krylov::bicgstab(A, x, *const_cast<Vec *>(&b), monitor);
		return monitor.iteration_count();
	}
}


int bi_cgstab_p(const LinearOperator &A, Vec &x, const Vec &b, const LinearOperator &M,
		const int iteration_limit,
		const myreal relative_tolerance,
		bool verbose,
		string save_prefix) {
	const int &lim = iteration_limit;
	const myreal &tol = relative_tolerance;

	if ((OUTPUT & OUTPUT_BICG_STEPS) && globals::h5file)
		bicg_save_prefix = save_prefix;
	else
		bicg_save_prefix = "";


	if (verbose) {
		cusp::verbose_monitor<myreal> monitor(b, lim, tol);
		cusp::krylov::bicgstab(A, x, *const_cast<Vec *>(&b), monitor,M);
		return monitor.iteration_count();
	} else {
		cusp::default_monitor<myreal> monitor(b, lim, tol);
		cusp::krylov::bicgstab(A, x, *const_cast<Vec *>(&b), monitor,M);
		return monitor.iteration_count();
	}
}

