#pragma once
#include "grid.h"
#include <cmath>
#include <cusp/krylov/bicgstab.h>



typedef cusp::linear_operator<myreal, cusp::device_memory> LinearOperatorBase;

class LinearOperator : public LinearOperatorBase {
	public:
		LinearOperator(int num_rows, int num_cols)
			: LinearOperatorBase(num_rows, num_cols) {}
		virtual void get_RHS(Vec &rhs) const = 0;
		virtual void operator()(const Vec &x, Vec &y) const = 0;
};



int bi_cgstab(const LinearOperator &A, Vec &x, const Vec &b,
		const int iteration_limit,
		const myreal relative_tolerance,
		bool verbose=false,
		const std::string save_prefix="");

int bi_cgstab_p(const LinearOperator &A, Vec &x, const Vec &b, const LinearOperator &M,
		const int iteration_limit,
		const myreal relative_tolerance,
		bool verbose=false,
		const std::string save_prefix="");
