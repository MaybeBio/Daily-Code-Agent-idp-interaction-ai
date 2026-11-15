#pragma once

#include "grid.h"
#include "nondim.h"


const int Nphi = 1;


class StateBase {
	public:
		Parameters CP;

		Vec  phi,  phi_old,  phi_bar, phi_np1;
		Vec  qq, qq_old, qq_np1, parqq;

		Vec  W, W_old, W_bar, W_np1;
		Vec  R, R_old, R_bar, R_np1;

		void init(size_t Ndata, int Ndim) {
			// These 0's are for padpad, so they don't mess up min/max/sum.

			assert(Nphi == 1);
			phi     .resize(Ndata*Nphi, 0.);
			phi_old .resize(Ndata*Nphi, 0.);
			phi_bar .resize(Ndata*Nphi, 0.);
			phi_np1 .resize(Ndata*Nphi, 0.);

			qq	.resize(Ndata*Nphi, 0.);
			qq_old	.resize(Ndata*Nphi, 0.);
			qq_np1	.resize(Ndata*Nphi, 0.);
			parqq	.resize(Ndata*Nphi, 0.);

			W	.resize(Ndata, 0.);
			W_old	.resize(Ndata, 0.);
			W_bar	.resize(Ndata, 0.);
			W_np1	.resize(Ndata, 0.);

			R	.resize(Ndata, 0.);
			R_old	.resize(Ndata, 0.);
			R_bar	.resize(Ndata, 0.);
			R_np1	.resize(Ndata, 0.);
		}
};


class StateBaseRaw {
	public:
		myreal *phi,  *phi_old,  *phi_bar,  *phi_np1;
		myreal *qq,  *qq_old, *qq_np1, *parqq;
		myreal *W, *W_old, *W_bar, *W_np1;
		myreal *R, *R_old, *R_bar, *R_np1;

		myreal dt, idt, idt3_2, idt_2;

		void set(StateBase &s, double dt_) {
			phi     = cast(s.phi);
			phi_old = cast(s.phi_old);
			phi_bar = cast(s.phi_bar);
			phi_np1 = cast(s.phi_np1);

			qq	= cast(s.qq);
			qq_old  = cast(s.qq_old);
			qq_np1  = cast(s.qq_np1);
			parqq	= cast(s.parqq);

			W	= cast(s.W);
			W_old	= cast(s.W_old);
			W_bar	= cast(s.W_bar);
			W_np1	= cast(s.W_np1);

			R	= cast(s.R);
			R_old	= cast(s.R_old);
			R_bar	= cast(s.R_bar);
			R_np1	= cast(s.R_np1);

			dt     = dt_;
			idt    = 1./dt;
			idt3_2 = idt*3./2;
			idt_2  = idt/2.;
		}
};
