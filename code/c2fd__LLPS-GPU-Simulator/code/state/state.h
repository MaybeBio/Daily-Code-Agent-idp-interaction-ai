#pragma once

#include <cmath>
#include <functional>
#include <string>

#include "state_base.h"
#include "grid.h"
#include "grid.cuda/bicgstab.h"
#include "unary_function.h"


extern StateBaseRaw st_raw;



void print_max_velo(const std::string &name, Vec3 &v);


class State : public StateBase {
	public:
		NondimensionalizedParameters NP;

		int round;
		int round_initial;
		double time;
		double dt;

		//myreal E0,E_now,E_old;

		bool save_it;

		void init() {
			assert(CP.inited());
			// We should also check that "NP" and "gp" are inited.

			StateBase::init(gp.Ndata, Ndim);

			round = 0;
			round_initial = 0;
			time = 0.;
			dt = CP["dt"];
			assert (gp.Ndata);
		}

		void my_H5Fcreate_with_grid(const std::string &filename);

	private:
		void set_initial_condition__from_start();
	public:
		void set_initial_condition();
		void finallize_conditions();

		void copy_to_device();
		void copy_NP_to_device();
		void save_some_data() const;

		void simulate();
		void one_step();
		void my_bicgstab(const std::string &name, const LinearOperator &linop, Vec &xx);
		void my_bicgstab_p(const std::string &name, const LinearOperator &linop, Vec &xx, const LinearOperator &M);



		void update_qq();
		void update_parqq();
		void init_qq(Vec &myqq, const Vec &myphi);
		//myreal get_total_energy();
		void save_some_model_info();

};








class PhiBDF2Mat : public LinearOperator {
	public:
		State *s;
		PhiBDF2Mat(State *s)
			: LinearOperator(gp.Ndata*Nphi, gp.Ndata*Nphi), s(s) {}
		void get_RHS(Vec &rhs) const;
		void operator()(const Vec &x, Vec &y) const;
};



class PrePhi: public LinearOperator {
	public:
		std::string do_name;
		double sub_p;
		double sub_q;
		double sub_lam;
		PrePhi(const myreal sub_p, const myreal sub_q, const myreal sub_lam):
			LinearOperator(gp.Ndata,gp.Ndata),sub_p(sub_p),sub_q(sub_q),sub_lam(sub_lam){
			}
		void get_RHS(Vec &rhs) const;
		void operator()(const Vec &x, Vec &y) const;
};

class WMat : public LinearOperator {
	public:
		State *s;
		WMat(State *s)
			: LinearOperator(gp.Ndata, gp.Ndata), s(s) {}
		void get_RHS(Vec &rhs) const;
		void operator()(const Vec &x, Vec &y) const;
};
class RMat : public LinearOperator {
	public:
		State *s;
		RMat(State *s)
			: LinearOperator(gp.Ndata, gp.Ndata), s(s) {}
		void get_RHS(Vec &rhs) const;
		void operator()(const Vec &x, Vec &y) const;
};

