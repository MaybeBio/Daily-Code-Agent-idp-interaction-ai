// the solver 
#include "state/state.h"

#include <iostream>
#include <ctime>
#include <cctype>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <limits>

#include "grid.h"
#include "parameters.h"
#include "tools.h"
#include "my_hdf5_viz.h"

using namespace std;
string save_filename("survival.txt");


void State::my_H5Fcreate_with_grid(const string &filename) {
	if (OUTPUT) {
		globals::h5file = my_hdf5::H5Fcreate_with_grid(
				filename, gp,
				round, time, globals::write_pad);
		my_hdf5::save_nondim(globals::h5file, gp, NP);
	} else {
		globals::h5file = 0;
	}
}



void State::save_some_data() const {
	write_Vec("phi", phi);
	write_Vec("qq" , qq );
	write_Vec("W" , W );
	write_Vec("R" , R );
}


void State::set_initial_condition__from_start() {
	cout<<"Initializing phi, v, c, tau."<<endl;
	myreal width = CP["interface_width"];
	string shape_phi = CP.as_string("shape_phi");
	string shape_W = CP.as_string("shape_W");
	string shape_R = CP.as_string("shape_R");

	gallery::set(phi, shape_phi, CP["init_phi"], width); pad_phi(phi);
	gallery::set(W, shape_W, CP["init_W"], width); pad_phi(W);
	gallery::set(R, shape_R, CP["init_R"], width); pad_phi(R);
}

void State::set_initial_condition() {
	string input_file = "";
	set_initial_condition__from_start();


	pad_phi(phi);
	copy_to_device();
	init_qq(qq,phi);

	qq_old = qq;
	phi_old = phi;
	W_old	= W;
	R_old	= R;


	string output_file = globals::data_dir+globals::file_number
		+ (input_file!="" ? "_loaded" : "")
		+ ".h5";
	cout<<"Saving "<<output_file<<endl;
	my_H5Fcreate_with_grid(output_file);
	save_some_data();
	my_hdf5::close(globals::h5file);
	globals::h5file = 0;
	globals::file_number++;

	std::cout<<"Init Done \n";
}



//================


void State::my_bicgstab(const string &name, const LinearOperator &linop, Vec &xx) {
	// Solve AX = B  -->  linop xx = bb
	//
	// Two options for Bi-CGSTAB
	//   1) Boundary needs to be zero,
	//   2) Bi-CGSTAB's dot(,) needs exclude padding from norm calculation.
	// We go with (1).
	// 
	// 2011-11-13: The following is no longer true?
	// // It turns out that bicgstab only use A^n xx where n>=1.
	// // Junk pad in, junk pad out.
	pad_zeros(xx);

	Vec bb(xx.size());
	linop.get_RHS(bb);

	if (save_it && (OUTPUT & OUTPUT_BICG_STEPS)) {
		write_Vec(name+"_rhs", bb);
	}


	clock_t lin_s, lin_e;
	lin_s = clock();
	int rounds = bi_cgstab(
			linop, xx, bb,
			CP.as_int ("bicg_iteration_limit"),
			CP        ["bicg_relative_tolerance_" + name],
			CP.as_bool("bicg_verbose"),
			name);
	lin_e = clock();
	check();

	if (save_it) {
		double tdiff = (lin_e-lin_s) / double(CLOCKS_PER_SEC);
		double _max = max(xx);
		double _min = min(xx);
		double _sum = sum(xx,0) * gp.dy * gp.dz;

		cout << name << " solved in " << scientific << setprecision(3)
			<< tdiff << " secs (" << rounds << " rounds)"
			<< "  Max,min,sum = " << _max
			<< ", " << _min
			<< ", " << _sum
			<< endl;
		// When we have a NaN,
		// _max, _min might return 0.0 while _sum returns NaN.
		// This is because comparison with NaN yields false. 
		if (!isfinite(_sum))
			my_exit(-1);
	}

	if (CP.as_bool("do_clip"))
		clip_and_print(xx, 0., 1., name, save_it);

	// Note:  cc's padded values are not used.  So we might as well not pad it.
	pad_reflex(xx);
}

void State::my_bicgstab_p(const string &name, const LinearOperator &linop, Vec &xx, const LinearOperator &M) {

	pad_zeros(xx);

	Vec bb(xx.size());
	linop.get_RHS(bb);

	if (save_it && (OUTPUT & OUTPUT_BICG_STEPS)) {
		write_Vec(name+"_rhs", bb);
	}


	clock_t lin_s, lin_e;
	lin_s = clock();
	int rounds = bi_cgstab_p(
			linop, xx, bb, M,
			CP.as_int ("bicg_iteration_limit"),
			CP        ["bicg_relative_tolerance_" + name],
			CP.as_bool("bicg_verbose"),
			name);
	lin_e = clock();
	check();

	if (save_it) {
		double tdiff = (lin_e-lin_s) / double(CLOCKS_PER_SEC);
		double _max = max(xx);
		double _min = min(xx);
		double _sum = sum(xx,0) * gp.dy * gp.dz;

		cout << name << " solved in " << scientific << setprecision(3)
			<< tdiff << " secs (" << rounds << " rounds)"
			<< "  Max,min,sum = " << _max
			<< ", " << _min
			<< ", " << _sum
			<< endl;
		if (!isfinite(_sum))
			my_exit(-1);
	}

	if (CP.as_bool("do_clip") && (name.find("velo") == string::npos))
		clip_and_print(xx, 0., 1., name, save_it);

}






void print_max_velo(const string &name, Vec3 &v) {
	cout << "      maxabs " << name << " x/y/z  : "
		<< scientific << setprecision(3)
		<< maxabs(v, 0)
		<< ", " << maxabs(v, 1)
		<< ", " << maxabs(v, 2)
		<< endl;
}


//================================================================




void State::simulate() {
	cout<<"simulate()"<<endl;
	// Values for CP["save"]
	//    0 -- no save
	//   -5 -- save every 5 steps
	//   +5.0 -- save every t = 5 * t_0.
	int save_stride = 0;
	double save_interval = CP["save"];
	double next_save_time = time + save_interval - dt/2.;
	if (save_interval <= 0.) {
		save_interval = 0.;
		next_save_time = 1e30;	// We will never get there (tm).
		save_stride  = -CP.as_int("save");
		// save==0 --> save_[interval/stride] == 0
	}

	//E0 = get_total_energy();
	//E_now = E0, E_old = E0;
	//std::cout<<"Energy "<<E_now-E_old<<"\t "<<E_now-E0<<"\t "<<E_now<<"\n";
	save_some_model_info();


	double t_end = CP["t_end"];
	int round_end = numeric_limits<int>::max() / 2;
	if (t_end<0.) {
		t_end = numeric_limits<double>::max() / 2.;
		round_end = -CP.as_int("t_end");
	}

	while(time + dt/2 < t_end && round < round_end) {
		round++;
		time += dt;

		if (save_stride) {
			save_it = ((round-round_initial)%save_stride == 0);
		} else {
			save_it = false;
			if (time >= next_save_time) {
				save_it = true;
				next_save_time += save_interval;
			}
		}
		if (round<CP["first_save_round"])
			save_it = false;

		if (save_it) {
			if (system("date"))
				cout<<" ^^^^ (>_<)~~~~ "<<endl;
			cerr<<flush;        // date output
			cout<<"\nround = "<<round<<",  time = "<<fixed<<time<<endl;
			my_H5Fcreate_with_grid(globals::data_dir+globals::file_number+".h5");
		}

		clock_t start;
		if (save_it)  start = clock();
		one_step();

		if (save_it) {
			clock_t end = clock();
			print_memory_usage();
			double tdiff = (end - start) / double(CLOCKS_PER_SEC);
			cout << "\n   took " << tdiff << " seconds."
				<< "\n-----------------------------------------------------\n" <<endl;
		}
		if (save_it) {
			// if(time < CP.t_cutoff_c) 
			//     don't write c_n
			my_hdf5::close(globals::h5file);
			globals::h5file = 0;
			globals::file_number++;
		}
	}


	// // All's well that ends well
	// globals::h5file = my_H5Fcreate_with_grid(globals::data_dir+globals::file_number+".h5");
	// save_some_data();
	// my_H5Fclose(globals::h5file);
}






//================================================================



void State::one_step() {
	clock_t lin_s, lin_e;
	if (save_it)
		lin_s = clock();

	copy_to_device();    // Do we still need this?

	project(phi,phi_old,2.,phi_bar); pad_phi(phi_bar);
	update_parqq();

	project(W,W_old,2.,W_bar);
	project(R,R_old,2.,R_bar);
	project(phi,phi_old,2.,phi_np1);



	{
		PhiBDF2Mat linop(this);
		//PrePhi PrP(NP.Lambda*NP.Gamma_1,0.,3./2./dt);
		//my_bicgstab_p("phi", linop, phi_np1,PrP);
		my_bicgstab("phi", linop, phi_np1);

		check("phi_np1", phi_np1);
		pad_phi(phi_np1);
		update_qq();

		swap(phi_old, phi);
		swap(phi, phi_np1);
	}

	{
		blas::copy(W_bar,W_np1);
		WMat linop(this);
		my_bicgstab("W",linop,W_np1);
		pad_phi(W_np1);
		check("WMat");

		swap(W_old, W);
		swap(W, W_np1);
	}

	{
		blas::copy(R_bar,R_np1);
		RMat linop(this);
		my_bicgstab("R",linop,R_np1);
		pad_phi(R_np1);
		check("RMat");

		swap(R_old, R);
		swap(R, R_np1);
	}






	if (save_it) {
		lin_e = clock();
		double tdiff = (lin_e-lin_s) / double(CLOCKS_PER_SEC);
		Vec tmp(gp.Ndata,0.);blas::copy(phi,tmp);
		pad_zeros(tmp);
		double volume = sum(tmp)*gp.dy*gp.dz;
		cout << "CH solved in " << scientific << setprecision(3) << tdiff<<" \t "<<volume<<" \n";
	}
	if(save_it){
		save_some_data();
		save_some_model_info();
	}
}

void State::save_some_model_info(){
	//string data_file_name(globals::data_dir+save_filename);
	//ofstream o_file(data_file_name.c_str(),ios::app);

	//E_now = get_total_energy();
	//std::cout<<"Energy "<<E_now-E_old<<"\t "<<E_now-E0<<"\t "<<E_now<<"\n";

	//o_file<<scientific<<setprecision(6)<<"Time: E2-E1/E2-E0/E2: \t "<<time<<" \t "<<E_now-E_old<<" \t "<<E_now-E0<<" \t "<<E_now<<" \n ";
	//o_file.close();
	//E_old = E_now;
}
