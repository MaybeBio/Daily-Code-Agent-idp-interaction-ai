#include "grid.cuda/io.h"
#include "grid.cuda/grid_params.h"
#include "grid.cuda/pad.h"
#include "grid.cuda/debug.h"
#include "my_hdf5_viz.h"
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;


namespace globals {
    std::string data_dir = "";
    int file_number = 0;
    hid_t h5file = 0;
    bool write_pad = false;
}


// template<int Dim> void write_array(const string &filename, const myreal *u,
// 				   bool transpose=false);

// template<>
// void write_array<2>(const string &filename, const myreal *u, bool transpose) {

//     // Important: 2d-coord is (y,z), not (x,y)
//     const int ni = gp.Nyp, si = gp.sy;
//     const int nj = gp.Nzp, sj = gp.sz;

//     ofstream ofile(filename.c_str(), ios::out);
//     ofile.flags(ios::scientific | ios::showpos | ios::showpoint );
//     ofile.precision(16);   
//     if(!ofile) {
// 	cout << "Cannot open output file : " << filename <<endl;
// 	my_exit(-1);
//     }

//     if (transpose) {
// 	for(int j=0; j<nj; j++) {
// 	    for(int i=0; i<ni; i++)
// 		ofile << u[i*si+j*sj] << "  ";
// 	    ofile << endl;
// 	}
//     } else {
// 	for(int i=0; i<ni; i++) {
// 	    for(int j=0; j<nj; j++)
// 		ofile << u[i*si+j*sj] << "  ";
// 	    ofile << endl;
// 	}
//     }
//     ofile.close();
// }


// template<>
// void write_array<3>(const string &filename, const myreal *u, bool transpose) {

//     // write out only one 2d-slice u[x,*,*]
//     int x = gp.Nxp/2;
//     write_array<2>(filename, u + x*gp.sx, transpose);
// }

//     // string filename = globals::data_dir + data_name + "_" + globals::file_number;
//     // write_array<Ndim>(filename, cast(hv));



//================================================================


// 2010-09-23
// The Visit's Vs reader doesn't work on compMajorC/F mode.
// Stick with scalar var.

// void write_Vec(const string &data_name, const Vec &v, int n) {
//     H_Vec hv(v);
//     if (n==0)
// 	write_h5(globals::h5file, data_name, cast(hv), gp.Nxp, gp.Nyp, gp.sy);
//     else
// 	write_h5(globals::h5file, data_name, cast(hv), n, gp.Nxp, gp.Nyp, gp.sy);
// }

// void read_Vec(const string &data_name, Vec &v, int n) {
//     // Note: n is infered from the hdf5 file.
//     cout<<"Reading "<<data_name<<endl;
//     H_Vec hv(v.size());
//     read_h5(globals::h5file, data_name, cast(hv));
//     v = hv;
// }


void write_Vec(const string &data_name, const Vec &v, int n) {
    assert(n);
    int ny, nz;
    if (globals::write_pad) {
	ny = gp.Nyp;
	nz = gp.Nzp;
    } else {
	ny = gp.Ny;
	nz = gp.Nz;
    }
    int ndata = ny*nz;

    for (int ii=0; ii<n; ii++) {
	Vec stripped(ndata);
	int action = globals::write_pad ? -1 : -2;
	padpad_vec(cast(v)+ii*gp.Ndata, cast(stripped), action);
	check();

	H_Vec hv(stripped);
	string name = data_name + (n==1 ? "" : itoa(ii));
        my_hdf5::write(globals::h5file, name, cast(hv), ny, nz);
    }
}


void read_Vec(const string &data_name, Vec &v, int n) {
    assert(n);
    int ndata = gp.NypNzp;
    H_Vec hv(ndata);

    cout<<"Reading "<<data_name<<" x "<<n<<endl;
    int pad_size = my_hdf5::get_pad_size(globals::h5file);
    cout<<"... pad_size = "<<pad_size<<endl;

    assert(pad_size==0 || pad_size==1);
    int action = (pad_size==1) ? 1 : 2;

    for (int ii=0; ii<n; ii++) {
	string name = data_name + (n==1 ? "" : itoa(ii));
        my_hdf5::read(globals::h5file, name, cast(hv));

	Vec stripped(hv);
	padpad_vec(cast(stripped), cast(v)+ii*gp.Ndata, action);
    }
}


void write_Vec3(const string &data_name, const Vec3 &v) {
    write_Vec(data_name, v, 3);
}
void read_Vec3(const string &data_name, Vec3 &v) {
    read_Vec(data_name, v, 3);
}

