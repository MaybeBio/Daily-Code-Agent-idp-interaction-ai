#pragma once
#include "grid.cuda/vec.h"
#include "my_hdf5.h"


namespace globals {
    // These are managed by state.cu
    extern std::string data_dir;
    extern int file_number;
    extern hid_t h5file;
    extern bool write_pad;		// write out ghost cells
}


void write_Vec (const std::string &data_name, const Vec &v, int n=1);
void write_Vec3(const std::string &data_name, const Vec3 &vec);
void read_Vec  (const std::string &data_name, Vec &v, int n=1);
void read_Vec3 (const std::string &data_name, Vec &v);

inline void save_and_exit() {
    my_hdf5::close(globals::h5file);
    my_exit(-1);
}
