#pragma once
#include <string>
#include "my_hdf5.h"
#include "grid.cuda/grid_params.h"
#include "nondim.h"


namespace my_hdf5 {
    hid_t H5Fcreate_with_grid(std::string filename,
                              const GridParams &gp,
                              int round, double time, bool with_pad);

    void save_nondim(hid_t file_id, const GridParams &gp,
                     const NondimensionalizedParameters &np);

    bool   is_vizschema_file(hid_t file_id);
    int    get_pad_size(hid_t file_id);
    double get_t_0(hid_t file_id);
    int    get_step(hid_t file_id);
    double get_time(hid_t file_id);

    double get_shear_velocity(hid_t file_id);
    void set_shear_velocity(hid_t file_id, double v);
}
