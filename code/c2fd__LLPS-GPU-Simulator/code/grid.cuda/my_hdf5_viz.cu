#include "my_hdf5_viz.h"
using namespace std;

namespace my_hdf5 {



    void set_attribute(hid_t id, const char *attr_name, const char *val) {
        H5LTset_attribute_string(id, ".", attr_name, val);
    }
    template <typename T>
    void set_attribute(hid_t id, const char *attr_name, const T val) {
        hid_t space_id  = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(id, attr_name, h5_type(&val),
                                space_id, H5P_DEFAULT, H5P_DEFAULT);
        assert(attr >= 0);
        herr_t ret = H5Awrite(attr, h5_type(&val), &val);
        assert(ret >= 0);
        ret = H5Aclose(attr);
        ret = H5Sclose(space_id);
    }
    void set_attribute(hid_t id, const char *attr_name,
                       int v1, int v2) {
        int val[2] = {v1,v2};
        H5LTset_attribute_int(id, ".", attr_name, val, 2);
    }
    void set_attribute(hid_t id, const char *attr_name,
                       double v1, double v2) {
        double val[2] = {v1,v2};
        H5LTset_attribute_double(id, ".", attr_name, val, 2);
    }


    template <typename T>
    void get_attribute(hid_t id, const char *attr_name, T &val) {
        hid_t aid  = H5Screate(H5S_SCALAR);
        hid_t attr = H5Aopen(id, attr_name, H5P_DEFAULT);
        assert(attr>=0);
        herr_t status = H5Aread(attr, h5_type(&val), &val);
        assert(status>=0);
        status = H5Aclose(attr); assert(status>=0);
        status = H5Sclose(aid);  assert(status>=0);
    }
    template <typename T>
    void get_attribute(hid_t base_id, const char *path,
                       const char *attr_name, T &val) {
        hid_t id = H5Gopen2(base_id, path, H5P_DEFAULT);
        assert(id >= 0);
        get_attribute(id, attr_name, val);
        herr_t status = H5Gclose(id);
        assert(status >= 0);
    }



    // ================

    template <typename T>
    void set_data(hid_t id, const char *attr_name, const T val) {
        hid_t space_id  = H5Screate(H5S_SCALAR);
        hid_t attr = H5Dcreate2(id, attr_name, h5_type(&val),
                                space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        assert(attr >= 0);
        herr_t ret = H5Dwrite(attr, h5_type(&val),
                              space_id, space_id,
                              H5P_DEFAULT, &val);
        assert(ret >= 0);
        ret = H5Dclose(attr);
        ret = H5Sclose(space_id);
    }
    template <typename T>
    void get_data(hid_t id, const char *attr_name, T& val) {
        hid_t space_id  = H5Screate(H5S_SCALAR);
        hid_t attr = H5Dopen2(id, attr_name, H5P_DEFAULT);
        assert(attr >= 0);
        herr_t ret = H5Dread(attr, h5_type(&val),
                             space_id, space_id,
                             H5P_DEFAULT, &val);
        assert(ret >= 0);
        ret = H5Dclose(attr);
        ret = H5Sclose(space_id);
    }



    // ================================================================




    hid_t H5Fcreate_with_grid(string filename,
                              const GridParams &gp,
                              int round, double time, bool with_pad) {
        hid_t file_id = create(filename);
        hid_t id;

        int    ny = gp.Ny, nz = gp.Nz;
        double dy = gp.dy, dz = gp.dz;

        // using namespace my_h5;
        id = H5Gcreate2(file_id, "mygrid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // size_t size_hint = 65 + 9; // including '\0' -- (6+6+11+10+13+13)+6 = 65
        set_attribute(id, "vsType", "mesh");
        set_attribute(id, "vsKind", "uniform");
        set_attribute(id, "vsStartCell"  , 0, 0);
        if (with_pad) {
            set_attribute(id, "vsNumCells"   , ny+2, nz+2);
            set_attribute(id, "vsLowerBounds", -dy, -dz);
           set_attribute(id, "vsUpperBounds", (ny+1)*dy, (nz+1)*dz);
        } else {
            set_attribute(id, "vsNumCells"   , ny, nz);
            set_attribute(id, "vsLowerBounds", 0., 0.);
            set_attribute(id, "vsUpperBounds", ny*dy, nz*dz);
        }
        set_attribute(id, "pad_size", with_pad ? 1 : 0);
        H5Gclose(id);


        id = H5Gcreate2(file_id, "mytime", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        // size_hint = 28 + 4;		// 7*4 = 28
        // https://ice.txcorp.com/trac/vizschema/wiki/OtherMetaData
        // 2011-02-18: Verified this with visit/src/databases/Vs/*
        set_attribute(id, "vsType", "time");
        set_attribute(id, "vsStep", round);
        set_attribute(id, "vsTime", time);
        H5Gclose(id);

        // // Don't know how to make time show up in VisIt
        // // through VisSchema reader.  We will just save it as a data.
        // set_data(file_id, "time", time);

        return file_id;
    }

    void save_nondim(hid_t file_id, const GridParams &gp,
                     const NondimensionalizedParameters &np) {
        hid_t id = H5Gcreate2(file_id, "mynondim", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        set_attribute(id, "Ny", gp.Ny);
        set_attribute(id, "Nz", gp.Nz);
        set_attribute(id, "dy", gp.dy);
        set_attribute(id, "dz", gp.dz);

	/*
        set_attribute(id, "Gamma_1", np.Gamma_1);
        set_attribute(id, "Gamma_2", np.Gamma_2);
        set_attribute(id, "Lambda", np.Lambda);
        set_attribute(id, "Kai_CH", np.Kai_CH);
        set_attribute(id, "NP_inv", np.NP_inv);

        set_attribute(id, "eta_net", np.eta_net);
        set_attribute(id, "eta_sol", np.eta_sol);
        set_attribute(id, "rho", np.rho);
	*/
        
        H5Gclose(id);
    }


    bool is_vizschema_file(hid_t file_id) {
        return exists(file_id, "mytime") && 
            exists(file_id, "mygrid") && 
            exists(file_id, "mynondim");
    }


    double get_t_0(hid_t file_id) {
        double ret;
        get_attribute(file_id, "/mynondim", "t_0", ret);
        return ret;
    }


    int get_pad_size(hid_t file_id) {
        if (!exists(file_id, "mygrid")) {
            cout << "The file doesn't have the attribute /mygrid.pad_size\n"
                 << "Assume that pad_size == 0" << endl;
            return 0;
        }
        int ret;
        get_attribute(file_id, "/mygrid", "pad_size", ret);
        return ret;
    }

    int get_step(hid_t file_id) {
        int ret;
        get_attribute(file_id, "/mytime", "vsStep", ret);
        return ret;
    }

    double get_time(hid_t file_id) {
        double ret;
        get_attribute(file_id, "/mytime", "vsTime", ret);
        return ret;
    }

    double get_shear_velocity(hid_t file_id) {
        double ret;
        get_data(file_id, "shear_velocity", ret);
        return ret;
    }
    void set_shear_velocity(hid_t file_id, double v) {
        set_data(file_id, "shear_velocity", v);
    }

} // namespace my_hdf5
