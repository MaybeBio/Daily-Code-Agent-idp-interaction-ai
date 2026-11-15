#include "my_hdf5.h"
// Hdf5 Lite (high level functions).  Link to hdf5_hl.a.
#include <hdf5.h>
#include <H5LTpublic.h>
using namespace std;

// hid_t h5file = my_H5Fcreate("data/", n);
// write_h5(h5file, "/10 Fx_before", FORTRAN_ARRAY, Fx, Nxp, Nyp, Nzp);
// my_H5Fclose(h5file);


namespace my_hdf5 {

    hid_t create(string filename) {
        return H5Fcreate(filename.c_str(),
                         H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t open(string filename) {
        hid_t ret = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (ret<0) {                // failed
            cout << "Failed to open h5 file : " << filename << endl;
            my_exit(-1);
        }
        return ret;
    }

    void close(hid_t file_id) {
        if (!file_id)
            return;
        herr_t status = H5Fclose(file_id);
        if (status<0)
            throw "Error while closing hdf5 (.h5) file.";
    }


    herr_t _list_all_links_helper(hid_t g_id, const char *name,
                                  const H5L_info_t *info, void *op_data) {
        vector<string> *v = (vector<string>*)op_data;
        v->push_back(name);
        return 0;
    }

    vector<string> list_all_links(hid_t file_id) {
        vector<string> ret;
        herr_t error = H5Literate(file_id, H5_INDEX_NAME, H5_ITER_INC,
                                  NULL, &_list_all_links_helper, &ret);
        return ret;
    }

    bool exists(hid_t hid, const char *name) {
        // TRUE (positive), FALSE (zero), or negative value if error
        htri_t ret = H5Lexists(hid, name, H5P_DEFAULT) ;
        return ret>0;
    }

} // namespace my_hdf5
