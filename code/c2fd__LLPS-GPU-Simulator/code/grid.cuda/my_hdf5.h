// 
// Export data to hdf5 (.h5) file.
// Use VizScheme convention.
// 

#pragma once
#include "tools.h"		// for DEBUG

#include <hdf5.h>
#include <H5LTpublic.h>

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>


namespace my_hdf5 {
    hid_t create(std::string filename);
    hid_t open(std::string filename);
    void  close(hid_t file_id);

    std::vector<std::string> list_all_links(hid_t file_id);
    bool exists(hid_t hid, const char *name);


// inline hid_t h5_type(std::string data) { return H5T_C_S1; }
// inline hid_t h5_type(const char *data)  { return H5T_NATIVE_CHAR; }
    inline hid_t h5_type(const short *data)  { return H5T_NATIVE_SHORT; }
    inline hid_t h5_type(const int *data)    { return H5T_NATIVE_INT;   }
    inline hid_t h5_type(const long *data)   { return H5T_NATIVE_LONG;  }
    inline hid_t h5_type(const float *data)  { return H5T_NATIVE_FLOAT; }
    inline hid_t h5_type(const double *data) { return H5T_NATIVE_DOUBLE;}



    template <class T>
    void write(hid_t file_id, const std::string &dataset_name,
               const T *data,	// C-ordered
               size_t n0, size_t n1=0, size_t n2=0, size_t n3=0)
    {
        if (!file_id) {
            std::cout<<"write_h5: file_id is NULL"<<std::endl;
            my_exit(-1);
        }
        if (DEBUG & DEBUG_VERBOSE_H5)
            std::cout<<"write_h5: "<<dataset_name<<std::endl;

        hsize_t dims[4] = {n0, n1, n2, n3};

        int rank = 1;
        if (n1>0) rank = 2;
        if (n2>0) rank = 3;
        if (n3>0) rank = 4;

        const char *name = dataset_name.c_str();
    
        herr_t status = H5LTmake_dataset(
            file_id, name, rank, dims, h5_type(data), data);

        if (status<0)
            throw "Error while writing hdf5 (.h5) file.";

        H5LTset_attribute_string (file_id, name, "vsType", "variable");
        H5LTset_attribute_string (file_id, name, "vsMesh", "mygrid");
        H5LTset_attribute_string (file_id, name, "vsTimeGroup", "mytime");

        // Note: as of 2010-09-23
        // The Visit's Vs reader doesn't work on compMajorC/F mode.
        // Stick with scalar var.
        assert(rank<=3);
        // H5LTset_attribute_string (file_id, name, "vsIndexOrder", "compMajorC");


        // Note: velo, tau data aren't at grid center.
        //
        // https://ice.txcorp.com/trac/vizschema/wiki/Variables
        //
        // > In the future, we hope to support face and edge, but for now,
        // > any data that is not labeled specifically "zonal" is assumed
        // > to be nodal.
        //
        H5LTset_attribute_string (file_id, name, "vsCentering", "zonal");
    }



    template <class T>
    void read(hid_t file_id, const std::string &dataset_name,
              T *data)	// C-ordered
    {
        if (DEBUG & DEBUG_VERBOSE_H5)
            std::cout<<"read_h5: "<<dataset_name<<std::endl;

        herr_t status = H5LTread_dataset(
            file_id, dataset_name.c_str(), h5_type(data), (void*)(cast(data)));

        if (status<0)
            throw "Error while writing hdf5 (.h5) file.";
    }
} // namespace my_hdf5
