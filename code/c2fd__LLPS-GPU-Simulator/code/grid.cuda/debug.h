#pragma once
#include "grid.cuda/vec.h"
#include "tools.h"
#include <string>
#include <vector>

class GridException {
public:
    GridException(const std::string text="");
};

void print_memory_usage();

// Assert that last kernel invocation returned no error.
void check(std::string label="");
void check(const std::string &label, const Vec &v,
           std::vector<int> pads=std::vector<int>(6,1));
void list_nan_inf(const Vec &v, std::vector<int> pad_dir, bool quit);



// inline void dout(const string &s) {
//     cudaThreadSynchronize();
//     check();
//     cout<<s<<endl;
// }

