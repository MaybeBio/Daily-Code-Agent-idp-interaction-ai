#include "grid.cuda/grid_params.h"
using namespace std;


// globals
GridParams gp;


ostream& operator<<(ostream& os, const GridParams &gp) {
    os << "\n Nxyz       = " << ' ' << gp.Ny     << ' ' << gp.Nz
       << "\n dxyz       = " << ' ' << gp.dy     << ' ' << gp.dz
       << "\n Ndata,sxyz = " << gp.Ndata  << ' ' << ' ' << gp.sy << ' ' << gp.sz
       << "\n wrap_xyz   = " << ' ' << gp.wrap_y << ' ' << gp.wrap_z
       << "\n NxpNypNzp/Ndata = "<<double(gp.NypNzp)/gp.Ndata

       << "\n dimBlock   = " << gp.dimBlock.x << ' ' << gp.dimBlock.y << ' ' << gp.dimBlock.z
       << "\n dimGrid    = " << gp.dimGrid.x  << ' ' << gp.dimGrid.y  << ' ' << gp.dimGrid.z<<" \n";
    return os;
}
