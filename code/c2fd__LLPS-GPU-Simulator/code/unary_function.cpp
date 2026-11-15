#include "tools.h"
using namespace std;

void UnaryFunction::init(string s, double h_0, double t_0) {
    cout << "### " << s << endl;
    istringstream ss(s);
    ss.exceptions(ios_base::failbit);
    ss>>name;
    bool dim = false;           // input are dimensionalized values
    try {
        if (name=="dim") {
            dim = true;
            ss>>name;
        }
	if (name=="constant") {
	    // argv: a
	    argv.resize(1);
	    ss>>argv[0];
            if (dim) {
                argv[0] /= h_0/t_0;
            }
	    return;
	}
	if (name=="ramp") {
	    // argv: v t t_end
	    //
	    // time               velocity
	    // -------------------------------
	    // (0    ,t    )      0 -> v (linearly)
	    // (t    ,t_end)      v
	    // (t_end,t_end+t)    v -> 0
	    argv.resize(3);
	    ss.exceptions(ios_base::goodbit); // clear flags
	    argv[0] = next_arg(ss, 1.);  // velocity (v)
	    argv[1] = next_arg(ss, 1.);  // ramp time (t)
	    argv[2] = next_arg(ss, 1e9); // end time (t_end)
            if (dim) {
                argv[0] /= h_0/t_0;
                argv[1] /= t_0;
                argv[2] /= t_0;
            }
	    return;
	}
	if (name=="sine" || name=="square") {
	    // a*sin(time/t)  for n periods
	    argv.resize(3);
	    ss.exceptions(ios_base::goodbit); // clear flags
	    argv[0] = next_arg(ss, 1.);  // amplitude (a)
	    argv[1] = next_arg(ss, 40.); // period (t)
	    argv[2] = next_arg(ss, 1.);  // number of periods (n)
            if (dim) {
                argv[0] /= h_0/t_0;
                argv[1] /= t_0;
            }
	    return;
	}
    } catch (...) {
	cerr << "Error parsing UnaryFunction.\n"
	     << "... Did you provide enough function arguments?" <<endl;
	my_exit(-1);
    }
    cout << "Unknown UnaryFunction name : "<<name;
    my_exit(-1);
}


double square(double t, double tt) {
    // (almost) square wave with period 1, amplitude 1
    t -= floor(t);
    if (t < tt)  return (t/tt);
    if (t < .5-tt)  return 1.;
    if (t < .5+tt)  return (.5-t)/tt;
    if (t < 1.-tt)  return -1;
    return (t-1.)/tt;
}



double UnaryFunction::operator()(double time) {
    if (name=="constant") {
	return argv[0];
    } else if (name=="ramp") {
	if (time <= argv[1])
	    return argv[0] * time / argv[1];
	if (time <= argv[2])
	    return argv[0];
	if (time <= argv[2]+argv[1])
	    return argv[0] * (argv[2]+argv[1] - time)/argv[1];
	else
	    return 0.;
    } else if (name=="sine") {
        double tt = min(time/argv[1], argv[2]);
        return argv[0]*sin(2.*PI*tt);
    } else if (name=="square") {
        double tt = min(time/argv[1], argv[2]);
	return argv[0]*square(tt, .025); // rise time = 5% of the period.
    }
    cout << "UnaryFunction::operator() : unknown name : "<<name<<endl;
    my_exit(-1);
}
void UnaryFunction::print() const {
    cout.unsetf(ios_base::floatfield);
    cout<<"nondimensinalized shear : " << name << "  ";
    for (int i=0; i<argv.size(); i++)
	cout<<argv[i]<<' ';
    cout<<endl;
}


// int main() {
//     UnaryFunction u;
//     u.init("ramp 7.2 10000", 1e-3, 1e2);
//     cout<<u(10)<<' '<<u(100)<<' '<<u(200)<<endl;
//     return 0;
// }
