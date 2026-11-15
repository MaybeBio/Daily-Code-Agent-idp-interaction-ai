#include "nondim.h"
#include <iostream>
#include <cassert>


void NondimensionalizedParameters::init(
    const Parameters &CP, const Parameters &GP) {

    Gamma_1	= GP["gam_1"];
    Lambda	= GP["lam"];
    LambdaW	= GP["lamW"];
    LambdaR	= GP["lamR"];
    StableC0	= GP["StableC0"];
    c1		= GP["c1"];
    c2		= GP["c2"];

}




void NondimensionalizedParameters::display() const {
    using namespace std;
    cout << "\n\n ************* Non-dimensionalized global parameters *************\n";
    cout << "\n  Gamma_1 = " << Gamma_1
	 << "\n  Lambda  = " << Lambda
	 << "\n  LambdaW  = " << LambdaW
	 << "\n  LambdaR  = " << LambdaR
	 << "\n  StableC0     = " << StableC0
	 << "\n  c1     = " << c1
	 << "\n  c2     = " << c2

	 << endl;
    cout << endl;
}

