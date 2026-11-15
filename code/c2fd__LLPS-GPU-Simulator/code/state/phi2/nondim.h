#pragma once
#include "parameters.h"
#include "tools.h"


class NondimensionalizedParameters {
	public:
		myreal t_0;
		myreal h;
		myreal Gamma_1; //used in the chemican potential
		myreal Lambda;  //relaxation factor of the growth of phi
		myreal LambdaW, LambdaR;
		myreal StableC0;
		myreal c1,c2;

		NondimensionalizedParameters() {};
		void init(const Parameters &CP, const Parameters &GP);
		void display() const;
};

