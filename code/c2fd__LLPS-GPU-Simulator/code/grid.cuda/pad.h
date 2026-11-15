#pragma once
#include "grid.cuda/vec.h"


void pad_reflex(Vec &v);
void pad_phi(Vec &v);
void pad_zeros(Vec &v);
void padpad_vec(const myreal *src, myreal *target, int action, int n=1);
