#pragma once
#include "iostream"
#include <fstream>
#include "math.h"
#include <functional>

using namespace std;
typedef double real;

void SLAEsolve();
void SLAEgen();
void GaussNewton();
void FindEdgesAndCenter();
real Integrate(real xL, real xR, real yL, real yR, real zL, real zR, function<real(real, real, real)> f);
//void ReadRec(real* RecX, real* RecY, int& nRecX, int& nRecY);
