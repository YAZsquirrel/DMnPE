#pragma once
#define _USE_MATH_DEFINES
#include "iostream"
#include <fstream>
#include "math.h"
#include <functional>

using namespace std;
typedef double real;


void SLAEsolve(double** M, int N, double* b, double* q);
void SLAEgen(int nX, int nY, int nZ, int nRecX, int nRecY, real* X, real* Y, real* Z, real* RecX, real* RecY);
void GaussNewton();
void FindEdgesAndCenter();
real Integrate(real xL, real xR, real yL, real yR, real zL, real zR, function<real(real, real, real, real* args, int argNum)> f, real* args, int argNum);
//void ReadRec(real* RecX, real* RecY, int& nRecX, int& nRecY);
