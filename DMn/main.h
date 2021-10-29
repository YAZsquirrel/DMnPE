#pragma once
#define _USE_MATH_DEFINES
#include "iostream"
#include <fstream>
#include "math.h"
#include <functional>

using namespace std;
typedef double real;

real Integrate(real xL, real xR, real yL, real yR, real zL, real zR, function<real(real, real, real, real* args, int argNum)> f, real* args, int argNum);
//void ReadRec(real* RecX, real* RecY, int& nRecX, int& nRecY);
