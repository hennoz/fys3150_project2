#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <time.h>
#include <ctime>
#include "catch.hpp"
#include "Tridiag.h"
#include "jacobi.h"


using namespace std;
using namespace arma;


void Jacobi_rotate( mat &A, mat &R, int k, int l, int N);
double maxoffdiag ( mat A, int *k, int *l, int N);
inline double analytic_eigenvals(int a, int d, int N,int j){return d + 2*a*cos(j*acos(-1.0)/(N+1)) ;}

#endif // MAIN_H
