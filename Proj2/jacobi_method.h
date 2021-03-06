#include <armadillo>

using namespace arma;

#ifndef JACOBI_METHOD_H
#define JACOBI_METHOD_H

void jacobi_method( int N, mat &A, mat &R, vec &lambdas);

#endif // JACOBI_METHOD_H
