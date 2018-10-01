#include <armadillo>
using namespace arma;

#ifndef TEST_H
#define TEST_H

void test_eigvals(int N, vec &lambdas, vec &a_lambdas, double a, double d);
void test_eigvals_quantum( vec &lambdas, double p_N);
void test_ortho( mat R, int N );


#endif // TEST_H
