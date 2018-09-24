//#include <iostream>
//#include <cmath>
#include <armadillo>
//#include "catch.hpp"
//#include "jacobi_rotate.h"
//#include "maxoffdiag.h"
#include "analytic_eigvals.h"


using namespace std;
using namespace arma;

void test_eigvals(int N, vec lambdas, vec a_lambdas)
{

    for ( int i = 0; i < N; i++)
    {
         a_lambdas[i] = analytic_eigenvals(-1,2,N,i+1);
//         cout << "   " << a_lambdas[i] << endl;
    }
    cout << "TEST_EIGVALS RUNS ...  ";
    int c = 0;
    double tol = 1e-10;
    for ( int i = 0; i < N; i++)
    {
        if ( fabs( lambdas[i] - a_lambdas[i] ) > tol )
        {
            c++;
        }
    }
    if ( c != 0)
    {
        cout << "Eigenvalue test did not PASS :( " << " --> Is rho_N ≠ 0 ?" << endl;

    } else {
        cout << "Eigenvalue test has PASSED :) " << endl;

    }
}
