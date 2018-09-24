#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <algorithm>
#include <sys/ioctl.h>
#include <time.h>
#include <unistd.h>
#include <ctime>
#include "catch.hpp"
#include "Tridiag.h"
#include "jacobi.h"


//  TEST EIGENVALUES
double tol = 1e-10;
for ( int i = 0; i < N; i++)
{
    if ( fabs( lambdas[i] - a_lambdas[i] ) > tol )
    {
//            cout << "WARNING!! Jacobi_rotate broke; lambda_" << i+1 << " does NOT correspond to the analytical solution" << endl;
        cout << "Eigenvalue test did not PASS :( " << endl;
    } else {
        cout << "Eigenvalue test has PASSED :) " << endl;

    }
}

// TEST ORTHOGONALITY OF EIGENVECTORS
for ( int k = 0; k < N; k++)
{
    double s = 0;
    for ( int j = k; j < N-1; j++)
    {
        for ( int i = 0; i < N; i++)
        {
            s += R(i,j)*R(i,j+1);
        }
    }

    if ( fabs(s) > tol )
    {
//            cout << "WARNING!! Jacobi_rotate broke; Eigenvector dotproduct ≠ 0" << endl;
        cout << "Eigenvector test did not PASS :( " << endl;
    } else {
        cout << "Eigenvector test has PASSED :) " << endl;

    }

}
