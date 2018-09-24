#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
//#include "catch.hpp"
#include "jacobi_rotate.h"
#include "maxoffdiag.h"
#include "test_ortho.h"
#include "test_eigvals.h"
#include "analytic_eigvals.h"
#include "jacobi_method.h"
#include "armadillo_eigpair.h"

using namespace std;
using namespace arma;

int main(int argc, const char * argv[])
{
    int N = 200;
    double p_N = 5;
    double w = 0.0;
    double w2 = w*w;
    int k, l;
    double h, Diag, NonDiag;

    vec p = zeros<vec>(N);
    vec V = zeros<vec>(N);

    if ( p_N != 0)
    {
        h = p_N / double(N + 1);
        Diag = 2.0/(h*h);
        NonDiag = -1.0/(h*h);
        for ( int i = 0; i < N; i++)
        {
            p[i] = (i + 1)*h;
            V[i] = p[i]*p[i];
//            cout << "V[i] = " << V[i] << endl;
        }
        if ( w > 0)
        {
            for ( int i = 0; i < N; i++ )
            {
                V[i] = V[i]*w2 + 1/p[i];
            }
        }
    } else {
        Diag = 2.0;
        NonDiag = -1.0;
    }

    // Create A and R
    mat R(N,N);
    mat A(N,N);
    for ( int i = 0; i < N; i++)
    {
        for ( int j = 0; j < N; j++)
        {
            if ( i == j )
            {
                A(i,j) = Diag + V[i];
                R(i,j) = 1.0;
            }
            else if ( i == j - 1 )
            {
                A(i,j) = NonDiag;
                R(i,j) = 0.0;
            }
            else if ( i == j + 1 )
            {
                A(i,j) = NonDiag;
                R(i,j) = 0.0;
            }
        }
    }

    vec lambdas(N);
    vec a_lambdas(N);
    test_ortho( R, N );
    test_eigvals( N, lambdas, a_lambdas );

    jacobi_method( N, A, R );
    armadillo_eigpair( N, A );


    return 0;

}

