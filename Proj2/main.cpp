#include <iostream>
#include <fstream>
#include <cmath>
#include <armadillo>
//#include "catch.hpp"
#include "jacobi_rotate.h"
#include "maxoffdiag.h"
#include "test.h"
#include "analytic_eigvals.h"
#include "jacobi_method.h"
#include "armadillo_eigpair.h"

using namespace std;
using namespace arma;

int main(int argc, const char * argv[])
{
    int N = 300;
    double p_N = 15;
    double w = 0.25;
    double w2 = w*w;
    int k, l;
    double h, d, a;

    vec p = zeros<vec>(N);
    vec V = zeros<vec>(N);

    if ( p_N != 0)
    {
        h = p_N / double(N + 1);
        d = 2.0/(h*h);
        a = -1.0/(h*h);
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
        d = 2.0;
        a = -1.0;
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
                A(i,j) = d + V[i];
                R(i,j) = 1.0;
            }
            else if ( i == j - 1 )
            {
                A(i,j) = a;
                R(i,j) = 0.0;
            }
            else if ( i == j + 1 )
            {
                A(i,j) = a;
                R(i,j) = 0.0;
            }
        }
    }

    vec lambdas(N);
    vec a_lambdas(N);
    jacobi_method( N, A, R, lambdas );
    armadillo_eigpair( N, A );

    // TESTS
    if ( p_N == 0)
    {
        test_eigvals( N, lambdas, a_lambdas, a, d );
    }

    if ( p_N != 0 && w == 0.0){
        test_eigvals_quantum( lambdas, p_N);
    }
    test_ortho( R, N );

    cout << "using N = " << N << " and p_N = " << p_N << " and w = " << w <<  " the four first eigvals are " << endl;
    for ( int i = 0; i < 4; i++ )
    {
        cout << lambdas[i] << endl;
    }


    return 0;
}

