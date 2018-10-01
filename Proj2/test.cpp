#include <armadillo>
#include "jacobi_rotate.h"
#include "maxoffdiag.h"
#include "analytic_eigvals.h"

using namespace std;
using namespace arma;


void test_eigvals(int N, vec &lambdas, vec &a_lambdas, double a, double d)
{
    cout << "TEST_EIGVALS RUNS ...  ";
    for ( int i = 0; i < N; i++)
    {
         a_lambdas[i] = analytic_eigenvals(a,d,N,i+1);
    }
    int c = 0;
    double tol = 1e-6;
    for ( int i = 0; i < N; i++)
    {
        if ( fabs( lambdas[i] - a_lambdas[i] ) > tol )
        {
            c+=1;
        }
    }
    if ( c != 0)
    {
        cout << "ERROR --> Eigenvalue test did not PASS :( " << endl;
    } else {
        cout << "SUCCESS! --> Eigenvalue test has PASSED :) " << endl;
    }
}

void test_eigvals_quantum( vec &lambdas, double p_N)
{
    {
        cout << "TEST_EIGVALS (QUANTUM) RUNS ...  ";
        int c = 0;
        double tol = 1e-3;
        if ( fabs( lambdas[0] - 3 ) > tol &&
             fabs( lambdas[1] - 7 ) > tol &&
             fabs( lambdas[2] - 11 ) > tol &&
             fabs( lambdas[4] - 15 ) > tol)
        {
            c++;
        }
        if ( c > 0 )
        {
            cout << "ERROR --> The 4 first eigvals are not = 3,7,11,15 within a tolerance of " << tol << ", try larger N " << endl;
        } else {
            cout << "SUCCESS! --> Eigenvector test for the quantum case has passed :) " << endl;
        }
    }

}

void test_ortho( mat R, int N )
{
    double s;
    int c = 0;
    cout << "TEST_ORTHO RUNS ...  ";
    for ( int k = 0; k < N; k++)
    {
        s = 0;
        for ( int j = k; j < N - 1; j++)
        {
            for ( int i = 0; i < N; i++)
            {
                s += R(i,j)*R(i,j+1);
            }
        }
        double tol = 1e-10;
        if ( fabs(s) > tol )
        {
            c++;
        }
    }
    if ( c > 0)
    {
        cout << "ERROR --> Eigenvector test did not PASS :( ... dotprods \neq zero .... " << endl;
    } else {
        cout << "SUCCESS! --> Eigenvector test has passed :) " << endl;
    }

}
