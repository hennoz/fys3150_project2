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

int main(int argc, const char * argv[])
{
    int N = atoi(argv[1]);
    int l;
    int k;
    mat A(N,N); //  Now to fill matrix A with 2's and -1's
    for ( int i = 0; i < N; i++)
    {
        for ( int j = 0; j < N; j++)
        {
            if ( i == j)
            {
                A(i,j) = 2;
            }
            else if ( i == j - 1)
            {
                A(i,j) = - 1;
            }
            else if ( i == j + 1) {
                A(i,j) = - 1;
            }
        }
    }


//  FIND EIGENVALUES AND EIGENVECTORS USING ARMADILLO
    clock_t start,end;
    start = clock();

    vec eigval(N);
    mat eigvec(N,N);
    eig_sym(eigval, eigvec, A);
    cout << "Armadillo eigenvals " << endl << eigval << endl;
    cout << "Armadillo eigenvecs " << endl << eigvec << endl;
    end = clock();
    double t_s = (double(end - start)/CLOCKS_PER_SEC)*1000;

    //  Eigenvector matrix R
    mat R(N,N);
    for ( int i = 0; i < N; i++)
    {
        for ( int j = 0; j < N; j++)
        {
            if ( i == j)
            {
                R(i,j) = 1.0;
            } else {
                R(i,j) = 0.0;
            }
        }
    }

    double epsilon = 1e-13;
    double max_num_of_itera = (double) N* (double) N* (double) N;
    int iterations = 0;
    double max_off_diag = maxoffdiag ( A , &k, &l, N);


//  FIND EIGENVECTORS AND EIGENVECTORS USING JACOBI_ROTATE ALGORITHM
    clock_t start_,end_;
    start_ = clock();
    while ( fabs(max_off_diag) > epsilon && (double) iterations < max_num_of_itera)
    {
        max_off_diag = maxoffdiag ( A , &k, &l, N);
        Jacobi_rotate( A, R, k, l, N);
        iterations++;
    }
    end_ = clock();
    double t_s_ = (double(end_ - start_)/CLOCKS_PER_SEC)*1000;

    vec diags(N);
    diags = A.diag();
    vec lambdas = sort(diags);
//    cout << "Jacobi_rotate eigenvals " << endl << lambdas << endl;
    cout << "Jacobi_rotate eigenvecs " << endl << R << endl << endl;

    // Analytic solutions for eigenvals
    cout << "Analytical eigenvals = " << endl;
    double *a_lambdas = new double[N];
    for ( int i = 0; i < N; i++)
    {
         a_lambdas[i] = analytic_eigenvals(-1,2,N,i+1);
         cout << "   " << a_lambdas[i] << endl;
    }
    cout << endl;
    cout << "Time to run Armadillo; " << t_s << "ms, when N = " << N << endl<<endl;
    cout << "Time to run Jacobi_rotate; " << t_s_ << "ms when N = " << N << endl<<endl;
    cout << "Number of iterations, or similarity transformations = " << max_num_of_itera << endl << endl;


    //  TEST EIGENVALUES (STOPS PROGRAM IF NOT PASSED)
    double tol = 1e-10;
    for ( int i = 0; i < N; i++)
    {
        if ( fabs( lambdas[i] - a_lambdas[i] ) > tol )
        {
            cout << "WARNING!! Jacobi_rotate broke; lambda_" << i+1 << " does NOT correspond to the analytical solution" << endl;
            return 0;
        }
    }
    cout << "Eigenvalue test has PASSED :) " << endl;

    // TEST ORTHOGONALITY OF EIGENVECTORS (STOPS PROGRAM IF NOT PASSED)
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
//        cout << "dotproduct of the eigenvectors = " << s << endl;
        if ( fabs(s) > tol )
        {
            cout << "WARNING!! Jacobi_rotate broke; Eigenvector dotproduct ≠ 0" << endl;
            return 0;
        }

    }
    cout << "Eigenvector test has PASSED :) " << endl;


    delete [] a_lambdas;
    return 0;
}

double maxoffdiag ( mat A, int *k, int *l, int N)
{
    double max = 0;
    for ( int i = 0; i < N; i++ )
    {
        for ( int j = i+1; j < N; j++ )
        {
            if ( fabs( A(i,j) ) > max)
            {
                max = fabs( A(i,j) );
                cout << "MAX = " << max << endl;
                *k = i;
                *l = j;
            }
        }
    }

    return max;
}

void Jacobi_rotate( mat &A, mat &R, int k, int l, int N)
{
    double s,c, tau, t;
    //  c, s = cos\theta, sin\theta
    //  if a_kl ≠ 0, tau = (a_ll - a_kk)/(2*a_kl)
    if ( A(k,l) != 0.0)
    {
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        /*
         * t = (±sqrt(1 + tau^2) - tau) *[(sqrt(1 + tau^2) - tau) / (sqrt(1 + tau^2) - tau)]
         *   = ±1.0 / (sqrt(1 + tau^2) - tau)
         * => times conjugate up and down to prevent loss of data
         */

        //  if tau is positive
        if ( tau >= 0.0 )
        {
            t = 1.0 / (sqrt(1 + tau*tau) + tau);
        }
        //  if tau is negative
        else
        {
            t = -1.0 / (sqrt(1 + tau*tau) - tau);
        }
        c = 1.0 / sqrt(1 + t*t);
        s = c*t;
    //  Define c=1 and s=0 if a_kl = 0 (otherwise will cause division by zero in tau = (a_ll - a_kk)/(2a_kl))
    } else {
    c = 1.0;
    s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    //  Making values for the diagonal elements
    a_kk = A(k,k);
    a_ll = A(l,l);
    //  Jacobi's method
    //  Fill matrix A's diagonal elements and 0's for the rest
    A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
    A(k,l) = 0.0;
    A(l,k) = 0.0; //  Hardcode zeros to save computing time
    // Iteration to fill rest of matrix A (s216 lect notes)

    for ( int i = 0; i < N; i++ )
    {
        if ( i != k && i != l )
        {
            a_ik   = A(i,k);
            a_il   = A(i,l);
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = A(i,l);
        }
        //  Fill eigenvector matrix R
        r_ik = R(i,k);
        r_il = R(i,l);

        //  S^TAS
        R(i,k) = r_ik*c - r_il*s;
        R(i,l) = r_il*c + r_ik*s;
    }
}

