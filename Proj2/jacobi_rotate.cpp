//#include <iostream>
//#include <fstream>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;


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


