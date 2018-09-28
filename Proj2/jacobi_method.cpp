#include <armadillo>
#include "maxoffdiag.h"
#include "jacobi_rotate.h"

using namespace arma;

void jacobi_method( int N, mat &A, mat &R, vec &lambdas)
{
    int k,l;
    double epsilon = 1e-10;
    double max_num_of_itera = (double) N* (double) N* (double) N;
    int iterations = 0;
    double max_off_diag = maxoffdiag ( A , &k, &l, N);

//  FIND EIGENVALUES AND EIGENVECTORS USING THE JACOBI_METHOD
    clock_t start_,end_;
    start_ = clock();
    while ( fabs( max_off_diag ) > epsilon && (double) iterations < max_num_of_itera)
    {
        max_off_diag = maxoffdiag ( A , &k, &l, N);
        Jacobi_rotate( A, R, k, l, N);
        iterations++;
    }
    end_ = clock();
    double t_s_ = (double(end_ - start_)/CLOCKS_PER_SEC)*1000;

    end_ = clock();
    vec diags(N);
    diags = A.diag();
    lambdas = sort(diags);

//    cout << "Jacobi_rotate eigvals (cout from jacobi_method.cpp) " << endl;
//    for ( int i = 0; i < 4; i++)
//    {
//        cout << lambdas[i] << endl;
//    }



//    cout << "Jacobi_rotate eigvecs; R = " << endl << R << endl << endl;
    cout << "Number of iterations, or similarity transformations = " << iterations << endl << endl;
    cout << "Time to run Jacobi_rotate; " << t_s_ << "ms when N = " << N << endl<<endl;


}

