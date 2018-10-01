#include <armadillo>

using namespace arma;


void armadillo_eigpair( int N, mat A)
{
//  FIND EIGENVALUES AND EIGENVECTORS USING ARMADILLO
    clock_t start,end;
    start = clock();

    vec eigval(N);
    mat eigvec(N,N);
    eig_sym(eigval, eigvec, A);
//    cout << "Armadillo eigvals " << endl << eigval << endl;
//    cout << "Armadillo eigvecs " << endl << eigvec << endl;
    end = clock();
    double t_s = (double(end - start)/CLOCKS_PER_SEC)*1000;

    cout << "Time to run Armadillo; " << t_s << "ms, when N = " << N << endl;
}

