//#include <iostream>
//#include <cmath>
#include <armadillo>
//#include "catch.hpp"
#include "jacobi_rotate.h"
#include "maxoffdiag.h"

using namespace std;
using namespace arma;

// TEST ORTHOGONALITY OF EIGENVECTORS
//void test_ortho( mat R, int N )
//{
//    double s;
//    int c = 0;
//    cout << "TEST_ORTHO RUNS ...  ";
//    for ( int k = 0; k < N; k++)
//    {
//        s = 0;
//        for ( int j = k; j < N - 1; j++)
//        {
//            for ( int i = 0; i < N; i++)
//            {
//                s += R(i,j)*R(i,j+1);
//            }
//        }
//        double tol = 1e-10;
//        if ( fabs(s) > tol )
//        {
//            c++;
//        }
//    }
//    if ( c > 0)
//    {
//        cout << "Eigenvector test did not PASS :( " << endl;
//    } else {
//        cout << "Eigenvector test has passed :)" << endl;
//    }

//}
