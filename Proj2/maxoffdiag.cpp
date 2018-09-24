//#include <iostream>
//#include <fstream>
//#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

double maxoffdiag ( mat A, int *k, int *l, int N)
{
    double max = 0;
    double a_ij;
    for ( int i = 0; i < N; i++ )
    {
        for ( int j = i+1; j < N; j++ )
        {
            a_ij = fabs( A(i,j) );
            if ( a_ij > max)
            {
                max = a_ij;
//                cout << "MAX = " << max << endl;
                *k = i;
                *l = j;
            }
        }
    }
    return max;
}
