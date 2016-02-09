#include <iostream>
#include "armadillo"

using namespace arma;
using namespace std;

int main(int argc, char *argv[])
{
    /*
     Program to be used to solve the linear algebra equation
     A[n,n] * x[n] = w[n]
     when A is a tridiagonal matrix.
    */

    mat A; mat L; mat U;
    /*initialize Matrix to be solved*/
    A << 2 << 1 << 0 << 0 << endr
      << 1 << 2 << 1 << 0 << endr
      << 0 << 1 << 2 << 1 << endr
      << 0 << 0 << 1 << 2 << endr;
    /*initialize vector of unknown variables*/
    vec x(4); x.zeros(4);

    /*initialize right side of equation*/
    vec w(4); w.randu(4);

    /* performs L-U decomposition of matrix A*/
    lu(L, U, A);

    /*determine how many times the loop will run*/
    int n = x.n_rows;

    /*initialize place holder (y-vector) for Lower Matrix forward substitution*/
    vec y(n);
    y(1) = w(1);

    /*fill in the rest of the y-vector*/
    for (int i = 2; i <= n; i++)
    {
        y(i) = w(i) - L(i,i-1)*y(i-1);
    }

    /*fill in x-vector*/
    x(n) = y(n);

    for (int j = n-1;j <= 1; j--)
    {
         x(j) = (y(j)-(U(j,j+1)*x(j+1)))/U(j,j);
    }

    A.print("Matrix A");
    w.print("Vector w");
    x.print("Answer x:");
    return 0;
}
