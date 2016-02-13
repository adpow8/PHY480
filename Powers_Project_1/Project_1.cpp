#include <iostream>
#include "armadillo"

using namespace arma;
using namespace std;

double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

double f(double x) {return 100*exp(-10*x);}

int main(int argc, char *argv[])
{
    /*
     Program to be used to solve the linear algebra equation
     A[n,n] * x[n] = w[n]
     when A is a tridiagonal matrix.
    */

    // Declaration of initial variables:
    char *outfilename;
    int n;

    // Read in output file and n,
    // abort if there are too few command-line arguments:
    if( argc <= 2 ){
      cout << "Bad Usage: " << argv[0] <<
          " read also output file and n (int) on same line" << endl;
      exit(1);
    }
    else{
      outfilename = argv[1]; // first command line argument.
      n = atoi(argv[2])-1;  // second command line argument.
    }

    //Vectors making up the tridiagonal matrix A
    //Zeroth element included to make indexing easy
    vec a(n+1); vec b(n+1); vec c(n+1);

    //Temporary element needed for Gaussian elimination
    vec dTemp(n+1);

    //Constants within the problem
    double h = 1.0/(n+1.0);
    vec x(n+2)

    //initialize vector components of matrix to be solved
    for (int i = 1; i<= n; i++)
    {
        a(i) = 1;
        b(i) = 2;
        c(i) = 1;
    }

    a(1) = 0;
    c(n) = 0;

    //initialize right side of equation
    vec w(n+1);

    for (int i = 1; i <= n; i++)
    {
        w(i) = _________________;

    }


    //initialize place holder (y-vector) for Lower Matrix forward substitution
    vec y(n+1);
    y(0) = 0;
    w(0) = 0;

    //fill in the rest of the y-vector
    for (int i = 1; i <= n; i++)
    {
        y(i) = w(i) - a(i)*y(i);
    }

    //fill in x-vector
    x(n) = y(n);

    for (int i = n; i >= 1; i--)
    {
         x(i) = ( y(i) - c(n)*x(n+1) ) / b(n);
    }

    w.print("Vector w: ");
    x.print("Answer x:");
    u.print("Answer u: ")

    return 0;


}
