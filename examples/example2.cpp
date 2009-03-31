#include <iostream>

#include "armadillo"

using namespace arma;
using namespace std;


int main(int argc, char** argv)
  {
  // matrix construction from a string representation
  mat A = \
   "\
   0.555950  0.274690  0.540605  0.798938;\
   0.108929  0.830123  0.891726  0.895283;\
   0.948014  0.973234  0.216504  0.883152;\
   0.023787  0.675382  0.231751  0.450332;\
   ";

  // print to the cout stream
  // with an optional string before the contents of the matrix
  A.print("A =");

  // determinant
  cout << "det(A) = " << det(A) << endl;

  // inverse
  cout << "inv(A) = " << endl << inv(A) << endl;


  return 0;
  }

