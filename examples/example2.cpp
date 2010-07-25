#include <iostream>

#include "armadillo"

using namespace arma;
using namespace std;


int main(int argc, char** argv)
  {
  cout << arma_version::as_string() << endl;
  
  // matrix construction from a string representation
  mat A = \
   "\
   0.165300  0.454037  0.995795  0.124098  0.047084;\
   0.688782  0.036549  0.552848  0.937664  0.866401;\
   0.348740  0.479388  0.506228  0.145673  0.491547;\
   0.148678  0.682258  0.571154  0.874724  0.444632;\
   0.245726  0.595218  0.409327  0.367827  0.385736;\
   ";
  
  // print to the cout stream
  // with an optional string before the contents of the matrix
  A.print("A =");
  
  // determinant
  cout << "det(A) = " << det(A) << endl;
  
  // inverse
  cout << "inv(A) = " << endl << inv(A) << endl;
  
  
  //
  
  double k = 1.23;
  
  mat    B = randu<mat>(5,5);
  mat    C = randu<mat>(5,5);
  
  rowvec r = randu<rowvec>(5);
  colvec q = randu<colvec>(5);
  
  
  // examples of some expressions
  // for which optimised implementations exist
  
  // optimised implementation of a trinary expression
  // that results in a scalar
  cout << "as_scalar( r*inv(diagmat(B))*q ) = ";
  cout << as_scalar( r*inv(diagmat(B))*q ) << endl;
  
  // example of an expression which is optimised 
  // as a call to the dgemm() function in BLAS:
  cout << "k*trans(B)*C = " << endl << k*trans(B)*C;
  
  
  // If you want to see a trace of how Armadillo
  // evaluates expressions, compile with the
  // ARMA_EXTRA_DEBUG macro defined.
  // This was designed to work the the GCC compiler,
  // but it may also work with other compilers
  // if you have the Boost libraries installed
  // and told Armadillo to use them.
  // 
  // Example for GCC:
  // g++ -DARMA_EXTRA_DEBUG -o example2 example2.cpp -larmadillo
  // 
  // Running example2 will now produce a truckload of messages,
  // so you may want to redirect the output to a log file.
  
  return 0;
  }

