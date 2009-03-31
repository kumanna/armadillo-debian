// Copyright (C) 2009 NICTA
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup auxlib
//! @{


//! wrapper for accessing external functions defined in ATLAS, LAPACK or BLAS libraries
class auxlib
  {
  public:
  
  template<typename eT>
  inline static const eT& tmp_real(const eT& X)              { return X; }
  
  template<typename T>
  inline static const T&  tmp_real(const std::complex<T>& X) { return X.real(); }


  //
  // inv
  
  template<typename eT>
  inline static void inv_noalias(basic_mat<eT>& out, const basic_mat<eT>& X);
  
  template<typename eT>
  inline static void inv_inplace(basic_mat<eT>& X);
  
  
  //
  // det
  
  template<typename eT>
  inline static eT det(const basic_mat<eT>& X);
  
  
  //
  // lu
  
  template<typename eT>
  inline static void lu(basic_mat<eT>& L, basic_mat<eT>& U, podarray<int>& ipiv, const basic_mat<eT>& X_orig);
  
  template<typename eT>
  inline static void lu(basic_mat<eT>& L, basic_mat<eT>& U, basic_mat<eT>& P, const basic_mat<eT>& X);
  
  template<typename eT>
  inline static void lu(basic_mat<eT>& L, basic_mat<eT>& U, const basic_mat<eT>& X);
  
  
  //
  // eig
  
  template<typename eT> 
  inline static void eig(basic_colvec<eT>& eigval, const basic_mat<eT>& A);
  
  template<typename eT>
  inline static void eig(basic_colvec<eT>& eigval, basic_mat<eT>& eigvec, const basic_mat<eT>& A);
  
  
  //
  // chol
  
  template<typename eT> 
  inline static bool chol(basic_mat<eT>& out, const basic_mat<eT>& X);
  
  
  //
  // qr
  
  template<typename eT> 
  inline static bool qr(basic_mat<eT>& Q, basic_mat<eT>& R, const basic_mat<eT>& X);
  
  
  //
  // svd
  
  template<typename eT>
  inline static bool svd(basic_colvec<eT>& S, const basic_mat<eT>& X);
  
  template<typename T>
  inline static bool svd(basic_colvec<T>& S, const basic_mat< std::complex<T> >& X);
  
  template<typename eT>
  inline static bool svd(basic_mat<eT>& U, basic_colvec<eT>& S, basic_mat<eT>& V, const basic_mat<eT>& X);
  
  template<typename T> 
  inline static bool svd(basic_mat< std::complex<T> >& U, basic_colvec<T>& S, basic_mat< std::complex<T> >& V, const basic_mat< std::complex<T> >& X);
  
  
  // solve
  
  template<typename eT>
  inline static bool solve(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B);
  
  template<typename eT>
  inline static bool solve_od(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B);
  
  template<typename eT>
  inline static bool solve_ud(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B);
  
  };


//! @}
