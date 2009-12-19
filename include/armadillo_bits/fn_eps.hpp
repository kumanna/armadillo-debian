// Copyright (C) 2009 NICTA
// Copyright (C) 2009 Dimitrios Bouzas
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup fn_eps
//! @{



//! \brief
//! eps version for non-complex matrices and vectors
template<typename T1>
inline
const Op<T1, op_eps>
eps(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  isnt_fltpt<eT>::check();
  
  return Op<T1, op_eps>(X.get_ref());
  }



//! \brief
//! eps version for complex matrices and vectors
template<typename T1>
inline
Mat< typename T1::pod_type >
eps(const Base< std::complex<typename T1::pod_type>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type   T;
  typedef typename T1::elem_type eT;
  
  isnt_fltpt<T>::check();
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  Mat<T> out;
  op_eps::direct_eps(out, A);
  
  return out;
  }



//! \brief
//! eps version for scalars of type float
arma_inline
float
eps(const float x)
  {
  return op_eps::direct_eps(x);
  }



//! \brief
//! eps version for scalars of type double
arma_inline
double
eps(const double x)
  {
  return op_eps::direct_eps(x);
  }



//! \brief
//! eps version for std::complex<float>
arma_inline
float
eps(const std::complex<float>& x)
  {
  return op_eps::direct_eps(x);
  }



//! \brief
//! eps version for std::complex<double>
arma_inline
double
eps(const std::complex<double>& x)
  {
  return op_eps::direct_eps(x);
  }



//! @}
