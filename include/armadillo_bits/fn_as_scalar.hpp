// Copyright (C) 2010 NICTA and the authors listed below
// http://nicta.com.au
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


//! \addtogroup fn_as_scalar
//! @{



template<typename T1>
inline
typename T1::elem_type
as_scalar(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( (A.n_elem != 1), "as_scalar(): expression doesn't evaluate to exactly one element" );
  
  return A.mem[0];
  }



template<typename T1>
inline
typename T1::elem_type
as_scalar(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> tmp(X.get_ref());
  const Cube<eT>& A   = tmp.M;
  
  arma_debug_check( (A.n_elem != 1), "as_scalar(): expression doesn't evaluate to exactly one element" );
  
  return A.mem[0];
  }


template<typename T> struct arma_scalar_only { };

template<> struct arma_scalar_only<char>   { typedef char   result; };
template<> struct arma_scalar_only<short>  { typedef short  result; };
template<> struct arma_scalar_only<int>    { typedef int    result; };
template<> struct arma_scalar_only<long>   { typedef long   result; };
template<> struct arma_scalar_only<float>  { typedef float  result; };
template<> struct arma_scalar_only<double> { typedef double result; };

template<> struct arma_scalar_only<unsigned char>  { typedef unsigned char  result; };
template<> struct arma_scalar_only<unsigned short> { typedef unsigned short result; };
template<> struct arma_scalar_only<unsigned int>   { typedef unsigned int   result; };
template<> struct arma_scalar_only<unsigned long>  { typedef unsigned long  result; };

template<typename T> struct arma_scalar_only< std::complex<T> > { typedef std::complex<T> result; };


template<typename T>
arma_inline
const typename arma_scalar_only<T>::result &
as_scalar(const T& x)
  {
  return x;
  }



//! @}
