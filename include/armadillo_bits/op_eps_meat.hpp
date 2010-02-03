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



//! \addtogroup op_eps
//! @{



template<typename eT>
inline
eT
op_eps::direct_eps(const eT x)
  {
  //arma_extra_debug_sigprint();
  
  // acording to IEEE Standard for Floating-Point Arithmetic (IEEE 754)
  // the mantissa length for double is 53 bits = std::numeric_limits<double>::digits
  // the mantissa length for float  is 24 bits = std::numeric_limits<float >::digits
  
  //return std::pow( std::numeric_limits<eT>::radix, (std::floor(std::log10(std::abs(x))/std::log10(std::numeric_limits<eT>::radix))-(std::numeric_limits<eT>::digits-1)) );
  
  const eT radix_eT     = eT(std::numeric_limits<eT>::radix);
  const eT digits_m1_eT = eT(std::numeric_limits<eT>::digits - 1);
  
  return std::pow( radix_eT, eT(std::floor(std::log10(std::abs(x))/std::log10(radix_eT)) - digits_m1_eT) );
  }



template<typename T>
inline
T
op_eps::direct_eps(const std::complex<T>& x)
  {
  //arma_extra_debug_sigprint();
  
  //return std::pow( std::numeric_limits<T>::radix, (std::floor(std::log10(std::abs(x))/std::log10(std::numeric_limits<T>::radix))-(std::numeric_limits<T>::digits-1)) );
  
  const T radix_T     = T(std::numeric_limits<T>::radix);
  const T digits_m1_T = T(std::numeric_limits<T>::digits - 1);
  
  return std::pow( radix_T, T(std::floor(std::log10(std::abs(x))/std::log10(radix_T)) - digits_m1_T) );
  }



template<typename eT>
inline
void
op_eps::direct_eps(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  out.copy_size(A);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = op_eps::direct_eps( A_mem[i] );
    }
  
  }



template<typename T>
inline
void
op_eps::direct_eps(Mat<T>& out, const Mat< std::complex<T> >& A)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  out.copy_size(A);
  
         T* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = op_eps::direct_eps( A_mem[i] );
    }
  
  }



template<typename T1>
inline
void
op_eps::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_eps>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type  eT;
  
  // in this case it doesn't matter if there is aliasing
  // (i.e. &out = &in.m), hence we use unwrap rather than unwrap_check
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  op_eps::direct_eps(out, A);
  }



//! @}
