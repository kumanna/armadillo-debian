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


//! \addtogroup operator_div
//! @{



//! Base / scalar
template<typename T1>
arma_inline
const Op<T1, op_scalar_div_post>
operator/
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_div_post>(X.get_ref(), k);
  }



//! Op<mat,op_ones_full> / scalar
template<typename eT>
arma_inline
Mat<eT>
operator/
(const Op<Mat<eT>,op_ones_full>& X, const eT k)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> tmp(X.aux_u32_a, X.aux_u32_b);
  tmp.fill( eT(1)/k );
  
  return tmp;
  }



//! Op<mat,op_ones_diag> / scalar
template<typename eT>
arma_inline
Mat<eT>
operator/
(const Op<Mat<eT>,op_ones_diag>& X, const eT k)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.aux_u32_a, X.aux_u32_b);
  
  const eT inv_k = eT(1)/k;
  
  for(u32 i=0; i<out.n_rows; ++i)
    {
    out.at(i,i) = inv_k;
    }
  
  return out;
  }



//! scalar / Base
template<typename T1>
arma_inline
const Op<T1, op_scalar_div_pre>
operator/
(const typename T1::elem_type k, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_div_pre>(X.get_ref(), k);
  }



//
// element-wise division of Base objects with different element types
//



//! Base / Base
template<typename eT1, typename T1, typename eT2, typename T2>
arma_inline
Mat<typename promote_type<eT1,eT2>::result>
operator/
(const Base<eT1,T1>& X, const Base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT1>& A = tmp1.M;
  const Mat<eT2>& B = tmp2.M;
  
  Mat< typename promote_type<eT1,eT2>::result > out;
  glue_div::apply_mixed(out, A, B);
  
  return out;
  }



//
// element-wise division of Base objects with same element types
//



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<std::complex<double>,T1>& X, const Base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<std::complex<float>,T1>& X, const Base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<double,T1>& X, const Base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<float,T1>& X, const Base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<s32,T1>& X, const Base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<u32,T1>& X, const Base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<s16,T1>& X, const Base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<u16,T1>& X, const Base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<s8,T1>& X, const Base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_div>
operator/
(const Base<u8,T1>& X, const Base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_div>(X.get_ref(), Y.get_ref());
  }



//! @}
