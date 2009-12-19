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


//! \addtogroup operator_cube_minus
//! @{



//! unary -
template<typename T1>
arma_inline
const OpCube<T1, op_neg>
operator-
(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1,op_neg>(X.get_ref());
  }



//! cancellation of two consecutive negations: -(-T1)
template<typename T1>
arma_inline
const T1&
operator-
(const OpCube<T1, op_neg>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.m;
  }



//! Base - scalar
template<typename T1>
arma_inline
const OpCube<T1, op_scalar_minus_post>
operator-
(const BaseCube<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_minus_post>(X.get_ref(), k);
  }



//! scalar - Base
template<typename T1>
arma_inline
const OpCube<T1, op_scalar_minus_pre>
operator-
(const typename T1::elem_type k, const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_minus_pre>(X.get_ref(), k);
  }



//! Base - - Base = Base + Base
template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator-
  (
  const BaseCube<typename T1::elem_type, T1                >& X,
  const BaseCube<typename T1::elem_type, OpCube<T2,op_neg> >& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const T1& A = X.get_ref();
  const T2& B = (Y.get_ref()).m;
  
  return GlueCube<T1, T2, glue_cube_plus>(A,B);
  }



//! Base - OpCube<T2,op_neg> = Base + T2
template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator-
(const BaseCube<typename T2::elem_type,T1>& X, const OpCube<T2, op_neg>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.m);
  }



//
// subtraction of Base objects with different element types
//



//! Base - Base
template<typename eT1, typename T1, typename eT2, typename T2>
arma_inline
Cube<typename promote_type<eT1,eT2>::result>
operator-
(const BaseCube<eT1,T1>& X, const BaseCube<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();
  
  const unwrap_cube<T1> tmp1(X.get_ref());
  const unwrap_cube<T2> tmp2(Y.get_ref());
  
  const Cube<eT1>& A = tmp1.M;
  const Cube<eT2>& B = tmp2.M;
  
  Cube< typename promote_type<eT1,eT2>::result > out;
  glue_cube_minus::apply_mixed(out, A, B);
  
  return out;
  }



//
// subtraction of Base objects with same element types
//



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<std::complex<double>,T1>& X, const BaseCube<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<std::complex<float>,T1>& X, const BaseCube<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<double,T1>& X, const BaseCube<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<float,T1>& X, const BaseCube<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<s32,T1>& X, const BaseCube<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<u32,T1>& X, const BaseCube<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<s16,T1>& X, const BaseCube<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<u16,T1>& X, const BaseCube<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<s8,T1>& X, const BaseCube<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_minus>
operator-
(const BaseCube<u8,T1>& X, const BaseCube<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_minus>(X.get_ref(), Y.get_ref());
  }



//! @}
