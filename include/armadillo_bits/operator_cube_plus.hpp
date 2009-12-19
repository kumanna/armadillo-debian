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


//! \addtogroup operator_cube_plus
//! @{



//! unary plus operation (does nothing, but is required for completeness)
template<typename T1>
arma_inline
const BaseCube<typename T1::elem_type,T1>&
operator+
(const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X;
  }



//! Base + scalar
template<typename T1>
arma_inline
const OpCube<T1, op_scalar_plus>
operator+
(const BaseCube<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_plus>(X.get_ref(), k);
  }



//! op + scalar, level 2
template<typename T1>
arma_inline
const OpCube<T1,op_scalar_plus>
operator+
(const OpCube<T1,op_scalar_plus>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_plus>(X.m, X.aux + k);
  }



//! scalar + Base
template<typename T1>
arma_inline
const OpCube<T1, op_scalar_plus>
operator+
(const typename T1::elem_type k, const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_plus>(X.get_ref(), k);  // NOTE: order is swapped
  }



//
// addition of Base objects with different element types
//



//! Base + Base
template<typename eT1, typename T1, typename eT2, typename T2>
arma_inline
Cube<typename promote_type<eT1,eT2>::result>
operator+
(const BaseCube<eT1,T1>& X, const BaseCube<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();

  const unwrap_cube<T1> tmp1(X.get_ref());
  const unwrap_cube<T2> tmp2(Y.get_ref());
  
  const Cube<eT1>& A = tmp1.M;
  const Cube<eT2>& B = tmp2.M;
  
  Cube< typename promote_type<eT1,eT2>::result > out;
  glue_cube_plus::apply_mixed(out, A, B);
  
  return out;
  }



//
// addition of Base objects with same element types
//



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<std::complex<double>,T1>& X, const BaseCube<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<std::complex<float>,T1>& X, const BaseCube<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<double,T1>& X, const BaseCube<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<float,T1>& X, const BaseCube<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<s32,T1>& X, const BaseCube<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<u32,T1>& X, const BaseCube<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<s16,T1>& X, const BaseCube<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<u16,T1>& X, const BaseCube<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<s8,T1>& X, const BaseCube<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const GlueCube<T1, T2, glue_cube_plus>
operator+
(const BaseCube<u8,T1>& X, const BaseCube<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return GlueCube<T1, T2, glue_cube_plus>(X.get_ref(), Y.get_ref());
  }



//! @}
