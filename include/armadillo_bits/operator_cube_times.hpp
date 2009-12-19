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




//! \addtogroup operator_cube_times
//! @{



//! BaseCube * scalar
template<typename T1>
arma_inline
const OpCube<T1, op_scalar_times>
operator*
(const BaseCube<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_times>(X.get_ref(),k);
  }



//! op * scalar, level 2
template<typename T1>
arma_inline
const OpCube<T1,op_scalar_times>
operator*
(const OpCube<T1,op_scalar_times>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_times>(X.m, X.aux * k);
  }



//! OpCube<cube,op_ones_full> * scalar
template<typename eT>
arma_inline
Cube<eT>
operator*
(const OpCube<Cube<eT>,op_ones_full>& X, const eT k)
  {
  arma_extra_debug_sigprint();
    
  Cube<eT> tmp(X.aux_u32_a, X.aux_u32_b, X.aux_u32_c);
  tmp.fill(k);
  
  return tmp;
  }



//! scalar * Base
template<typename T1>
arma_inline
const OpCube<T1, op_scalar_times>
operator*
(const typename T1::elem_type k, const BaseCube<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return OpCube<T1, op_scalar_times>(X.get_ref(),k);  // NOTE: order is swapped
  }



//! scalar * OpCube<cube,op_ones_full>
template<typename eT>
arma_inline
Cube<eT>
operator*
(const eT k, const OpCube<Cube<eT>,op_ones_full>& X)
  {
  arma_extra_debug_sigprint();
    
  Cube<eT> tmp(X.aux_u32_a, X.aux_u32_b, X.aux_u32_c);
  tmp.fill(k);
  
  return tmp;
  }



//! @}
