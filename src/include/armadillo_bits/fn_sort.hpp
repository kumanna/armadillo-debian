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


//! \addtogroup fn_sort
//! @{


template<typename T1>
inline
const op_data<T1, op_sort>
sort(const arma_base<typename T1::elem_type,T1>& X, const u32 sort_type = 0, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();
  
  return op_data<T1, op_sort>(X.get_ref(), sort_type, dim);
  }



template<typename eT>
inline
const op_data<basic_colvec<eT>, op_sort>
sort(const basic_colvec<eT>& X, const u32 sort_type = 0)
  {
  arma_extra_debug_sigprint();
  
  const u32 dim = 0;
  
  return op_data<basic_colvec<eT>, op_sort>(X, sort_type, dim);
  }



template<typename eT>
inline
const op_data<basic_rowvec<eT>, op_sort>
sort(const basic_rowvec<eT>& X, const u32 sort_type = 0)
  {
  arma_extra_debug_sigprint();
  
  const u32 dim = 1;
  
  return op_data<basic_rowvec<eT>, op_sort>(X, sort_type, dim);
  }



//! @}
