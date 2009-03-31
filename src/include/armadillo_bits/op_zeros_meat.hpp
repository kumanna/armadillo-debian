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


//! \addtogroup op_zeros
//! @{


template<typename eT>
inline
void
op_zeros::apply(basic_mat<eT>& out, const op_data<basic_mat<eT>,op_zeros>& in)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_rows = in.aux_u32_a;
  const u32 n_cols = (in.aux_u32_b > 0) ? in.aux_u32_b : 1;
  
  out.zeros(n_rows, n_cols);
  }



template<typename eT>
inline
void
op_zeros::apply(basic_mat<eT>& out, const op_data<basic_colvec<eT>,op_zeros>& in)
  {
  arma_extra_debug_sigprint();
  
  out.zeros(in.aux_u32_a, 1);
  }



template<typename eT>
inline
void
op_zeros::apply(basic_mat<eT>& out, const op_data<basic_rowvec<eT>,op_zeros>& in)
  {
  arma_extra_debug_sigprint();
  
  out.zeros(1, in.aux_u32_a);
  }



template<typename eT>
inline
void
op_zeros::apply(basic_colvec<eT>& out, const op_data<basic_colvec<eT>,op_zeros>& in)
  {
  arma_extra_debug_sigprint();
  
  out.zeros(in.aux_u32_a);
  }



template<typename eT>
inline
void
op_zeros::apply(basic_rowvec<eT>& out, const op_data<basic_rowvec<eT>,op_zeros>& in)
  {
  arma_extra_debug_sigprint();
  
  out.zeros(in.aux_u32_a);
  }




//! @}
