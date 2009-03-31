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


//! \addtogroup op_ones
//! @{

template<typename eT>
inline
void
op_ones_full::apply(basic_mat<eT>& out, const op_data<basic_mat<eT>,op_ones_full>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(in.aux_u32_a, in.aux_u32_b);
  out.fill(eT(1));
  }



template<typename eT>
inline
void
op_ones_full::apply(basic_mat<eT>& out, const op_data<basic_colvec<eT>,op_ones_full>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(in.aux_u32_a, 1);
  out.fill(eT(1));
  }



template<typename eT>
inline
void
op_ones_full::apply(basic_mat<eT>& out, const op_data<basic_rowvec<eT>,op_ones_full>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(1, in.aux_u32_a);
  out.fill(eT(1));
  }



template<typename eT>
inline
void
op_ones_full::apply(basic_colvec<eT>& out, const op_data<basic_colvec<eT>,op_ones_full>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(in.aux_u32_a);
  out.fill(eT(1));
  }



template<typename eT>
inline
void
op_ones_full::apply(basic_rowvec<eT>& out, const op_data<basic_rowvec<eT>,op_ones_full>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(in.aux_u32_a);
  out.fill(eT(1));
  }



template<typename eT>
inline
void
op_ones_diag::apply(basic_mat<eT>& out, const op_data<basic_mat<eT>,op_ones_diag>& in)
  {
  arma_extra_debug_sigprint();
  
  out.zeros(in.aux_u32_a, in.aux_u32_b);
  
  for(u32 i=0; i<out.n_rows; ++i)
    {
    out.at(i,i) = eT(1);
    }
  
  }


//! @}
