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


template<typename T1>
inline
void
op_zeros::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_zeros>& in)
  {
  arma_extra_debug_sigprint();
  
  out.zeros(in.aux_u32_a, in.aux_u32_b);
  }



template<typename eT>
inline
void
op_zeros::apply(Cube<eT>& out, const OpCube<Cube<eT>,op_zeros>& in)
  {
  arma_extra_debug_sigprint();
  
  out.zeros(in.aux_u32_a, in.aux_u32_b, in.aux_u32_c);
  }



//! @}
