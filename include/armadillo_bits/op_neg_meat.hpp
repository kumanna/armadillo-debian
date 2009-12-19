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



//! \addtogroup op_neg
//! @{



//! Negate all elements of a dense matrix
template<typename T1>
inline
void
op_neg::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_neg> &in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;  
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = -A_mem[i];
    }
  }




//! Negate all elements of a dense cube
template<typename T1>
inline
void
op_neg::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_neg> &in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;  
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = -A_mem[i];
    }
  }



//! @}
