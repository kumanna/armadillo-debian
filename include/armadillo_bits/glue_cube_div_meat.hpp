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


//! \addtogroup glue_cube_div
//! @{



template<typename T1, typename T2>
inline
void
glue_cube_div::apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_cube_div>& X)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const u32 N_cube = 1 + depth_lhs_cube< glue_cube_div, GlueCube<T1,T2,glue_cube_div> >::num;
  arma_extra_debug_print( arma_boost::format("N_cube = %d") % N_cube );

  if(N_cube == 2)
    {
    const unwrap_cube<T1> tmp1(X.A);
    const unwrap_cube<T2> tmp2(X.B);
    
    glue_cube_div::apply(out, tmp1.M, tmp2.M);
    }
  else
    {
    const Cube<eT>* ptrs[N_cube];
    bool             del[N_cube];

    cube_ptrs<glue_cube_div, GlueCube<T1,T2,glue_cube_div> >::get_ptrs(ptrs, del, X);

    for(u32 i=0; i<N_cube; ++i)  arma_extra_debug_print( arma_boost::format("ptrs[%d] = %x") % i % ptrs[i] );
    for(u32 i=0; i<N_cube; ++i)  arma_extra_debug_print( arma_boost::format(" del[%d] = %d") % i %  del[i] );

    const Cube<eT>& tmp_cube = *(ptrs[0]);

    for(u32 i=1; i<N_cube; ++i)
      {
      arma_debug_assert_same_size(tmp_cube, *(ptrs[i]), "element-wise cube division");
      }
  
    
    const u32 n_rows   = ptrs[0]->n_rows;
    const u32 n_cols   = ptrs[0]->n_cols;
    const u32 n_slices = ptrs[0]->n_slices;

    // no aliasing problem
    out.set_size(n_rows, n_cols, n_slices);
    
    const u32 n_elem = ptrs[0]->n_elem;
    
    for(u32 j=0; j<n_elem; ++j)
      {
      eT acc = ptrs[0]->mem[j];
      
      for(u32 i=1; i<N_cube; ++i)
        {
        acc /= ptrs[i]->mem[j];
        }
      
      out[j] = acc;
      }
    
    
    for(u32 i=0; i<N_cube; ++i)
      {
      if(del[i] == true)
        {
        arma_extra_debug_print( arma_boost::format("delete ptrs[%d]") % i );
        delete ptrs[i];
        }
      }

    }
  }



template<typename T1>
inline
void
glue_cube_div::apply_inplace(Cube<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> tmp(X);
  const Cube<eT>& B   = tmp.M;
  
  arma_debug_assert_same_size(out, B, "element-wise cube division");
  
  const u32 n_elem = out.n_elem;
  
        eT* out_mem = out.memptr();
  const eT* B_mem   = B.mem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] /= B_mem[i];
    }
  
  }



//! element-wise division with different element types
template<typename eT1, typename eT2>
inline
void
glue_cube_div::apply_mixed(Cube<typename promote_type<eT1,eT2>::result>& out, const Cube<eT1>& X, const Cube<eT2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  arma_debug_assert_same_size(X,Y, "element-wise cube division");
  
  //out.set_size(X.n_rows, X.n_cols, X.n_slices);
  out.copy_size(X);
  
        out_eT* out_mem = out.memptr();
  const eT1*    X_mem   = X.mem;
  const eT2*    Y_mem   = Y.mem;
  
  const u32 n_elem = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(X_mem[i]) / upgrade_val<eT1,eT2>::apply(Y_mem[i]);
    }
  }



template<typename eT>
inline
void
glue_cube_div::apply(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "element-wise cube division");
    
  // no aliasing problem
  //out.set_size(A.n_rows, A.n_cols, A.n_slices);
  out.copy_size(A);
  
  const u32 n_elem = A.n_elem;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A_mem[i] / B_mem[i];
    }
  
  }



template<typename eT>
inline
void
glue_cube_div::apply(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B, const Cube<eT>& C)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "element-wise cube division");
  arma_debug_assert_same_size(B, C, "element-wise cube division");
  
  // no aliasing problem
  //out.set_size(A.n_rows, A.n_cols, A.n_slices);
  out.copy_size(A);
    
  const u32 n_elem = A.n_elem;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  const eT* C_mem   = C.mem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A_mem[i] / B_mem[i] / C_mem[i];
    }
  
  }



#if defined(ARMA_GOOD_COMPILER)



template<typename eT>
inline
void
glue_cube_div::apply(Cube<eT>& out, const GlueCube<Cube<eT>,Cube<eT>,glue_cube_div>& X)
  {
  glue_cube_div::apply(out, X.A, X.B);
  }



template<typename eT>
inline
void
glue_cube_div::apply(Cube<eT>& out, const GlueCube< GlueCube<Cube<eT>,Cube<eT>,glue_cube_div>, Cube<eT>,glue_cube_div>& X)
  {
  glue_cube_div::apply(out, X.A.A, X.A.B, X.B);
  }



#endif



//! @}
