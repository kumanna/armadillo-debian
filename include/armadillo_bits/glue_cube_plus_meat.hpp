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


//! \addtogroup glue_cube_plus
//! @{



template<typename T1, typename T2>
inline
void
glue_cube_plus::apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_cube_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 N_cube = 1 + depth_lhs_cube< glue_cube_plus, GlueCube<T1,T2,glue_cube_plus> >::num;
  arma_extra_debug_print( arma_boost::format("N_cube = %d") % N_cube );

  if(N_cube == 2)
    {
    const unwrap_cube<T1> tmp1(X.A);
    const unwrap_cube<T2> tmp2(X.B);
    
    glue_cube_plus::apply(out, tmp1.M, tmp2.M);
    }
  else
    {
    const Cube<eT>* ptrs[N_cube];
    bool             del[N_cube];

    cube_ptrs<glue_cube_plus, GlueCube<T1,T2,glue_cube_plus> >::get_ptrs(ptrs, del, X);

    for(u32 i=0; i<N_cube; ++i)  arma_extra_debug_print( arma_boost::format("ptrs[%d] = %x") % i % ptrs[i] );
    for(u32 i=0; i<N_cube; ++i)  arma_extra_debug_print( arma_boost::format(" del[%d] = %d") % i %  del[i] );

    const u32 n_rows   = ptrs[0]->n_rows;
    const u32 n_cols   = ptrs[0]->n_cols;
    const u32 n_slices = ptrs[0]->n_slices;
  
    const Cube<eT>& tmp_cube = *(ptrs[0]);
    
    for(u32 i=1; i<N_cube; ++i)
      {
      arma_debug_assert_same_size(tmp_cube, *(ptrs[i]), "cube addition");
      }
  
  
    // no aliasing problem
    out.set_size(n_rows, n_cols, n_slices);
    
    const u32 n_elem = ptrs[0]->n_elem;
    
    for(u32 j=0; j<n_elem; ++j)
      {
      eT acc = ptrs[0]->mem[j];
    
      for(u32 i=1; i < N_cube; ++i)
        {
        acc += ptrs[i]->mem[j];
        }
    
      out[j] = acc;
      }
    
    // TODO:CUBE the cube_ptrs class should have a function to do the deletion
    // TODO:MAT in a similar vein, the mat_ptrs should also have a function to do the deletion
    
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
glue_cube_plus::apply_inplace(Cube<typename T1::elem_type>& out, const T1& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> tmp(X);
  const Cube<eT>& B   = tmp.M;

  arma_debug_assert_same_size(out, B, "cube addition");
  
  
        eT* out_mem = out.memptr();
  const eT* B_mem   = B.mem;
  
  const u32 n_elem  = B.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] += B_mem[i];
    }
  
  }



//! cube addition with different element types
template<typename eT1, typename eT2>
inline
void
glue_cube_plus::apply_mixed(Cube<typename promote_type<eT1,eT2>::result>& out, const Cube<eT1>& X, const Cube<eT2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  arma_debug_assert_same_size(X,Y, "cube addition");  
  
  //out.set_size(X.n_rows, X.n_cols, X.n_slices);
  out.copy_size(X);
  
        out_eT* out_mem = out.memptr();
  const eT1*    X_mem   = X.mem;
  const eT2*    Y_mem   = Y.mem;
  
  const u32 n_elem = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(X_mem[i]) + upgrade_val<eT1,eT2>::apply(Y_mem[i]);
    }
  }



template<typename eT>
inline
void
glue_cube_plus::apply(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "cube addition");
  
  // no aliasing problem
  //out.set_size(A.n_rows, A.n_cols, A.n_slices);
  out.copy_size(A);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
    
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A_mem[i] + B_mem[i];
    }
    
  }



template<typename eT>
inline
void
glue_cube_plus::apply(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B, const Cube<eT>& C)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "cube addition");
  arma_debug_assert_same_size(B, C, "cube addition");
  
  // no aliasing problem
  //out.set_size(A);
  out.copy_size(A);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  const eT* C_mem   = C.mem;
  
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A_mem[i] + B_mem[i] + C_mem[i];
    }
    
  }



#if defined(ARMA_GOOD_COMPILER)



template<typename eT>
inline
void
glue_cube_plus::apply(Cube<eT>& out, const GlueCube<Cube<eT>,Cube<eT>,glue_cube_plus>& X)
  {
  glue_cube_plus::apply(out, X.A, X.B);
  }



template<typename eT>
inline
void
glue_cube_plus::apply(Cube<eT>& out, const GlueCube< GlueCube<Cube<eT>,Cube<eT>,glue_cube_plus>, Cube<eT>, glue_cube_plus>& X)
  {
  glue_cube_plus::apply(out, X.A.A, X.A.B, X.B);
  }



template<typename eT>
inline
void
glue_cube_plus::apply(Cube<eT>& out, const GlueCube<Cube<eT>, subview_cube<eT>, glue_cube_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT>& orig_A = X.A;
  const Cube<eT>& orig_B = X.B.m;
  
  if( &out != &orig_B )
    {
    arma_debug_assert_same_size(X.A, X.B, "cube addition");
    
    //out.set_size(orig_A.n_rows, orig_A.n_cols, orig_A.n_slices);
    out.copy_size(orig_A);


    for(u32 slice = 0; slice < orig_A.n_slices; ++slice)
      {
      const u32 B_slice_mod = X.B.aux_slice1 + slice;
      
      for(u32 col = 0; col < orig_A.n_cols; ++col)
        {
        const u32 B_col_mod = X.B.aux_col1 + col;
        
        for(u32 row = 0; row < orig_A.n_rows; ++row)
          {
          const u32 B_row_mod = X.B.aux_row1 + row;
          
          out.at(row,col,slice) =  orig_A.at(row, col,slice) + orig_B.at(B_row_mod, B_col_mod, B_slice_mod);
          }
        }
      }
    
    }
  else
    {
    const Cube<eT> processed_B(X.B);  // create a cube out of subview_cube
    glue_cube_plus::apply(out, orig_A, processed_B);
    }
     
  }



template<typename eT>
inline
void
glue_cube_plus::apply(Cube<eT>& out, const GlueCube< subview_cube<eT>, Cube<eT>, glue_cube_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT>& orig_A = X.A.m;
  
  const unwrap_cube_check< Cube<eT> > tmp(X.B, out);
  const Cube<eT>& orig_B = tmp.M;
  
  if( &out != &orig_A )
    {
    const u32 sub_A_n_rows   = X.A.n_rows;
    const u32 sub_A_n_cols   = X.A.n_cols;
    const u32 sub_A_n_slices = X.A.n_slices;
    
    arma_debug_assert_same_size(X.A, X.B, "cube addition");
      
    out.set_size(sub_A_n_rows, sub_A_n_cols, sub_A_n_slices);

    for(u32 slice = 0; slice < sub_A_n_slices; ++slice)
      {
      const u32 A_slice_mod = X.A.aux_slice1 + slice;
        
      for(u32 col = 0; col < sub_A_n_cols; ++col)
        {
        const u32 A_col_mod = X.A.aux_col1 + col;
        
        for(u32 row = 0; row < sub_A_n_rows; ++row)
          {
          const u32 A_row_mod = X.A.aux_row1 + row;
          
          out.at(row,col,slice) =  orig_A.at(A_row_mod, A_col_mod, A_slice_mod) + orig_B.at(row, col, slice);
          }
        }
      }
    }
  else
    {
    const Cube<eT> processed_A(X.A);
    glue_cube_plus::apply(out, processed_A, orig_B);
    }
  
  
  }



template<typename eT>
inline
void
glue_cube_plus::apply(Cube<eT>& out, const GlueCube< subview_cube<eT>, subview_cube<eT>, glue_cube_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  const Cube<eT>& orig_A = X.A.m;
  const Cube<eT>& orig_B = X.B.m;
  
  if( (&out != &orig_A) && (&out != &orig_B) )
    {
    const u32 sub_A_n_rows   = X.A.n_rows;
    const u32 sub_A_n_cols   = X.A.n_cols;
    const u32 sub_A_n_slices = X.A.n_slices;
    
    //const u32 sub_B_n_rows = X.B.n_rows;
    //const u32 sub_B_n_cols = X.B.n_cols;
    
    arma_debug_assert_same_size(X.A, X.B, "cube addition");
      
    out.set_size(sub_A_n_rows, sub_A_n_cols, sub_A_n_slices);

    for(u32 slice = 0; slice < sub_A_n_slices; ++slice)
      {
      const u32 A_slice_mod = X.A.aux_slice1 + slice;
      const u32 B_slice_mod = X.B.aux_slice1 + slice;
        
      for(u32 col = 0; col < sub_A_n_cols; ++col)
        {
        const u32 A_col_mod = X.A.aux_col1 + col;
        const u32 B_col_mod = X.B.aux_col1 + col;
        
        for(u32 row = 0; row < sub_A_n_rows; ++row)
          {
          const u32 A_row_mod = X.A.aux_row1 + row;
          const u32 B_row_mod = X.B.aux_row1 + row;
          
          out.at(row,col,slice) =  orig_A.at(A_row_mod, A_col_mod, A_slice_mod) + orig_B.at(B_row_mod, B_col_mod, B_slice_mod);
          }
        }
      }
    }
  else
    {
    const Cube<eT> processed_A(X.A);
    const Cube<eT> processed_B(X.B);
    
    glue_cube_plus::apply(out, processed_A, processed_B);
    }
  }



template<typename T1>
inline
void
glue_cube_plus::apply_inplace(Cube<typename T1::elem_type>& out, const OpCube<T1, op_square>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube<T1> tmp(X.m);
  const Cube<eT>& B   = tmp.M;
  
  arma_debug_assert_same_size(out, B, "cube addition");
    
        eT* out_mem = out.memptr();
  const eT* B_mem   = B.mem;
  
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const eT tmp_val = B_mem[i];
    out_mem[i] += tmp_val*tmp_val;
    }
  }



template<typename T1, typename T2>
inline
void
glue_cube_plus::apply_inplace(Cube<typename T1::elem_type>& out, const GlueCube<T1, T2, glue_cube_plus>& X)
  {
  arma_extra_debug_sigprint();
    
  out = X + out;
  }



#endif



//! @}
