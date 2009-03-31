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


//! \addtogroup operartor_times_dot
//! @{



//! Row vector multiplied by a column vector (i.e. a dot product)
template<typename eT>
inline
eT
operator*
  (
  const basic_rowvec<eT>& A,
  const basic_colvec<eT>& B
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  
  }



//! Transpose of a column vector multiplied by a column vector (i.e. a dot product)
template<typename eT>
inline
eT
operator*
  (
  const op_data<basic_colvec<eT>, op_trans>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.m;
  const basic_colvec<eT>& B = Y;
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  
  }



//! Row vector multiplied by the transpose of a row vector (i.e. a dot product)
template<typename eT>
inline
eT
operator*
  (
  const basic_rowvec<eT>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_rowvec<eT>& A = X;
  const basic_rowvec<eT>& B = Y.m;
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  }
  

  
//! Transpose of a column vector multiplied by the transpose of a row vector (i.e. a dot product)
template<typename eT>
inline
eT
operator*
  (
  const op_data<basic_colvec<eT>, op_trans>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.m;
  const basic_rowvec<eT>& B = Y.m;
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  }



//
//
//



//! rowvec * mat * colvec 
template<typename eT>
inline
eT
operator*
  (
  const glue_data<basic_rowvec<eT>, basic_mat<eT>, glue_times_vec>& X,
  const basic_colvec<eT>& Y)
  {
  arma_extra_debug_sigprint();
  
  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * mat * colvec 
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>,op_trans>, basic_mat<eT>, glue_times>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! rowvec * mat * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data<basic_rowvec<eT>, basic_mat<eT>, glue_times_vec>& X,
  const op_data<basic_rowvec<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * mat * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>,op_trans>, basic_mat<eT>, glue_times>& X,
  const op_data<basic_rowvec<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//
//
// 


//! rowvec * diagmat(mat) * colvec 
template<typename eT>
inline
eT
operator*
  (
  const glue_data<basic_rowvec<eT>, op_data<basic_mat<eT>,op_diagmat>, glue_times_diag>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * diagmat(mat) * colvec 
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>,op_trans>, op_data<basic_mat<eT>,op_diagmat>, glue_times_diag>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//! rowvec * diagmat(mat) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data<basic_rowvec<eT>, op_data<basic_mat<eT>,op_diagmat>, glue_times_diag>& X,
  const op_data<basic_rowvec<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * diagmat(mat) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>,op_trans>, op_data<basic_mat<eT>,op_diagmat>, glue_times_diag>& X,
  const op_data<basic_rowvec<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//
//
//

//! rowvec * diagmat(colvec or rowvec) * colvec
template<typename eT>
inline
eT
operator*
  (
  const glue_data< basic_rowvec<eT>, op_data<basic_mat<eT>,op_diagmat_vec>, glue_times_vec>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//! trans(colvec) * diagmat(colvec or rowvec) * colvec
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>, op_trans>, op_data<basic_mat<eT>,op_diagmat_vec>, glue_times>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//! rowvec * diagmat(colvec or rowvec) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< basic_rowvec<eT>, op_data<basic_mat<eT>,op_diagmat_vec>, glue_times_vec>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();

  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//! trans(colvec) * diagmat(colvec or rowvec) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>, op_trans>, op_data<basic_mat<eT>,op_diagmat_vec>, glue_times>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//
// 
// 



//! rowvec * inv(T1) * colvec
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const glue_data< basic_rowvec<typename T1::elem_type>, op_data<T1, op_inv>, glue_times_vec>& X,
  const basic_colvec<typename T1::elem_type>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const basic_rowvec<eT>& A = X.A;
        basic_mat<eT>     B; op_inv::apply(B, X.B.m);
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(T1) * colvec
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const glue_data< op_data<basic_colvec<typename T1::elem_type>, op_trans>, op_data<T1, op_inv>, glue_times>& X,
  const basic_colvec<typename T1::elem_type>& Y
  )
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  const basic_colvec<eT>& A = X.A.m;
  basic_mat<eT>           B; op_inv::apply(B, X.B.m);
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! rowvec * inv(T1) * trans(rowvec)
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const glue_data< basic_rowvec<typename T1::elem_type>, op_data<T1, op_inv>, glue_times_vec>& X,
  const op_data<basic_rowvec<typename T1::elem_type>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;

  const basic_rowvec<eT>& A = X.A;
  basic_mat<eT> B;  op_inv::apply(B, X.B.m);
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(T1) * trans(rowvec)
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const glue_data< op_data<basic_colvec<typename T1::elem_type>, op_trans>, op_data<T1, op_inv>, glue_times>& X,
  const op_data<basic_rowvec<typename T1::elem_type>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const basic_colvec<eT>& A = X.A.m;
  basic_mat<eT>           B;  op_inv::apply(B, X.B.m);
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



// 
//
//



//! rowvec * inv(diagmat(mat)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const glue_data< basic_rowvec<eT>, op_data< op_data<basic_mat<eT>,op_diagmat>, op_inv>, glue_times_vec>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(mat)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>, op_trans>, op_data< op_data<basic_mat<eT>,op_diagmat>, op_inv>, glue_times>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//! rowvec * inv(diagmat(mat)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< basic_rowvec<eT>, op_data< op_data<basic_mat<eT>,op_diagmat>, op_inv>, glue_times_vec>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();

  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(mat)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>, op_trans>, op_data< op_data<basic_mat<eT>,op_diagmat>, op_inv>, glue_times>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//
//
//



//! rowvec * inv(diagmat(colvec or rowvec)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const glue_data< basic_rowvec<eT>, op_data< op_data<basic_mat<eT>,op_diagmat_vec>, op_inv>, glue_times_vec>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(colvec or rowvec)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>, op_trans>, op_data< op_data<basic_mat<eT>,op_diagmat_vec>, op_inv>, glue_times>& X,
  const basic_colvec<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_colvec<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//! rowvec * inv(diagmat(colvec or rowvec)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< basic_rowvec<eT>, op_data< op_data<basic_mat<eT>,op_diagmat_vec>, op_inv>, glue_times_vec>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();

  const basic_rowvec<eT>& A = X.A;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(colvec or rowvec)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data< op_data<basic_colvec<eT>, op_trans>, op_data< op_data<basic_mat<eT>,op_diagmat_vec>, op_inv>, glue_times>& X,
  const op_data<basic_rowvec<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A = X.A.m;
  const basic_mat<eT>&    B = X.B.m.m;
  const basic_rowvec<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//
//
//



//! trans(colvec-colvec) * inv(diagmat(colvec or rowvec)) * (colvec-colvec)
template<typename eT>
inline
eT
operator*
  (
  const glue_data
    <
    op_data< glue_data<basic_colvec<eT>, basic_colvec<eT>, glue_minus>, op_trans>,
    op_data< op_data<basic_mat<eT>,op_diagmat_vec>, op_inv>,
    glue_times
    >& X,
  const glue_data<basic_colvec<eT>, basic_colvec<eT>, glue_minus>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const basic_colvec<eT>& A1 = X.A.m.A;
  const basic_colvec<eT>& A2 = X.A.m.B;
  
  const basic_mat<eT>&    B  = X.B.m.m;
  
  const basic_colvec<eT>& C1 = Y.A;
  const basic_colvec<eT>& C2 = Y.B;
  
  arma_debug_check
    (
    (A1.n_rows != A2.n_rows)
    ||
    (A2.n_rows != B.n_rows)
    ||
    (B.n_elem != C1.n_rows)
    ||
    (C1.n_rows != C2.n_rows)
    ,
    "operator*(): incompatible matrix dimensions"
    );
  
  
  if( (&A1 == &C1) && (&A2 == &C2) )
    {
    eT val = eT(0);
    for(u32 i=0; i<A1.n_rows; ++i)
      {
      const eT tmp = (A1.mem[i] - A2.mem[i]);
      val += (tmp*tmp) / B.mem[i];
      }
    return val;
    }
  else
    {
    eT val = eT(0);
    for(u32 i=0; i<A1.n_rows; ++i)
      {
      val += ( (A1.mem[i] - A2.mem[i]) * (C1.mem[i] - C2.mem[i]) ) / B.mem[i];
      }
    return val;
    }

  }



//! @}
