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


//! \addtogroup glue_times
//! @{


//! Class which implements the immediate multiplication of two or more matrices
class glue_times
  {
  public:
  
  
  template<typename eT>
  inline static u32  mul_storage_cost(const basic_mat<eT>& X, const basic_mat<eT>& Y);
  
  template<typename eT>
  inline static void apply_noalias(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B, const basic_mat<eT>& C);
  
  
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data<basic_mat<eT>,basic_mat<eT>,glue_times>& X);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data< glue_data<basic_mat<eT>,basic_mat<eT>, glue_times>, basic_mat<eT>, glue_times>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1,T2,glue_times>& X);
  
  
  
  template<typename eT>
  inline static void apply_inplace(basic_mat<eT>& out, const basic_mat<eT>& B);
  
  template<typename T1, typename op_type>
  inline static void apply_inplace(basic_mat<typename T1::elem_type>& out, const op_data<T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply_inplace(basic_mat<typename T1::elem_type>& out, const glue_data<T1, T2, glue_type>& X);
  
  
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1,op_data<T2,op_trans>,glue_times >& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<op_data<T1,op_trans>,T2,glue_times>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<op_data<T1,op_trans>,op_data<T2,op_trans>,glue_times>& X);
  
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data< op_data<T1, op_neg>, T2, glue_times>& X);
  
  
  template<typename eT>
  inline static eT direct_rowvec_mat_colvec(const eT* A_mem, const basic_mat<eT>& B, const eT* C_mem);

  template<typename eT>
  inline static eT direct_rowvec_diagmat_colvec(const eT* A_mem, const basic_mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagmat_colvec(const eT* A_mem, const basic_mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagvec_colvec(const eT* A_mem, const basic_mat<eT>& B, const eT* C_mem);
  
  
  // matrix multiplication with different element types
  
  template<typename eT1, typename eT2>
  inline static void apply_mixed(basic_mat<typename promote_type<eT1,eT2>::result>& out, const basic_mat<eT1>& X, const basic_mat<eT2>& Y);
  
  };



class glue_times_diag
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const T1& A, const op_data<T2,op_diagmat>& B);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const op_data<T1,op_diagmat>& A, const T2& B);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const op_data<T1,op_diagmat>& A, const op_data<T2,op_diagmat>& B);
  
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1, op_data<T2,op_diagmat>, glue_times_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<op_data<T1,op_diagmat>, T2, glue_times_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_times_diag>& X);
  
  };



class glue_times_vec
  {
  public:
  
  template<typename eT>
  inline static void mul_col_row(basic_mat<eT>& out, const eT* A_mem, const eT* B_mem);

  template<typename T1>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1, basic_colvec<typename T1::elem_type>, glue_times_vec>& X);
  
  template<typename T1>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1, basic_rowvec<typename T1::elem_type>, glue_times_vec>& X);
  
  template<typename T1>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<basic_colvec<typename T1::elem_type>, T1, glue_times_vec>& X);
  
  template<typename T1>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<basic_rowvec<typename T1::elem_type>, T1, glue_times_vec>& X);
  
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data<basic_colvec<eT>, basic_rowvec<eT>, glue_times_vec>& X);

  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data< op_data<basic_rowvec<eT>, op_trans>, basic_rowvec<eT>, glue_times_vec>& X);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data< basic_colvec<eT>, op_data<basic_colvec<eT>, op_trans>, glue_times_vec>& X);
  
  
  
  template<typename T1>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<op_data<T1, op_trans>, basic_colvec<typename T1::elem_type>, glue_times_vec>& X);
  
  };


//! @}

