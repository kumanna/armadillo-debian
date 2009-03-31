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


//! \addtogroup glue_minus
//! @{



//! Class for the minus operation, where the result is always a dense matrix
class glue_minus
  {
  public:

  
  // mat

  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B, const basic_mat<eT>& C);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data< basic_mat<eT>, basic_mat<eT>, glue_minus>& X);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data< glue_data< basic_mat<eT>,basic_mat<eT>,glue_minus>, basic_mat<eT>, glue_minus>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1,T2,glue_minus>& X);

  
  // mat, inplace
  
  template<typename eT>
  inline static void apply_inplace(basic_mat<eT>& out, const basic_mat<eT>& B);
  
  template<typename T1, typename op_type>
  inline static void apply_inplace(basic_mat<typename T1::elem_type>& out, const op_data<T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply_inplace(basic_mat<typename T1::elem_type>& out, const glue_data<T1, T2, glue_type>& X);
  
  
  // matrix subtraction with different element types
  
  template<typename eT1, typename eT2>
  inline static void apply_mixed(basic_mat<typename promote_type<eT1,eT2>::result>& out, const basic_mat<eT1>& X, const basic_mat<eT2>& Y);
  
  };



//! Class for the minus operation, where one of the operands is a diagonal matrix
class glue_minus_diag
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const T1& A, const op_data<T2,op_diagmat>& B);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const op_data<T1,op_diagmat>& A, const T2& B);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const op_data<T1,op_diagmat>& A, const op_data<T2,op_diagmat>& B);
  
  //
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1, op_data<T2,op_diagmat>, glue_minus_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<op_data<T1,op_diagmat>, T2, glue_minus_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_minus_diag>& X);
  
  
  };

//! @}
