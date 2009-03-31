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


//! \addtogroup glue_div
//! @{


//! Class which implements the immediate element-wise division of two or more matrices
class glue_div
  {
  public:
  
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const basic_mat<eT>& A, const basic_mat<eT>& B, const basic_mat<eT>& C);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data<basic_mat<eT>,basic_mat<eT>,glue_div>& X);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const glue_data< glue_data<basic_mat<eT>,basic_mat<eT>, glue_div>, basic_mat<eT>, glue_div>& X);
  
  template<typename T1, typename T2>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const glue_data<T1,T2,glue_div>& X);
  
  
  template<typename eT>
  inline static void apply_inplace(basic_mat<eT>& out, const basic_mat<eT>& B);
  
  template<typename T1, typename op_type>
  inline static void apply_inplace(basic_mat<typename T1::elem_type>& out, const op_data<T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply_inplace(basic_mat<typename T1::elem_type>& out, const glue_data<T1, T2, glue_type>& X);


  // element-wise division with different element types
  
  template<typename eT1, typename eT2>
  inline static void apply_mixed(basic_mat<typename promote_type<eT1,eT2>::result>& out, const basic_mat<eT1>& X, const basic_mat<eT2>& Y);
  
  };



//! @}

