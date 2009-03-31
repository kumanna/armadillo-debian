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


//! \addtogroup op_inv
//! @{



//! 'invert matrix' operation

class op_inv
  {
  public:
  
  // mat

  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const basic_mat<eT>& A);
  
  template<typename T1>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const op_data<T1,op_inv>& in);

  template<typename T1>
  inline static void apply(basic_mat<typename T1::elem_type>& out, const op_data< op_data<T1,op_diagmat>, op_inv>& in);

  //

  template<typename eT>
  inline static void apply_diagvec(basic_mat<eT>& out, const basic_mat<eT>& X);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>& out, const op_data< op_data<basic_mat<eT>,op_diagmat_vec>, op_inv>& in);
  
  };

//! @}
