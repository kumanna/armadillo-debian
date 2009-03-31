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


//! generate matrix/vector with all elements set to zero

class op_zeros
  {
  public:
  
  template<typename eT>
  inline static void apply(basic_mat<eT>&    out, const op_data<basic_mat<eT>,   op_zeros>& in);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>&    out, const op_data<basic_colvec<eT>,op_zeros>& in);
  
  template<typename eT>
  inline static void apply(basic_mat<eT>&    out, const op_data<basic_rowvec<eT>,op_zeros>& in);
  
  template<typename eT>
  inline static void apply(basic_colvec<eT>& out, const op_data<basic_colvec<eT>,op_zeros>& in);
  
  template<typename eT>
  inline static void apply(basic_rowvec<eT>& out, const op_data<basic_rowvec<eT>,op_zeros>& in);
  };


//! @}
