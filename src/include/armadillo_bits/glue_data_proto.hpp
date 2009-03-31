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


//! \addtogroup glue_data
//! @{


//! Class for storing data required for delayed binary operations,
//! such as the operands (e.g. two matrices) and the binary operator (e.g. addition).
//! The operands are stored as references (which can be optimised away),
//! while the operator is "stored" through the template definition (glue_type).
//! The operands can be 'mat', 'rowvec', 'colvec', 'op_data', and 'glue_data'.
//! Note that as 'glue_data' can be one of the operands, more than two matrices can be stored.
//!
//! For example, we could have: glue_data< glue_data<mat, mat, glue_times>, mat, glue_plus>
//! 
//! Another example is: glue_data< op_data<mat, op_trans>, op_data<mat, op_inv>, glue_times >
//!
//! More complicated example: glue_data< op_data< glue_data<mat, mat, glue_plus>, op_trans>, op_data<mat, op_inv>, glue_times >
//! 

template<typename T1, typename T2, typename glue_type>
class glue_data : public arma_base<typename T1::elem_type, glue_data<T1, T2, glue_type> >
  {
  public:
  
  typedef typename T1::elem_type elem_type;
  typedef typename get_pod_type<elem_type>::pod_type pod_type;

  inline  glue_data(const T1& in_A, const T2& in_B);
  inline ~glue_data();
  
  const T1& A;  //!< first operand
  const T2& B;  //!< second operand
  
  };


//! @}
