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



//! \addtogroup op_scalar_misc
//! @{


//! 'add scalar to a matrix/cube' operation
class op_scalar_plus
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_scalar_plus>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_scalar_plus>& in);
  };



//! 'subtract matrix/cube from a scalar' operation
class op_scalar_minus_pre
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_scalar_minus_pre>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_scalar_minus_pre>& in);
  };



//! 'subtract scalar from a matrix/cube' operation
class op_scalar_minus_post
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_scalar_minus_post>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_scalar_minus_post>& in);
  };


//! 'multiply matrix/cube by a scalar' operation
class op_scalar_times
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_scalar_times>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_scalar_times>& in);
  
  #if defined(ARMA_GOOD_COMPILER)
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_plus>,  op_scalar_times>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_plus>,  op_scalar_times>& in);
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_minus>, op_scalar_times>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_minus>, op_scalar_times>& in);
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_schur>, op_scalar_times>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_schur>, op_scalar_times>& in);
  
  #endif
  };



//! 'divide scalar by a matrix/cube' operation
class op_scalar_div_pre
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_scalar_div_pre>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_scalar_div_pre>& in);
  
  #if defined(ARMA_GOOD_COMPILER)

  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_plus>,  op_scalar_div_pre>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_plus>,  op_scalar_div_pre>& in);
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_minus>, op_scalar_div_pre>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_minus>, op_scalar_div_pre>& in);
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_schur>, op_scalar_div_pre>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_schur>, op_scalar_div_pre>& in);

  #endif
  };



//! 'divide matrix/cube by a scalar' operation
class op_scalar_div_post
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_scalar_div_post>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_scalar_div_post>& in);
  
  
  #if defined(ARMA_GOOD_COMPILER)
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_plus>,  op_scalar_div_post>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_plus>,  op_scalar_div_post>& in);
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_minus>, op_scalar_div_post>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_minus>, op_scalar_div_post>& in);
  
  template<typename T1, typename T2> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<    Glue<T1,T2,     glue_schur>, op_scalar_div_post>& in);
  template<typename T1, typename T2> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<GlueCube<T1,T2,glue_cube_schur>, op_scalar_div_post>& in);
  
  #endif
  };



//! @}
