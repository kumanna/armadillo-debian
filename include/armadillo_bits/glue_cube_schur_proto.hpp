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


//! \addtogroup glue_cube_schur
//! @{


//! Class which implements the immediate Schur product (element-wise multiplication) of two or more cubes
class glue_cube_schur
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Cube<typename T1::elem_type>& out, const GlueCube<T1,T2,glue_cube_schur>& X);
  
  
  template<typename T1>
  inline static void apply_inplace(Cube<typename T1::elem_type>& out, const T1& X);
  
  
  template<typename eT1, typename eT2>
  inline static void apply_mixed(Cube<typename promote_type<eT1,eT2>::result>& out, const Cube<eT1>& X, const Cube<eT2>& Y);
  
  
  template<typename eT>
  inline static void apply(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B);
  
  template<typename eT>
  inline static void apply(Cube<eT>& out, const Cube<eT>& A, const Cube<eT>& B, const Cube<eT>& C);
  
  
  #if defined(ARMA_GOOD_COMPILER)
  
  
  template<typename eT>
  inline static void apply(Cube<eT>& out, const GlueCube<Cube<eT>,Cube<eT>,glue_cube_schur>& X);
  
  template<typename eT>
  inline static void apply(Cube<eT>& out, const GlueCube< GlueCube<Cube<eT>,Cube<eT>,glue_cube_schur>, Cube<eT>, glue_cube_schur>& X);
  
  template<typename T1, typename T2>
  inline static void apply_inplace(Cube<typename T1::elem_type>& out, const GlueCube<T1, T2, glue_cube_schur>& X);
  
  
  #endif
  
  };



//! @}

