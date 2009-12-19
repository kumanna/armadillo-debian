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


//! \addtogroup op_ones
//! @{


//! Class for creation of a dense matrix/vector/cube with all elements set to one
class op_ones_full
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_ones_full>& in);
  
  template<typename eT>
  inline static void apply(Cube<eT>& out, const OpCube<Cube<eT>,op_ones_full>& in);
  };



class op_ones_diag
  {
  public:
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Op<Mat<eT>,op_ones_diag>& in);
  };


//! @}
