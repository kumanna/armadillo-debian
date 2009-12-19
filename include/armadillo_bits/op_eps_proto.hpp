// Copyright (C) 2009 NICTA
// Copyright (C) 2009 Dimitrios Bouzas
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// - Dimitrios Bouzas (dimitris dot mpouzas at gmail dot com)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_eps
//! @{



class op_eps
  {
  public:
  
  template<typename eT> inline static   eT direct_eps(const eT               x);
  template<typename  T> inline static    T direct_eps(const std::complex<T>& x);
  
  template<typename eT> inline static void direct_eps(Mat<eT>& out, const Mat< eT              >& A);
  template<typename  T> inline static void direct_eps(Mat< T>& out, const Mat< std::complex<T> >& A);
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_eps>& in);
  };



//! @}
