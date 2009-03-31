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


//! \addtogroup fn_eig
//! @{


template<typename eT, typename T1>
inline
void
eig(basic_colvec<eT>& eigval, const arma_base<eT,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(A.get_ref());

  auxlib::eig(eigval, tmp.M);
  }



template<typename eT, typename T1>
inline
void
eig(basic_colvec<eT>& eigval, basic_mat<eT>& eigvec, const arma_base<eT,T1>& A)
  {
  arma_extra_debug_sigprint();

  const unwrap<T1> tmp(A.get_ref());
  
  auxlib::eig(eigval, eigvec, tmp.M);
  }


//! @}
