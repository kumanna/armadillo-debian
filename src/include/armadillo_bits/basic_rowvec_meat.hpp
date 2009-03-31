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


//! \addtogroup basic_rowvec
//! @{

template<typename eT>
inline
basic_rowvec<eT>::basic_rowvec()
  : basic_mat<eT>()
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
basic_rowvec<eT>::basic_rowvec(const u32 in_n_elem)
  : basic_mat<eT>(1,in_n_elem)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
basic_rowvec<eT>::basic_rowvec(const char* text)
  : basic_mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  }
  


template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(text);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
basic_rowvec<eT>::basic_rowvec(const basic_rowvec<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator=(const basic_rowvec<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  return *this;
  }



template<typename eT>
inline basic_rowvec<eT>::basic_rowvec(const basic_mat<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  }



template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator=(const basic_mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator*=(const basic_mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



//! construct a row vector from a given auxillary array
template<typename eT>
inline
basic_rowvec<eT>::basic_rowvec(const eT* aux_mem, const u32 aux_length)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::set_size(1, aux_length);
  arma_check( (basic_mat<eT>::n_elem != aux_length), "basic_rowvec(): don't know how to handle the given array" );

  syslib::copy_elem( basic_mat<eT>::memptr(), aux_mem, basic_mat<eT>::n_elem );
  }



template<typename eT>
template<typename T1, typename T2>
inline
basic_rowvec<eT>::basic_rowvec
  (
  const arma_base<typename basic_rowvec<eT>::pod_type, T1>& A,
  const arma_base<typename basic_rowvec<eT>::pod_type, T2>& B
  )
  : basic_mat<eT>(A,B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_colvec(): incompatible dimensions" );
  }



template<typename eT>
inline
basic_rowvec<eT>::basic_rowvec(const subview<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  }



template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



//! construct a row vector from given a diagview
template<typename eT>
inline basic_rowvec<eT>::basic_rowvec(const diagview<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  std::swap( access::rw(basic_mat<eT>::n_rows), access::rw(basic_mat<eT>::n_cols) );
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  }



//! construct a row vector from given a diagview
template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  
  //std::swap( access::rw(basic_mat<eT>::n_rows), access::rw(basic_mat<eT>::n_cols) );
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator*=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
basic_rowvec<eT>::basic_rowvec(const op_data<T1, op_type>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator=(const op_data<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator*=(const op_data<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
basic_rowvec<eT>::basic_rowvec(const glue_data<T1, T2, glue_type>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator=(const glue_data<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const basic_rowvec<eT>&
basic_rowvec<eT>::operator*=(const glue_data<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
void
basic_rowvec<eT>::set_size(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::set_size(1,in_n_elem);
  }



template<typename eT>
inline
void
basic_rowvec<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::set_size( (std::min)( u32(1), in_n_rows), in_n_cols );
  
  arma_debug_check( (in_n_rows > 1), "basic_rowvec::set_size(): incompatible dimensions" );
  }



template<typename eT>
inline
void
basic_rowvec<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::zeros();
  }



template<typename eT>
inline
void
basic_rowvec<eT>::zeros(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::zeros(1,in_n_elem);
  }



template<typename eT>
inline
void
basic_rowvec<eT>::zeros(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::zeros( (std::min)( u32(1), in_n_rows), in_n_cols );
  arma_debug_check( (in_n_rows > 1), "basic_rowvec<eT>::zeros(): incompatible dimensions" );
  }



template<typename eT>
inline
void
basic_rowvec<eT>::load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::load(name,type);
  arma_debug_check( (basic_mat<eT>::n_rows > 1), "basic_rowvec(): incompatible dimensions" );
  }



//! @}
