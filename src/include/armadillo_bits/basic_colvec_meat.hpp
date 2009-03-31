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


//! \addtogroup basic_colvec
//! @{


//! construct an empty column vector
template<typename eT>
inline
basic_colvec<eT>::basic_colvec()
  : basic_mat<eT>()
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector with the specified number of n_elem
template<typename eT>
inline
basic_colvec<eT>::basic_colvec(const u32 in_n_elem)
  : basic_mat<eT>(in_n_elem,1)
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector from specified text
template<typename eT>
inline
basic_colvec<eT>::basic_colvec(const char* text)
  : basic_mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  std::swap( access::rw(basic_mat<eT>::n_rows), access::rw(basic_mat<eT>::n_cols) );
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }



//! construct a column vector from specified text
template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(text);
  std::swap( access::rw(basic_mat<eT>::n_rows), access::rw(basic_mat<eT>::n_cols) );
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from a given column vector
template<typename eT>
inline
basic_colvec<eT>::basic_colvec(const basic_colvec<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector from a given column vector
template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator=(const basic_colvec<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  return *this;
  }



//! construct a column vector from a given matrix; the matrix must have exactly one column
template<typename eT>
inline
basic_colvec<eT>::basic_colvec(const basic_mat<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }



//! construct a column vector from a given matrix; the matrix must have exactly one column
template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator=(const basic_mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator*=(const basic_mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from a given auxillary array of eTs
template<typename eT>
inline
basic_colvec<eT>::basic_colvec(const eT* aux_mem, const u32 aux_length)
  {
  arma_extra_debug_sigprint();
  
  set_size(aux_length, 1);

  arma_check( (basic_mat<eT>::n_elem != aux_length), "basic_colvec::basic_colvec(): don't know how to handle the given array" );

  syslib::copy_elem( basic_mat<eT>::memptr(), aux_mem, basic_mat<eT>::n_elem );
  }



template<typename eT>
template<typename T1, typename T2>
inline
basic_colvec<eT>::basic_colvec
  (
  const arma_base<typename basic_colvec<eT>::pod_type, T1>& A,
  const arma_base<typename basic_colvec<eT>::pod_type, T2>& B
  )
  : basic_mat<eT>(A,B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }



//! construct a column vector from given a submatrix; the submatrix must have exactly one column
template<typename eT>
inline
basic_colvec<eT>::basic_colvec(const subview<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }



//! construct a column vector from given a submatrix; the submatrix must have exactly one column
template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from given a diagview
template<typename eT>
inline
basic_colvec<eT>::basic_colvec(const diagview<eT>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }



//! construct a column vector from given a diagview
template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator*=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  return *this;
  }



//! construct a column vector from op_data, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename op_type>
inline
basic_colvec<eT>::basic_colvec(const op_data<T1, op_type>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }



//! construct a column vector from op_data, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename op_type>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator=(const op_data<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec::operator=(): given matrix can't be interpreted as a column vector" );
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator*=(const op_data<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec::operator=(): incompatible dimensions" );
  return *this;
  }



//! construct a column vector from glue_data, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
basic_colvec<eT>::basic_colvec(const glue_data<T1, T2, glue_type>& X)
  : basic_mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }



//! construct a column vector from glue_data, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator=(const glue_data<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const basic_colvec<eT>&
basic_colvec<eT>::operator*=(const glue_data<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::operator*=(X);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  
  return *this;
  }



//! change the number of n_rows
template<typename eT>
inline
void
basic_colvec<eT>::set_size(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::set_size(in_n_elem,1);
  }



//! change the number of n_rows  (this function re-implements mat::set_size() in order to check the number of columns)
template<typename eT>
inline
void
basic_colvec<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::set_size( in_n_rows, (std::min)( u32(1), in_n_cols ) );
  arma_debug_check( (in_n_cols > 1), "basic_colvec::set_size(): incompatible dimensions" );
  }



template<typename eT>
inline
void
basic_colvec<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::zeros();
  }



template<typename eT>
inline
void
basic_colvec<eT>::zeros(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::zeros(in_n_elem,1);
  }



template<typename eT>
inline
void
basic_colvec<eT>::zeros(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::zeros( in_n_rows, (std::min)( u32(1), in_n_cols ) );
  arma_debug_check( (in_n_cols > 1), "basic_colvec::zeros(): incompatible dimensions" );
  }



template<typename eT>
inline
void
basic_colvec<eT>::load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<eT>::load(name,type);
  arma_debug_check( (basic_mat<eT>::n_cols > 1), "basic_colvec(): incompatible dimensions" );
  }


//! @}
