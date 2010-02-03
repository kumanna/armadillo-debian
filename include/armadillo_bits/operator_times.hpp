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




//! \addtogroup operator_times
//! @{

// "rowvec * colvec" and other combinations which result
// in a scalar are declared in "operator_times_dot.hpp"

//
// new operators

//! Base * scalar
template<typename T1>
arma_inline
const Op<T1, op_scalar_times>
operator*
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_times>(X.get_ref(),k);
  }



//! op * scalar, level 2
template<typename T1>
arma_inline
const Op<T1,op_scalar_times>
operator*
(const Op<T1,op_scalar_times>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_times>(X.m, X.aux * k);
  }



//! Op<mat,op_ones_full> * scalar
template<typename eT>
arma_inline
Mat<eT>
operator*
(const Op<Mat<eT>,op_ones_full>& X, const eT k)
  {
  arma_extra_debug_sigprint();
    
  Mat<eT> tmp(X.aux_u32_a, X.aux_u32_b);
  tmp.fill(k);
  
  return tmp;
  }



//! Op<mat,op_ones_diag> * scalar
template<typename eT>
arma_inline
Mat<eT>
operator*
(const Op<Mat<eT>,op_ones_diag>& X, const eT k)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.aux_u32_a, X.aux_u32_b);
  
  for(u32 i=0; i<out.n_rows; ++i)
    {
    out.at(i,i) = k;
    }
  
  return out;
  }



//! scalar * Base
template<typename T1>
arma_inline
const Op<T1, op_scalar_times>
operator*
(const typename T1::elem_type k, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_times>(X.get_ref(),k);  // NOTE: order is swapped
  }



//! scalar * Op<mat,op_ones_full>
template<typename eT>
arma_inline
Mat<eT>
operator*
(const eT k, const Op<Mat<eT>,op_ones_full>& X)
  {
  arma_extra_debug_sigprint();
    
  Mat<eT> tmp(X.aux_u32_a, X.aux_u32_b);
  tmp.fill(k);
  
  return tmp;
  }



//! scalar * Op<mat,op_ones_diag>
template<typename eT>
arma_inline
Mat<eT>
operator*
(const eT k, const Op<Mat<eT>,op_ones_diag>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.aux_u32_a, X.aux_u32_b);
  
  for(u32 i=0; i<out.n_rows; ++i)
    out.at(i,i) = k;
  
  return out;
  }



//! Base * diagmat
template<typename T1, typename T2>
arma_inline
const Glue<T1, Op<T2,op_diagmat>, glue_times_diag>
operator*
(const Base<typename T2::elem_type, T1>& X, const Op<T2,op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Op<T2,op_diagmat>, glue_times_diag>(X.get_ref(), Y);
  }



//! diagmat * Base
template<typename T1, typename T2>
arma_inline
const Glue<Op<T1,op_diagmat>, T2, glue_times_diag>
operator*
(const Op<T1,op_diagmat>& X, const Base<typename T1::elem_type, T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1,op_diagmat>, T2, glue_times_diag>(X, Y.get_ref());
  }



//! diagmat * diagmat
template<typename T1, typename T2>
inline
Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
operator*
(const Op<T1, op_diagmat>& X, const Op<T2, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT1;
  typedef typename T2::elem_type eT2;
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.m);
  const unwrap<T2> tmp2(Y.m);
  
  const Mat<eT1> A(tmp1.M);
  const Mat<eT2> B(tmp2.M);
  
  arma_debug_assert_same_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  Mat<out_eT> out(A.n_rows, A.n_rows);
  out.zeros();
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    out.at(i,i) = upgrade_val<eT1,eT2>::apply( A.at(i,i) ) * upgrade_val<eT1,eT2>::apply( B.at(i,i) );
    }
  
  return out;
  }



//! colvec * rowvec
template<typename eT>
arma_inline
const Glue<Col<eT>, Row<eT>, glue_times_vec>
operator*
(const Col<eT>& X, const Row<eT>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Col<eT>, Row<eT>, glue_times_vec>(X, Y);
  }



//! Base * colvec
template<typename T1>
arma_inline
const Glue<T1, Col<typename T1::elem_type>, glue_times_vec>
operator*
(const Base<typename T1::elem_type,T1>& X, const Col<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Col<typename T1::elem_type>, glue_times_vec>(X.get_ref(), Y);
  }



//! Base * rowvec
template<typename T1>
arma_inline
const Glue<T1, Row<typename T1::elem_type>, glue_times_vec>
operator*
(const Base<typename T1::elem_type,T1>& X, const Row<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Row<typename T1::elem_type>, glue_times_vec>(X.get_ref(), Y);
  }



//! colvec * Base
template<typename T1>
arma_inline
const Glue<Col<typename T1::elem_type>, T1, glue_times_vec>
operator*
(const Col<typename T1::elem_type>& X, const Base<typename T1::elem_type,T1>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Col<typename T1::elem_type>, T1, glue_times_vec>(X, Y.get_ref());
  }



//! rowvec * Base
template<typename T1>
arma_inline
const Glue<Row<typename T1::elem_type>, T1, glue_times_vec>
operator*
(const Row<typename T1::elem_type>& X, const Base<typename T1::elem_type,T1>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Row<typename T1::elem_type>, T1, glue_times_vec>(X, Y.get_ref());
  }



//! diagmat * colvec
template<typename T1>
arma_inline
const Glue<Op<T1, op_diagmat>, Col<typename T1::elem_type>, glue_times_diag>
operator*
(const Op<T1, op_diagmat>& X, const Col<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1, op_diagmat>, Col<typename T1::elem_type>, glue_times_diag>(X, Y);
  }



//! diagmat * rowvec
template<typename T1>
arma_inline
const Glue<Op<T1, op_diagmat>, Row<typename T1::elem_type>, glue_times_diag>
operator*
(const Op<T1, op_diagmat>& X, const Row<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1, op_diagmat>, Row<typename T1::elem_type>, glue_times_diag>(X, Y);
  }



//! colvec * diagmat
template<typename T1>
arma_inline
const Glue<Col<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>
operator*
(const Col<typename T1::elem_type>& X, const Op<T1, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Col<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>(X, Y);
  }



//! rowvec * diagmat
template<typename T1>
arma_inline
const Glue<Row<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>
operator*
(const Row<typename T1::elem_type>& X, const Op<T1, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Row<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>(X, Y);
  }



//
// multiplication of Base objects with different element types
//



//! Base * Base
template<typename eT1, typename T1, typename eT2, typename T2>
arma_inline
Mat<typename promote_type<eT1,eT2>::result>
operator*
(const Base<eT1,T1>& X, const Base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();

  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat< eT1 >& A = tmp1.M;
  const Mat< eT2 >& B = tmp2.M;
  
  Mat< typename promote_type<eT1,eT2>::result > out;
  
  glue_times::apply_mixed(out, A, B);
  
  return out;
  }



//
// multiplication of Base objects with same element types
//



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<std::complex<double>,T1>& X, const Base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<std::complex<float>,T1>& X, const Base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<double,T1>& X, const Base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<float,T1>& X, const Base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<s32,T1>& X, const Base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<u32,T1>& X, const Base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<s16,T1>& X, const Base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<u16,T1>& X, const Base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<s8,T1>& X, const Base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<u8,T1>& X, const Base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



//! @}
