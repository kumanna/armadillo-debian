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


//! \addtogroup operator_schur
//! @{


// operator %, which we define it to do a schur product (element-wise multiplication)


//! Base % Base
template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<typename T1::elem_type,T1>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



//! Base % diagmat
template<typename T1, typename T2>
arma_inline
const Glue<T1, Op<T2,op_diagmat>, glue_schur_diag>
operator%
(const Base<typename T2::elem_type,T1>& X, const Op<T2,op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Op<T2,op_diagmat>, glue_schur_diag>(X.get_ref(), Y);
  }



//! diagmat % Base
template<typename T1, typename T2>
arma_inline
const Glue< Op<T1,op_diagmat>, T2, glue_schur_diag>
operator%
(const Op<T1,op_diagmat>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue< Op<T1,op_diagmat>, T2, glue_schur_diag>(X, Y.get_ref());
  }



//! diagmat % diagmat
template<typename T1, typename T2>
inline
Mat< typename promote_type<typename T1::elem_type, typename T2::elem_type>::result >
operator%
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
  
  arma_debug_assert_same_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "element-wise matrix multiplication");
  
  Mat<out_eT> out(A.n_rows, A.n_rows);
  out.zeros();
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    out.at(i,i) = upgrade_val<eT1,eT2>::apply( A.at(i,i) ) * upgrade_val<eT1,eT2>::apply( B.at(i,i) );
    }
  
  return out;
  }



//
// schur product of Base objects with different element types
//



//! Base % Base
template<typename eT1, typename T1, typename eT2, typename T2>
arma_inline
Mat<typename promote_type<eT1,eT2>::result>
operator%
(const Base<eT1,T1>& X, const Base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT1>& A = tmp1.M;
  const Mat<eT2>& B = tmp2.M;
  
  Mat< typename promote_type<eT1,eT2>::result > out;
  glue_schur::apply_mixed(out, A, B);
  
  return out;
  }



//
// schur product of Base objects with same element types
//



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<std::complex<double>,T1>& X, const Base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<std::complex<float>,T1>& X, const Base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<double,T1>& X, const Base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<float,T1>& X, const Base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<s32,T1>& X, const Base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<u32,T1>& X, const Base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<s16,T1>& X, const Base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<u16,T1>& X, const Base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<s8,T1>& X, const Base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_schur>
operator%
(const Base<u8,T1>& X, const Base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



//! @}
