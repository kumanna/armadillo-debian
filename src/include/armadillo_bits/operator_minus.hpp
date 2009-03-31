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


//! \addtogroup operator_minus
//! @{



//! unary -
template<typename T1>
inline
const op_data<T1, op_neg>
operator-
(const arma_base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return op_data<T1,op_neg>(X.get_ref());
  }



//! cancellation of two consecutive negations: -(-T1)
template<typename T1>
inline
const T1&
operator-
(const op_data<T1, op_neg>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.m;
  }



//! arma_base - scalar
template<typename T1>
inline
const op_data<T1, op_scalar_minus_post>
operator-
(const arma_base<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return op_data<T1, op_scalar_minus_post>(X.get_ref(), k);
  }



//! scalar - arma_base
template<typename T1>
inline
const op_data<T1, op_scalar_minus_pre>
operator-
(const typename T1::elem_type k, const arma_base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return op_data<T1, op_scalar_minus_pre>(X.get_ref(), k);
  }



//! arma_base - - arma_base = arma_base + arma_base
template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator-
  (
  const arma_base<typename T1::elem_type, T1                 >& X,
  const arma_base<typename T1::elem_type, op_data<T2,op_neg> >& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const T1& A = X.get_ref();
  const T2& B = (Y.get_ref()).m;
  
  return glue_data<T1, T2, glue_plus>(A,B);
  }



//! arma_base - diagmat
template<typename T1, typename T2>
inline
const glue_data<T1, op_data<T2,op_diagmat>, glue_minus_diag>
operator-
(const arma_base<typename T2::elem_type,T1>& X, const op_data<T2,op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, op_data<T1,op_diagmat>, glue_minus_diag>(X.get_ref(), Y);
  }



//! diagmat - arma_base
template<typename T1, typename T2>
inline
const glue_data< op_data<T1,op_diagmat>, T2, glue_minus_diag>
operator-
(const op_data<T1,op_diagmat>& X, const arma_base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data< op_data<T1,op_diagmat>, T2, glue_minus_diag>(X, Y.get_ref());
  }



//! arma_base - op_data<T2,op_neg> = arma_base + T2
template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator-
(const arma_base<typename T2::elem_type,T1>& X, const op_data<T2, op_neg>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.m);
  }



//! diagmat - op_data<T2,op_neg> = diagmat + T2
template<typename T1, typename T2>
inline
const glue_data< op_data<T1,op_diagmat>, T2, glue_plus_diag>
operator-
  (
  const arma_base<typename T1::elem_type, op_data<T1,op_diagmat> >& X,
  const arma_base<typename T1::elem_type, op_data<T2,op_neg    > >& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return glue_data< op_data<T1,op_diagmat>, T2, glue_plus_diag>(X.get_ref(), (Y.get_ref()).m);
  }



//
// subtraction of arma_base objects with different element types
//



//! arma_base - arma_base
template<typename eT1, typename T1, typename eT2, typename T2>
inline
basic_mat<typename promote_type<eT1,eT2>::result>
operator-
(const arma_base<eT1,T1>& X, const arma_base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const basic_mat<eT1>& A = tmp1.M;
  const basic_mat<eT2>& B = tmp2.M;
  
  basic_mat< typename promote_type<eT1,eT2>::result > out;
  glue_minus::apply_mixed(out, A, B);
  
  return out;
  }



//
// subtraction of arma_base objects with same element types
//



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<std::complex<double>,T1>& X, const arma_base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<std::complex<float>,T1>& X, const arma_base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<double,T1>& X, const arma_base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<float,T1>& X, const arma_base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<s32,T1>& X, const arma_base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<u32,T1>& X, const arma_base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<s16,T1>& X, const arma_base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<u16,T1>& X, const arma_base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<s8,T1>& X, const arma_base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_minus>
operator-
(const arma_base<u8,T1>& X, const arma_base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



//
// old operators
//


// //! unary -
// template<typename mT>
// inline
// const op_data<basic_mat<mT>, op_neg>
// operator-
// (const basic_mat<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<basic_mat<mT>,op_neg>(X);
//   }
// 
// 
// 
// //! unary -
// template<typename T1, typename op_type>
// inline
// const op_data< op_data<T1, op_type>, op_neg >
// operator-
// (const op_data<T1, op_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<op_data<T1, op_type>,op_neg>(X);
//   }
// 
// 
// 
// //! cancellation of two consecutive negations: -(-T1)
// template<typename T1>
// inline
// const T1&
// operator-
// (const op_data<T1, op_neg>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return X.m;
//   }
// 
// 
// 
// //! -(glue_data<T1,T2,glue_type>)
// template<typename T1, typename T2, typename glue_type>
// inline
// const op_data< glue_data<T1,T2,glue_type>, op_neg >
// operator-
// (const glue_data<T1,T2,glue_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data< glue_data<T1,T2,glue_type>,op_neg >(X);
//   }
// 
// 
// 
// //! mat - scalar
// template<typename mT>
// inline
// const op_data<basic_mat<mT>, op_scalar_minus_post>
// operator-
// (const basic_mat<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<basic_mat<mT>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! op - scalar
// template<typename T1, typename op_type>
// inline
// const op_data<op_data<T1,op_type>, op_scalar_minus_post>
// operator-
// (const op_data<T1,op_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<op_data<T1,op_type>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! op - scalar, level 2
// template<typename T1>
// inline
// const op_data<T1, op_scalar_minus_post>
// operator-
// (const op_data<T1,op_scalar_plus>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<T1,op_scalar_minus_post>(X.m, X.aux - k);
//   }
// 
// 
// 
// //! glue - scalar
// template<typename T1, typename T2, typename glue_type>
// inline
// const op_data<glue_data<T1,T2,glue_type>, op_scalar_minus_post>
// operator-
// (const glue_data<T1,T2,glue_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<glue_data<T1,T2,glue_type>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! scalar - mat
// template<typename mT>
// inline
// const op_data<basic_mat<mT>, op_scalar_minus_pre>
// operator-
// (const mT k, const basic_mat<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<basic_mat<mT>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! scalar - op
// template<typename T1, typename op_type>
// inline
// const op_data<op_data<T1,op_type>, op_scalar_minus_pre>
// operator-
// (const typename T1::elem_type k, const op_data<T1,op_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<op_data<T1,op_type>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! scalar - glue
// template<typename T1, typename T2, typename glue_type>
// inline
// const op_data<glue_data<T1,T2,glue_type>, op_scalar_minus_pre>
// operator-
// (const typename T1::elem_type k, const glue_data<T1,T2,glue_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<glue_data<T1,T2,glue_type>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! mat - mat
// template<typename mT>
// inline
// const glue_data<basic_mat<mT>, basic_mat<mT>, glue_minus>
// operator-
// (const basic_mat<mT>& X, const basic_mat<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<mT>, basic_mat<mT>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! mat - diagmat(T1)
// template<typename T1>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_diagmat>, glue_minus_diag>
// operator-
// (const basic_mat<typename T1::elem_type>& X, const op_data<T1,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - mat
// template<typename T1>
// inline
// const glue_data< op_data<T1,op_diagmat>, basic_mat<typename T1::elem_type>, glue_minus_diag>
// operator-
// (const op_data<T1,op_diagmat>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, basic_mat<typename T1::elem_type>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - diagmat(T2)
// template<typename T1, typename T2>
// inline
// const glue_data< op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_minus_diag>
// operator-
// (const op_data<T1,op_diagmat>& X, const op_data<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! mat - op
// template<typename T1, typename op_type>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_type>, glue_minus>
// operator-
// (const basic_mat<typename T1::elem_type>& X, const op_data<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - op
// template<typename T1, typename T2, typename op_type>
// inline
// const glue_data< op_data<T1,op_diagmat>, op_data<T2,op_type>, glue_minus_diag>
// operator-
// (const op_data<T1,op_diagmat>& X, const op_data<T2,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, op_data<T2,op_type>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! op - mat
// template<typename T1, typename op_type>
// inline
// const glue_data<op_data<T1,op_type>, basic_mat<typename T1::elem_type>, glue_minus>
// operator-
// (const op_data<T1,op_type>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, basic_mat<typename T1::elem_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! op - diagmat(T2)
// template<typename T1, typename op_type, typename T2>
// inline
// const glue_data<op_data<T1,op_type>, op_data<T2,op_diagmat>, glue_minus_diag>
// operator-
// (const op_data<T1,op_type>& X, const op_data<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, op_data<T2,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! op - glue
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// inline
// const glue_data<op_data<T1,op_type>, glue_data<T2, T3, glue_type>, glue_minus>
// operator-
// (const op_data<T1,op_type>& X, const glue_data<T2, T3, glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, glue_data<T2, T3, glue_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! glue - op
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// inline
// const glue_data<glue_data<T2, T3, glue_type>, op_data<T1,op_type>, glue_minus>
// operator-
// (const glue_data<T2, T3, glue_type>& X, const op_data<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T2, T3, glue_type>, op_data<T1,op_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! mat - glue
// template<typename T1, typename T2, typename glue_type>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, glue_data<T1, T2,glue_type>, glue_minus>
// operator-
// (const basic_mat<typename T1::elem_type>& X, const glue_data<T1, T2,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, glue_data<T1, T2,glue_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - glue
// template<typename T1, typename T2, typename T3, typename glue_type>
// inline
// const glue_data<op_data<T1,op_diagmat>, glue_data<T2,T3,glue_type>, glue_minus_diag>
// operator-
// (const op_data<T1,op_diagmat>& X, const glue_data<T2,T3,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_diagmat>, glue_data<T2,T3,glue_type>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! glue - mat
// template<typename T1, typename T2, typename glue_type>
// inline
// const glue_data< glue_data<T1,T2,glue_type>, basic_mat<typename T1::elem_type>, glue_minus>
// operator-
// (const glue_data<T1,T2,glue_type>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< glue_data<T1,T2,glue_type>, basic_mat<typename T1::elem_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! glue - diagmat(T1)
// template<typename T1, typename T2, typename glue_type, typename T3>
// inline
// const glue_data< glue_data<T1,T2,glue_type>, op_data<T3,op_diagmat>, glue_minus_diag>
// operator-
// (const glue_data<T1,T2,glue_type>& X, const op_data<T3,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< glue_data<T1,T2,glue_type>, op_data<T3,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! op - op
// template<typename T1, typename op_type1, typename T2, typename op_type2>
// inline
// const glue_data<op_data<T1,op_type1>, op_data<T2,op_type2>, glue_minus>
// operator-
// (const op_data<T1,op_type1>& X, const op_data<T2,op_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type1>, op_data<T2,op_type2>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! glue - glue
// template<typename T1, typename T2, typename glue_type1, typename T3, typename T4, typename glue_type2>
// inline
// const glue_data<glue_data<T1,T2,glue_type1>, glue_data<T3,T4,glue_type2>, glue_minus>
// operator-
// (const glue_data<T1,T2,glue_type1>& X, const glue_data<T3,T4,glue_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T1,T2,glue_type1>, glue_data<T3,T4,glue_type2>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! mat - op_data<T2,op_neg> = mat + T2
// template<typename T2>
// inline
// const glue_data<basic_mat<typename T2::elem_type>, T2, glue_plus>
// operator-
// (const basic_mat<typename T2::elem_type>& X, const op_data<T2, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T2::elem_type>, T2, glue_plus>(X,Y.m);
//   }
// 
// 
// 
// //! diagmat - op_data<T2,op_neg> = diagmat + T2
// template<typename T1, typename T2>
// inline
// const glue_data< op_data<T1,op_diagmat>, T2, glue_plus_diag>
// operator-
// (const op_data<T1,op_diagmat>& X, const op_data<T2, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, T2, glue_plus_diag>(X,Y.m);
//   }
// 
// 
// 
// //! op - op_data<T2,op_neg> = op + T2
// template<typename T1, typename op_type, typename T2>
// inline
// const glue_data<op_data<T1, op_type>, T2, glue_plus>
// operator-
// (const op_data<T1, op_type>& X, const op_data<T2, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1, op_type>, T2, glue_plus>(X,Y.m);
//   }
// 
// 
// 
// //! glue - op_data<T3,op_neg> = glue + T3
// template<typename T1, typename T2, typename glue_type, typename T3>
// inline
// const glue_data<glue_data<T1,T2,glue_type>, T3, glue_plus>
// operator-
// (const glue_data<T1,T2,glue_type>& X, const op_data<T3, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< glue_data<T1,T2,glue_type>, T3, glue_plus>(X,Y.m);
//   }
// 
// 
// 
// //! arma_base - subview
// template<typename T1>
// inline
// const glue_data<T1, subview<typename T1::elem_type>, glue_minus>
// operator-
// (const arma_base<T1>& X, const subview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<T1, subview<typename T1::elem_type>, glue_minus>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! subview - arma_base
// template<typename T1>
// inline
// const glue_data<subview<typename T1::elem_type>, T1, glue_minus>
// operator-
// (const subview<typename T1::elem_type>& X, const arma_base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<subview<typename T1::elem_type>, T1, glue_minus>(X,Y.get_ref());
//   }
// 
// 
// 
// //! arma_base - diagview
// template<typename T1>
// inline
// const glue_data<T1, diagview<typename T1::elem_type>, glue_minus>
// operator-
// (const arma_base<T1>& X, const diagview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<T1, diagview<typename T1::elem_type>, glue_minus>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! diagview - arma_base
// template<typename T1>
// inline
// const glue_data<diagview<typename T1::elem_type>, T1, glue_minus>
// operator-
// (const diagview<typename T1::elem_type>& X, const arma_base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<diagview<typename T1::elem_type>, T1, glue_minus>(X,Y.get_ref());
//   }
// 
// 
// 
// //! scalar - subview
// template<typename mT>
// inline
// const op_data<subview<mT>, op_scalar_minus_pre>
// operator-
// (const mT k, const subview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<subview<mT>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! scalar - diagview
// template<typename mT>
// inline
// const op_data<diagview<mT>, op_scalar_minus_pre>
// operator-
// (const mT k, const diagview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<diagview<mT>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! subview - scalar
// template<typename mT>
// inline
// const op_data<subview<mT>, op_scalar_minus_post>
// operator-
// (const subview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<subview<mT>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! diagview - scalar
// template<typename mT>
// inline
// const op_data<diagview<mT>, op_scalar_minus_post>
// operator-
// (const diagview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<diagview<mT>, op_scalar_minus_post>(X,k);
//   }



//! @}
