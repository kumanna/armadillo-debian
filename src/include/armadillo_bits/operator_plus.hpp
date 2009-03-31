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


//! \addtogroup operator_plus
//! @{



//! unary plus operation (does nothing, but is required for completeness)
template<typename T1>
inline
const arma_base<typename T1::elem_type,T1>&
operator+
(const arma_base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return X;
  }



//! arma_base + scalar
template<typename T1>
inline
const op_data<T1, op_scalar_plus>
operator+
(const arma_base<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return op_data<T1, op_scalar_plus>(X.get_ref(), k);
  }



//! op + scalar, level 2
template<typename T1>
inline
const op_data<T1,op_scalar_plus>
operator+
(const op_data<T1,op_scalar_plus>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return op_data<T1, op_scalar_plus>(X.m, X.aux + k);
  }



//! scalar + arma_base
template<typename T1>
inline
const op_data<T1, op_scalar_plus>
operator+
(const typename T1::elem_type k, const arma_base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return op_data<T1, op_scalar_plus>(X.get_ref(), k);  // NOTE: order is swapped
  }



//! arma_base + diagmat
template<typename T1, typename T2>
inline
const glue_data<T1, op_data<T2,op_diagmat>, glue_plus_diag>
operator+
(const arma_base<typename T2::elem_type, T1>& X, const op_data<T2,op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, op_data<T2,op_diagmat>, glue_plus_diag>(X.get_ref(), Y);
  }



//! diagmat + arma_base
template<typename T1, typename T2>
inline
const glue_data<T1, op_data<T2,op_diagmat>, glue_plus_diag>
operator+
(const op_data<T2,op_diagmat>& Y, const arma_base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, op_data<T2,op_diagmat>, glue_plus_diag>(X.get_ref(), Y);  // NOTE: order is swapped
  }



//
// addition of arma_base objects with different element types
//



//! arma_base + arma_base
template<typename eT1, typename T1, typename eT2, typename T2>
inline
basic_mat<typename promote_type<eT1,eT2>::result>
operator+
(const arma_base<eT1,T1>& X, const arma_base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();

  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const basic_mat<eT1>& A = tmp1.M;
  const basic_mat<eT2>& B = tmp2.M;
  
  basic_mat< typename promote_type<eT1,eT2>::result > out;
  glue_plus::apply_mixed(out, A, B);
  
  return out;
  }



//
// addition of arma_base objects with same element types
//



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<std::complex<double>,T1>& X, const arma_base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<std::complex<float>,T1>& X, const arma_base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<double,T1>& X, const arma_base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<float,T1>& X, const arma_base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<s32,T1>& X, const arma_base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<u32,T1>& X, const arma_base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<s16,T1>& X, const arma_base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<u16,T1>& X, const arma_base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<s8,T1>& X, const arma_base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_plus>
operator+
(const arma_base<u8,T1>& X, const arma_base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_plus>(X.get_ref(), Y.get_ref());
  }



//
// old operators
// 

// //! unary plus operation (does nothing, but is required for completeness)
// template<typename mT>
// inline
// const basic_mat<mT>&
// operator+
// (const basic_mat<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return X;
//   }
// 
// 
// 
// //! unary plus operation (does nothing, but is required for completeness)
// template<typename T1, typename op_type>
// inline
// const op_data<T1, op_type>&
// operator+
// (const op_data<T1, op_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return X;
//   }
// 
// //! unary plus operation (does nothing, but is required for completeness)
// template<typename T1, typename T2, typename glue_type>
// inline
// const glue_data<T1, T2, glue_type>&
// operator+
// (const glue_data<T1, T2, glue_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return X;
//   }
// 
// 
// 
// //! mat + scalar
// template<typename mT>
// inline
// const op_data<basic_mat<mT>, op_scalar_plus>
// operator+
// (const basic_mat<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<basic_mat<mT>, op_scalar_plus>(X,k);
//   }
// 
// 
// 
// //! op + scalar
// template<typename T1, typename op_type>
// inline
// const op_data<op_data<T1,op_type>, op_scalar_plus>
// operator+
// (const op_data<T1,op_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<op_data<T1,op_type>, op_scalar_plus>(X,k);
//   }
// 
// 
// 
// //! op + scalar, level 2
// template<typename T1>
// inline
// const op_data<T1,op_scalar_plus>
// operator+
// (const op_data<T1,op_scalar_plus>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<T1, op_scalar_plus>(X.m, X.aux + k);
//   }
// 
// 
// //! glue + scalar
// template<typename T1, typename T2, typename glue_type>
// inline
// const op_data<glue_data<T1,T2,glue_type>, op_scalar_plus>
// operator+
// (const glue_data<T1,T2,glue_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<glue_data<T1,T2,glue_type>, op_scalar_plus>(X,k);
//   }
// 
// 
// //! scalar + mat
// template<typename mT>
// inline
// const op_data<basic_mat<mT>, op_scalar_plus>
// operator+
// (const mT k, const basic_mat<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<basic_mat<mT>, op_scalar_plus>(X,k);  // NOTE: order is swapped
//   }
// 
// 
// 
// //! scalar + op
// template<typename T1, typename op_type>
// inline
// const op_data<op_data<T1,op_type>, op_scalar_plus>
// operator+
// (const typename T1::elem_type k, const op_data<T1,op_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<op_data<T1,op_type>, op_scalar_plus>(X,k);  // NOTE: order is swapped
//   }
// 
// 
// 
// //! scalar + glue
// template<typename T1, typename T2, typename glue_type>
// inline
// const op_data<glue_data<T1,T2,glue_type>, op_scalar_plus>
// operator+
// (const typename T1::elem_type k, const glue_data<T1,T2,glue_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<glue_data<T1,T2,glue_type>, op_scalar_plus>(X,k);  // NOTE: order is swapped
//   }
// 
// 
// //! mat + mat
// template<typename mT>
// inline
// const glue_data<basic_mat<mT>, basic_mat<mT>, glue_plus>
// operator+
// (const basic_mat<mT>& X, const basic_mat<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<mT>, basic_mat<mT>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! mat + diagmat(T1)
// template<typename T1>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_diagmat>, glue_plus_diag>
// operator+
// (const basic_mat<typename T1::elem_type>& X, const op_data<T1,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_diagmat>, glue_plus_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) + mat
// template<typename T1>
// inline
// const glue_data<op_data<T1,op_diagmat>, basic_mat<typename T1::elem_type>, glue_plus_diag>
// operator+
// (const op_data<T1,op_diagmat>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_diagmat>, basic_mat<typename T1::elem_type>, glue_plus_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) + diagmat(T2)
// template<typename T1, typename T2>
// inline
// const glue_data< op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_plus_diag>
// operator+
// (const op_data<T1,op_diagmat>& X, const op_data<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_plus_diag>(X,Y);
//   }
// 
// 
// 
// //! mat + op
// template<typename T1, typename op_type>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_type>, glue_plus>
// operator+
// (const basic_mat<typename T1::elem_type>& X, const op_data<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_type>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) + op
// template<typename T1, typename T2, typename op_type>
// inline
// const glue_data< op_data<T1,op_diagmat>, op_data<T2,op_type>, glue_plus_diag>
// operator+
// (const op_data<T1,op_diagmat>& X, const op_data<T2,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, op_data<T2,op_type>, glue_plus_diag>(X,Y);
//   }
// 
// 
// //! op + mat
// template<typename T1, typename op_type>
// inline
// const glue_data<op_data<T1,op_type>, basic_mat<typename T1::elem_type>, glue_plus>
// operator+
// (const op_data<T1,op_type>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, basic_mat<typename T1::elem_type>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! op + diagmat(T2)
// template<typename T1, typename op_type, typename T2>
// inline
// const glue_data<op_data<T1,op_type>, op_data<T2,op_diagmat>, glue_plus_diag>
// operator+
// (const op_data<T1,op_type>& X, const op_data<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, op_data<T2,op_diagmat>, glue_plus_diag>(X,Y);
//   }
// 
// 
// 
// //! op + glue
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// inline
// const glue_data<op_data<T1,op_type>, glue_data<T2, T3, glue_type>, glue_plus>
// operator+
// (const op_data<T1,op_type>& X, const glue_data<T2, T3, glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, glue_data<T2, T3, glue_type>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! glue + op
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// inline
// const glue_data<glue_data<T2, T3, glue_type>, op_data<T1,op_type>, glue_plus>
// operator+
// (const glue_data<T2, T3, glue_type>& X, const op_data<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T2, T3, glue_type>, op_data<T1,op_type>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! mat + glue
// template<typename T1, typename T2, typename glue_type>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, glue_data<T1, T2,glue_type>, glue_plus>
// operator+
// (const basic_mat<typename T1::elem_type>& X, const glue_data<T1, T2,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, glue_data<T1, T2,glue_type>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) + glue
// template<typename T1, typename T2, typename T3, typename glue_type>
// inline
// const glue_data< op_data<T1,op_diagmat>, glue_data<T2,T3,glue_type>, glue_plus_diag>
// operator+
// (const op_data<T1,op_diagmat>& X, const glue_data<T2,T3,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, glue_data<T2,T3,glue_type>, glue_plus_diag>(X,Y);
//   }
// 
// 
// 
// //! glue + mat
// template<typename T1, typename T2, typename glue_type>
// inline
// const glue_data<glue_data<T1, T2,glue_type>, basic_mat<typename T1::elem_type>, glue_plus>
// operator+
// (const glue_data<T1, T2,glue_type>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T1, T2,glue_type>, basic_mat<typename T1::elem_type>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! glue + diagmat(T3)
// template<typename T1, typename T2, typename glue_type, typename T3>
// inline
// const glue_data< glue_data<T1,T2,glue_type>, op_data<T3,op_diagmat>, glue_plus_diag>
// operator+
// (const glue_data<T1, T2,glue_type>& X, const op_data<T3,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< glue_data<T1,T2,glue_type>, op_data<T3,op_diagmat>, glue_plus_diag>(X,Y);
//   }
// 
// 
// 
// //! op + op
// template<typename T1, typename op_type1, typename T2, typename op_type2>
// inline
// const glue_data<op_data<T1,op_type1>, op_data<T2,op_type2>, glue_plus>
// operator+
// (const op_data<T1,op_type1>& X, const op_data<T2,op_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type1>, op_data<T2,op_type2>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! glue + glue
// template<typename T1, typename T2, typename glue_type1, typename T3, typename T4, typename glue_type2>
// inline
// const glue_data<glue_data<T1,T2,glue_type1>, glue_data<T3,T4,glue_type2>, glue_plus>
// operator+
// (const glue_data<T1,T2,glue_type1>& X, const glue_data<T3,T4,glue_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T1,T2,glue_type1>, glue_data<T3,T4,glue_type2>, glue_plus>(X,Y);
//   }
// 
// 
// 
// //! arma_base + subview
// template<typename T1>
// inline
// const glue_data<T1, subview<typename T1::elem_type>, glue_plus>
// operator+
// (const arma_base<T1>& X, const subview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<T1, subview<typename T1::elem_type>, glue_plus>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! subview + arma_base
// template<typename T1>
// inline
// const glue_data<subview<typename T1::elem_type>, T1, glue_plus>
// operator+
// (const subview<typename T1::elem_type>& X, const arma_base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<subview<typename T1::elem_type>, T1, glue_plus>(X,Y.get_ref());
//   }
// 
// 
// 
// //! arma_base + diagview
// template<typename T1>
// inline
// const glue_data<T1, diagview<typename T1::elem_type>, glue_plus>
// operator+
// (const arma_base<T1>& X, const diagview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<T1, diagview<typename T1::elem_type>, glue_plus>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! diagview + arma_base
// template<typename T1>
// inline
// const glue_data<diagview<typename T1::elem_type>, T1, glue_plus>
// operator+
// (const diagview<typename T1::elem_type>& X, const T1& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<diagview<typename T1::elem_type>, T1, glue_plus>(X,Y.get_ref());
//   }
// 
// 
// 
// //! scalar + subview
// template<typename mT>
// inline
// const op_data<subview<mT>, op_scalar_plus>
// operator+
// (const mT k, const subview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<subview<mT>, op_scalar_plus>(X,k);
//   }
// 
// 
// 
// //! scalar + diagview
// template<typename mT>
// inline
// const op_data<diagview<mT>, op_scalar_plus>
// operator+
// (const mT k, const diagview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<diagview<mT>, op_scalar_plus>(X,k);
//   }
// 
// 
// 
// //! subview + scalar
// template<typename mT>
// inline
// const op_data<subview<mT>, op_scalar_plus>
// operator+
// (const subview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<subview<mT>, op_scalar_plus>(X,k);
//   }
// 
// 
// 
// //! diagview + scalar
// template<typename mT>
// inline
// const op_data<diagview<mT>, op_scalar_plus>
// operator+
// (const diagview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return op_data<diagview<mT>, op_scalar_plus>(X,k);
//   }



//! @}
