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


//! arma_base % arma_base
template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<typename T1::elem_type,T1>& X, const arma_base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



//! arma_base % diagmat
template<typename T1, typename T2>
inline
const glue_data<T1, op_data<T2,op_diagmat>, glue_schur_diag>
operator%
(const arma_base<typename T2::elem_type,T1>& X, const op_data<T2,op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, op_data<T2,op_diagmat>, glue_schur_diag>(X.get_ref(), Y);
  }



//! diagmat % arma_base
template<typename T1, typename T2>
inline
const glue_data< op_data<T1,op_diagmat>, T2, glue_schur_diag>
operator%
(const op_data<T1,op_diagmat>& X, const arma_base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data< op_data<T1,op_diagmat>, T2, glue_schur_diag>(X, Y.get_ref());
  }



//
// schur product of arma_base objects with different element types
//



//! arma_base % arma_base
template<typename eT1, typename T1, typename eT2, typename T2>
inline
basic_mat<typename promote_type<eT1,eT2>::result>
operator%
(const arma_base<eT1,T1>& X, const arma_base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const basic_mat<eT1>& A = tmp1.M;
  const basic_mat<eT2>& B = tmp2.M;
  
  basic_mat< typename promote_type<eT1,eT2>::result > out;
  glue_schur::apply_mixed(out, A, B);
  
  return out;
  }



//
// schur product of arma_base objects with same element types
//



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<std::complex<double>,T1>& X, const arma_base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<std::complex<float>,T1>& X, const arma_base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<double,T1>& X, const arma_base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<float,T1>& X, const arma_base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<s32,T1>& X, const arma_base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<u32,T1>& X, const arma_base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<s16,T1>& X, const arma_base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<u16,T1>& X, const arma_base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<s8,T1>& X, const arma_base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
inline
const glue_data<T1, T2, glue_schur>
operator%
(const arma_base<u8,T1>& X, const arma_base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return glue_data<T1, T2, glue_schur>(X.get_ref(), Y.get_ref());
  }



//
// old operators
//


// //! mat % mat
// template<typename mT>
// inline
// const glue_data<basic_mat<mT>, basic_mat<mT>, glue_schur>
// operator%
// (const basic_mat<mT>& X, const basic_mat<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<mT>, basic_mat<mT>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! mat % diagmat(T1)
// template<typename T1>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_diagmat>, glue_schur_diag>
// operator%
// (const basic_mat<typename T1::elem_type>& X, const op_data<T1,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_diagmat>, glue_schur_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) % mat
// template<typename mT, typename T1>
// inline
// const glue_data< op_data<T1,op_diagmat>, basic_mat<typename T1::elem_type>, glue_schur_diag>
// operator%
// (const op_data<T1,op_diagmat>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, basic_mat<typename T1::elem_type>, glue_schur_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) % diagmat(T2)
// template<typename T1, typename T2>
// inline
// const glue_data< op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_schur_diag>
// operator%
// (const op_data<T1,op_diagmat>& X, const op_data<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, op_data<T2,op_diagmat>, glue_schur_diag>(X,Y);
//   }
// 
// 
// 
// //! mat % op
// template<typename T1, typename op_type>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_type>, glue_schur>
// operator%
// (const basic_mat<typename T1::elem_type>& X, const op_data<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, op_data<T1,op_type>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! op % mat
// template<typename T1, typename op_type>
// inline
// const glue_data<op_data<T1,op_type>, basic_mat<typename T1::elem_type>, glue_schur>
// operator%
// (const op_data<T1,op_type>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, basic_mat<typename T1::elem_type>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! op % diagmat(T2)
// template<typename T1, typename op_type, typename T2>
// inline
// const glue_data<op_data<T1,op_type>, op_data<T2,op_diagmat>, glue_schur_diag>
// operator%
// (const op_data<T1,op_type>& X, const op_data<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, op_data<T2,op_diagmat>, glue_schur_diag>(X,Y);
//   }
// 
// 
// 
// //! op % glue
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// inline
// const glue_data<op_data<T1,op_type>, glue_data<T2, T3, glue_type>, glue_schur>
// operator%
// (const op_data<T1,op_type>& X, const glue_data<T2, T3, glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type>, glue_data<T2, T3, glue_type>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! glue % op
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// inline
// const glue_data<glue_data<T2, T3, glue_type>, op_data<T1,op_type>, glue_schur>
// operator%
// (const glue_data<T2, T3, glue_type>& X, const op_data<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T2, T3, glue_type>, op_data<T1,op_type>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! mat % glue
// template<typename T1, typename T2, typename glue_type>
// inline
// const glue_data<basic_mat<typename T1::elem_type>, glue_data<T1,T2,glue_type>, glue_schur>
// operator%
// (const basic_mat<typename T1::elem_type>& X, const glue_data<T1,T2,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<basic_mat<typename T1::elem_type>, glue_data<T1,T2,glue_type>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) % glue
// template<typename T1, typename T2, typename T3, typename glue_type>
// inline
// const glue_data< op_data<T1,op_diagmat>, glue_data<T2,T3,glue_type>, glue_schur_diag>
// operator%
// (const op_data<T1,op_diagmat>& X, const glue_data<T2,T3,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data< op_data<T1,op_diagmat>, glue_data<T2,T3,glue_type>, glue_schur_diag>(X,Y);
//   }
// 
// 
// 
// //! glue % mat
// template<typename T1, typename T2, typename glue_type>
// inline
// const glue_data<glue_data<T1, T2,glue_type>, basic_mat<typename T1::elem_type>, glue_schur>
// operator%
// (const glue_data<T1,T2,glue_type>& X, const basic_mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T1,T2,glue_type>, basic_mat<typename T1::elem_type>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! glue % diagmat(T3)
// template<typename T1, typename T2, typename glue_type, typename T3>
// inline
// const glue_data<glue_data<T1, T2,glue_type>, op_data<T3,op_diagmat>, glue_schur_diag>
// operator%
// (const glue_data<T1,T2,glue_type>& X, const op_data<T3,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T1,T2,glue_type>, op_data<T3,op_diagmat>, glue_schur_diag>(X,Y);
//   }
// 
// 
// 
// //! op % op
// template<typename T1, typename op_type1, typename T2, typename op_type2>
// inline
// const glue_data<op_data<T1,op_type1>, op_data<T2,op_type2>, glue_schur>
// operator%
// (const op_data<T1,op_type1>& X, const op_data<T2,op_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<op_data<T1,op_type1>, op_data<T2,op_type2>, glue_schur>(X,Y);
//   }
// 
// 
// //! glue % glue
// template<typename T1, typename T2, typename glue_type1, typename T3, typename T4, typename glue_type2>
// inline
// const glue_data<glue_data<T1,T2,glue_type1>, glue_data<T3,T4,glue_type2>, glue_schur>
// operator%
// (const glue_data<T1,T2,glue_type1>& X, const glue_data<T3,T4,glue_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<glue_data<T1,T2,glue_type1>, glue_data<T3,T4,glue_type2>, glue_schur>(X,Y);
//   }
// 
// 
// 
// //! arma_base % subview
// template<typename T1>
// inline
// const glue_data<T1, subview<typename T1::elem_type>, glue_schur>
// operator%
// (const arma_base<T1>& X, const subview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<T1, subview<typename T1::elem_type>, glue_schur>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! subview % arma_base
// template<typename T1>
// inline
// const glue_data<subview<typename T1::elem_type>, T1, glue_schur>
// operator%
// (const subview<typename T1::elem_type>& X, const arma_base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<subview<typename T1::elem_type>, T1, glue_schur>(X,Y.get_ref());
//   }
// 
// 
// 
// //! arma_base % diagview
// template<typename T1>
// inline
// const glue_data<T1, diagview<typename T1::elem_type>, glue_schur>
// operator%
// (const arma_base<T1>& X, const diagview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<T1, diagview<typename T1::elem_type>, glue_schur>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! diagview % arma_base
// template<typename T1>
// inline
// const glue_data<diagview<typename T1::elem_type>, T1, glue_schur>
// operator%
// (const diagview<typename T1::elem_type>& X, const arma_base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return glue_data<diagview<typename T1::elem_type>, T1, glue_schur>(X,Y.get_ref());
//   }



//! @}
