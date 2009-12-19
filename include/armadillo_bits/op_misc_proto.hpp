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



//! \addtogroup op_misc
//! @{


class op_log
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_log>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_log>& in);
  };
  


class op_trunc_log
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_trunc_log>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_trunc_log>& in);
  };
  


class op_log10
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_log10>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_log10>& in);
  };
  


class op_exp
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_exp>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_exp>& in);
  };
  


class op_trunc_exp
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_trunc_exp>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_trunc_exp>& in);
  };
  


class op_sqrt
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_sqrt>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_sqrt>& in);
  };
  


class op_square
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_square>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_square>& in);
  };



class op_pow
  {
  public:
  
  template<typename T1, typename T2> static inline T1 internal_pow(const T1 base, const T2 exponent);
  
  template<typename T2> static inline          char internal_pow(const          char base, const T2 exponent);
  template<typename T2> static inline unsigned char internal_pow(const unsigned char base, const T2 exponent);
  
  template<typename T2> static inline          int  internal_pow(const          int  base, const T2 exponent);
  template<typename T2> static inline unsigned int  internal_pow(const unsigned int  base, const T2 exponent);
  
  
  template<typename T1> inline static void apply( Mat<typename T1::pod_type>& out, const     Op<T1,op_pow>& in);
  template<typename T1> inline static void apply(Cube<typename T1::pod_type>& out, const OpCube<T1,op_pow>& in);
  
  template<typename T1> inline static void apply( Mat< std::complex<typename T1::pod_type> >& out, const     Op<T1,op_pow>& in);
  template<typename T1> inline static void apply(Cube< std::complex<typename T1::pod_type> >& out, const OpCube<T1,op_pow>& in);
  };



class op_pow_s32
  {
  public:
  
  template<typename T1> static inline T1 internal_pow(const T1 base, const int exponent);
  
  static inline          char internal_pow(const          char base, const int exponent);
  static inline unsigned char internal_pow(const unsigned char base, const int exponent);
  
  static inline          int  internal_pow(const          int  base, const int exponent);
  static inline unsigned int  internal_pow(const unsigned int  base, const int exponent);
  

  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_pow_s32>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_pow_s32>& in);
  };



class op_conj
  {
  public:
  
  template<typename T1> inline static void apply( Mat<typename T1::elem_type>& out, const     Op<T1,op_conj>& in);
  template<typename T1> inline static void apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_conj>& in);
  };


//! @}
