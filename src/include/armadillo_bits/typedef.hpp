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


//! \addtogroup typedef
//! @{


// //
// // only supported by C++0x, via #include <cstdint>, or by C99, via #include <stdint.h>
// 
// //! unsigned 8 bit type
// typedef uint8_t u8;
// 
// //! unsigned 16 bit type  
// typedef uint16_t u16;
// 
// //! unsigned 32 bit type
// typedef uint32_t u32;
//
// //! signed 32 bit type
// typedef int32_t s32;

#if UCHAR_MAX == 0xff
  //! unsigned 8 bit type
  typedef unsigned char u8;
  typedef          char s8;
#else
  #error "don't know how to typedef 'u8' on this system"
#endif

// NOTE:
// "signed char" is not the same as "char". 
// see http://www.embedded.com/columns/programmingpointers/206107018

#if USHRT_MAX == 0xffff
  //! unsigned 16 bit type  
  typedef unsigned short u16;
  typedef          short s16;
#else
  #error "don't know how to typedef 'u16' on this system"
#endif


#if   UINT_MAX  == 0xffffffff
  typedef unsigned int  u32;
  typedef          int  s32;
#elif ULONG_MAX == 0xffffffff
  typedef unsigned long u32;
  typedef          long s32;
#else
  #error "don't know how to typedef 'u32' on this system"
#endif


typedef std::complex<float>  cx_float;
typedef std::complex<double> cx_double;


typedef basic_mat   <unsigned char> uchar_mat;
typedef basic_colvec<unsigned char> uchar_vec;
typedef basic_colvec<unsigned char> uchar_colvec;
typedef basic_rowvec<unsigned char> uchar_rowvec;

typedef basic_mat   <u32> umat;
typedef basic_colvec<u32> uvec;
typedef basic_colvec<u32> ucolvec;
typedef basic_rowvec<u32> urowvec;

typedef basic_mat   <s32> imat;
typedef basic_colvec<s32> ivec;
typedef basic_colvec<s32> icolvec;
typedef basic_rowvec<s32> irowvec;

typedef basic_mat   <float> fmat;
typedef basic_colvec<float> fvec;
typedef basic_colvec<float> fcolvec;
typedef basic_rowvec<float> frowvec;

typedef basic_mat   <double> mat;
typedef basic_colvec<double> vec;
typedef basic_colvec<double> colvec;
typedef basic_rowvec<double> rowvec;

typedef basic_mat   <cx_float> cx_fmat;
typedef basic_colvec<cx_float> cx_fvec;
typedef basic_colvec<cx_float> cx_fcolvec;
typedef basic_rowvec<cx_float> cx_frowvec;

typedef basic_mat   <cx_double> cx_mat;
typedef basic_colvec<cx_double> cx_vec;
typedef basic_colvec<cx_double> cx_colvec;
typedef basic_rowvec<cx_double> cx_rowvec;


namespace junk
  {
  struct arma_elem_size_test
    {
  
    arma_static_assert<sizeof(u8) == 1> ERROR___TYPE_U8_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(s8) == 1> ERROR___TYPE_S8_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(u16) == 2> ERROR___TYPE_U16_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(s16) == 2> ERROR___TYPE_S16_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(u32) == 4> ERROR___TYPE_U32_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(s32) == 4> ERROR___TYPE_S32_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(float)  == 4> ERROR___TYPE_FLOAT_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(double) == 8> ERROR___TYPE_DOUBLE_HAS_UNSUPPORTED_SIZE;
    
    arma_static_assert<sizeof(std::complex<float>)  == 8>  ERROR___TYPE_COMPLEX_FLOAT_HAS_UNSUPPORTED_SIZE;
    arma_static_assert<sizeof(std::complex<double>) == 16> ERROR___TYPE_COMPLEX_DOUBLE_HAS_UNSUPPORTED_SIZE;
  
    };
  }

//! @}
