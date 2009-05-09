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


//! \addtogroup constants
//! @{


// the long lengths of the constants are for future support of "long double"
// and any smart compiler that does high-precision computation at compile-time

template<typename eT>
class Math
  {
  public:
  
  //! ratio of any circle's circumference to its diameter
  static const eT pi()     { return eT(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679); }
  
  //! base of the natural logarithm
  static const eT e()      { return eT(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274); }
  
  //! Euler's constant, aka Euler-Mascheroni constant
  static const eT euler()  { return eT(0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495); }
  
  //! golden ratio
  static const eT gratio() { return eT(1.6180339887498948482045868343656381177203091798057628621354486227052604628189024497072072041893911374); }
  
  //! square root of 2
  static const eT sqrt2()  { return eT(1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727); }
  
  //! the difference between 1 and the least value greater than 1 that is representable
  static const eT eps()     { return std::numeric_limits<eT>::epsilon(); }
  
  //! log of the minimum representable value
  static const eT log_min() { static const eT out = std::log(std::numeric_limits<eT>::min()); return out; }
    
  //! log of the maximum representable value
  static const eT log_max() { static const eT out = std::log(std::numeric_limits<eT>::max()); return out; }
  };



typedef Math<float>  fmath;
typedef Math<double> math;



struct arma_version
  {
  static const unsigned int major = 0;
  static const unsigned int minor = 6;
  static const unsigned int patch = 10;
  };



struct arma_config
  {
  #if defined(ARMA_USE_ATLAS)
    static const bool atlas = true;
  #else
    static const bool atlas = false;
  #endif
  
  
  #if defined(ARMA_USE_LAPACK)
    static const bool lapack = true;
  #else
    static const bool lapack = false;
  #endif
  
  
  #if defined(ARMA_USE_BLAS)
    static const bool blas = true;
  #else
    static const bool blas = false;
  #endif


  #if defined(ARMA_USE_BOOST)
    static const bool boost = true;
  #else
    static const bool boost = false;
  #endif
  
  
  #if !defined(ARMA_NO_DEBUG) && !defined(NDEBUG)
    static const bool debug = true;
  #else
    static const bool debug = false;
  #endif
  
  
  #if defined(ARMA_EXTRA_DEBUG)
    static const bool extra_debug = true;
  #else
    static const bool extra_debug = false;
  #endif
  };

//! @}


