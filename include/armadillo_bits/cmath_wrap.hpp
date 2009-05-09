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



//! \addtogroup cmath_wrap
//! @{


// TODO:
// after upgrade to cmake >= 2.6,
// change to use ARMA_SYS_ISFINITE
// instead of checking for GCC specific macros


template<typename eT>
arma_inline
bool
arma_isfinite(eT val)
  {
  return true;
  }



arma_inline
bool
arma_isfinite(float x)
  {
  #if defined(_GLIBCXX_USE_C99_MATH)
    {
    return (std::isfinite(x) != 0);
    }
  #else
    {
    const bool x_is_inf = ( (x == x) && ((x - x) != float(0)) );
    const bool x_is_nan = (x != x);

    return ( (x_is_inf == false) && (x_is_nan == false) );
    }
  #endif
  }



arma_inline
bool
arma_isfinite(double x)
  {
  #if defined(_GLIBCXX_USE_C99_MATH)
    {
    return (std::isfinite(x) != 0);
    }
  #else
    {
    const bool x_is_inf = ( (x == x) && ((x - x) != double(0)) );
    const bool x_is_nan = (x != x);

    return ( (x_is_inf == false) && (x_is_nan == false) );
    }
  #endif
  }



template<typename T>
arma_inline
bool
arma_isfinite(const std::complex<T>& x)
  {
  if( (arma_isfinite(x.real()) == false) || (arma_isfinite(x.imag()) == false) )
    {
    return false;
    }
  else
    {
    return true;
    }
  }



//! @}
