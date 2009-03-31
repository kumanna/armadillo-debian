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


/* #undef ARMA_USE_ATLAS */
/* #undef ARMA_USE_LAPACK */
/* #undef ARMA_USE_BLAS */
/* #undef ARMA_USE_BOOST */
/* #undef ARMA_EXTRA_DEBUG */
/* #undef ARMA_NO_DEBUG */

#if defined(ARMA_USE_ATLAS)
  #if !defined(ARMA_ATLAS_INCLUDE_DIR)
    #define ARMA_ATLAS_INCLUDE_DIR
  #endif
#endif

