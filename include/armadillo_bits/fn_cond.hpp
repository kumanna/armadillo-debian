// Copyright (C) 2013 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_cond
//! @{


template<typename T1>
arma_warn_unused
inline
typename enable_if2<is_supported_blas_type<typename T1::elem_type>::value, typename T1::pod_type>::result
cond(const Base<typename T1::elem_type, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type T;
  
  Col<T> S;
  
  const bool status = auxlib::svd_dc(S, X);
  
  if(status == false)
    {
    arma_debug_warn("cond(): svd failed");
    
    return T(0);
    }
  
  if(S.n_elem > 0)
    {
    return T( max(S) / min(S) );
    }
  else
    {
    return T(0);
    }
  }




//! @}
