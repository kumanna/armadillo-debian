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


//! \addtogroup arma_base
//! @{


//! Class for static polymorphism, modelled after the "Curiously Recurring Template Pattern" (CRTP).
//! Used for type-safe downcasting in functions that restrict their input(s) to be classes that are
//! derived from arma_base (e.g. basic_mat, op_data, glue_data, diagview, subview).
//! An arma_base object can be converted to a basic_mat object by the unwrap class.

template<typename elem_type, typename arma_derived>
struct arma_base
  {
  
  inline
  const arma_derived&
  get_ref() const
    {
    return static_cast<const arma_derived&>(*this);
    }

  };



template<typename elem_type, typename arma_derived_vec>
struct arma_base_vec
  {
  
  inline
  const arma_derived_vec&
  get_ref() const
    {
    return static_cast<const arma_derived_vec&>(*this);
    }

  };



//! @}
