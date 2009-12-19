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


//! \addtogroup OpCube
//! @{


//! Analog of the Op class, intended for cubes

template<typename T1, typename op_type>
class OpCube : public BaseCube<typename T1::elem_type, OpCube<T1, op_type> >
  {
  public:
  
  typedef typename T1::elem_type elem_type;
  typedef typename get_pod_type<elem_type>::pod_type pod_type;
  
  
  inline explicit OpCube(const BaseCube<typename T1::elem_type, T1>& in_m);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline          OpCube(const BaseCube<typename T1::elem_type, T1>& in_m, const elem_type in_aux, const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline          OpCube(const u32 in_aux_u32_a, const u32 in_aux_u32_b, const u32 in_aux_u32_c);
  inline         ~OpCube();
  
  // TODO:MAT: if the restriction to inputs of type BaseCube works, the use of the "junk" parameter may no longer be necessary

  const T1&       m;          //!< storage of reference to the operand (e.g. a cube)
  const elem_type aux;        //!< storage of auxiliary data, user defined format
  const u32       aux_u32_a;  //!< storage of auxiliary data, u32 format
  const u32       aux_u32_b;  //!< storage of auxiliary data, u32 format
  const u32       aux_u32_c;  //!< storage of auxiliary data, u32 format
  
  };



//! @}
