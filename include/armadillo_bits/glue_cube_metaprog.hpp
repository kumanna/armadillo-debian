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


//! \addtogroup glue_cube_metaprog
//! @{



//! \brief
//! Template metaprogram depth_lhs_cube
//! calculates the number of GlueCube<Tx,Ty, glue_type> instances on the left hand side argument of GlueCube<Tx,Ty, glue_cube_type>
//! i.e. it recursively expands each Tx, until the type of Tx is not "GlueCube<..,.., glue_cube_type>"  (i.e the "glue_cube_type" changes)

template<typename glue_cube_type, typename T1>
struct depth_lhs_cube
  {
  static const u32 num = 0;
  };

template<typename glue_cube_type, typename T1, typename T2>
struct depth_lhs_cube< glue_cube_type, GlueCube<T1,T2,glue_cube_type> >
  {
  static const u32 num = 1 + depth_lhs_cube<glue_cube_type, T1>::num;
  };



//! \brief
//! Template metaprogram cube_ptrs
//! fills a given array with addresses of cubes from a recursive instance of GlueCube<Tx,Ty, glue_cube_type>.
//! While parsing the recursive instance, if encountered objects are of type OpCube<..>,
//! they are converted to type 'Cube' first

template<typename glue_type, typename T1>
struct cube_ptrs
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 0;

  inline
  static
  void
  get_ptrs
    (
    const Cube<elem_type>** ptrs,
          bool*             del,
    const T1&               X
    )
    {
    ptrs[0] = 
      (
      is_Cube<T1>::value ?
        reinterpret_cast<const Cube<elem_type>*>(&X)
      :
        new Cube<elem_type>(X)
      );

    del[0] = 
      (
      is_Cube<T1>::value ?
        false
      :
        true
      );

    
    }
  
  };



template<typename glue_type, typename T1, typename T2>
struct cube_ptrs<glue_type, GlueCube<T1,T2,glue_type> >
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 1 + cube_ptrs<glue_type, T1>::num;
  
  inline
  static
  void
  get_ptrs
    (
    const Cube<elem_type>**          in_ptrs,
          bool*                      del,
    const GlueCube<T1,T2,glue_type>& X
    )
    {
    isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
    
    cube_ptrs<glue_type, T1>::get_ptrs(in_ptrs, del, X.A);
    
    in_ptrs[num]  = 
      (
      is_Cube<T2>::value ?
        reinterpret_cast<const Cube<elem_type>*>(&X.B)
      :
        new Cube<elem_type>(X.B)
      );
    
    del[num] = 
      (
      is_Cube<T2>::value ?
        false
      :
        true
      );
    }
  
  };



//! template metaprogram cube_ptrs_outcheck
//! builds on 'cube_ptrs' by also checking whether any of the input cubes are aliases of the output cube

template<typename glue_type, typename T1>
struct cube_ptrs_outcheck
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 0;

  inline
  static
  void
  get_ptrs
    (
    const Cube<elem_type>** ptrs,
          bool*             del,
    const T1&               X,
    const Cube<elem_type>*  out_ptr
    )
    {

    const bool same_ptr = 
      (
      is_Cube<T1>::value ?
        (
        (out_ptr == reinterpret_cast<const Cube<elem_type>*>(&X)) ?
          true
        :
          false
        )
      :
        false
      );

    
    ptrs[0] = 
      (
      same_ptr ?
        new Cube<elem_type>(X)
      :
        (
        is_Cube<T1>::value ?
          reinterpret_cast<const Cube<elem_type>*>(&X)
        :
          new Cube<elem_type>(X)
        )
      );

    
    del[0] = 
      (
      same_ptr ?
        true
      :
        (
        is_Cube<T1>::value ?
          false
        :
          true
        )
      );

    
    }
  
  };



template<typename glue_type, typename T1, typename T2>
struct cube_ptrs_outcheck<glue_type, GlueCube<T1,T2,glue_type> >
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 1 + cube_ptrs_outcheck<glue_type, T1>::num;
  
  inline
  static
  void
  get_ptrs
    (
    const Cube<elem_type>**          ptrs,
          bool*                      del,
    const GlueCube<T1,T2,glue_type>& X,
    const Cube<elem_type>*           out_ptr
    )
    {
    isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
    
    cube_ptrs_outcheck<glue_type, T1>::get_ptrs(ptrs, del, X.A, out_ptr);
    
    const bool same_ptr =
      (
      is_Cube<T2>::value ?
        (
        (out_ptr == reinterpret_cast<const Cube<elem_type>*>(&X.B)) ?
          true
        :
          false
        )
      :
        false
      );
    
    
    ptrs[num]  = 
      (
      same_ptr ?
        new Cube<elem_type>(X.B)
      :
        (
        is_Cube<T2>::value ?
          reinterpret_cast<const Cube<elem_type>*>(&X.B)
        :
          new Cube<elem_type>(X.B)
        )
      );
    
    
    del[num] = 
      (
      same_ptr ?
        true
      :
        (
        is_Cube<T2>::value ?
          false
        :
          true
        )
      );
    }

  };



//! @}
