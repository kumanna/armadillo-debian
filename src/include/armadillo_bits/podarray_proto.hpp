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


//! \addtogroup podarray
//! @{



//! A lightweight array for POD types. If the amount of memory requested is small, the stack is used.

template<typename T1>
class podarray
  {
  public:
  
  const u32       n_elem;    //!< number of elements held
        u32       pad;
  
  const T1* const mem;       //!< pointer to memory used by the object
  
  
  protected:
  
  T1 mem_local[ 128 / sizeof(T1) ];
  //!< Internal memory, allocated on the stack.
  //!< Designed to avoid calling the 'new' operator for small amounts of memory.
  
  
  public:
  
  inline ~podarray();
  inline podarray();
  
  inline podarray(const podarray& x);
  inline const podarray& operator=(const podarray& x);
  
  inline podarray(const u32 new_N);

  inline T1& operator[] (const u32 i);
  inline T1  operator[] (const u32 i) const;
  
  inline T1& operator() (const u32 i);
  inline T1  operator() (const u32 i) const;

  void set_size(const u32 new_n_elem);
  void fill(const T1 val);

  void zeros();
  void zeros(const u32 new_n_elem);

  inline T1* memptr();
  
  
  protected:
  
  inline void init(const u32 new_n_elem);

  };

//! @}
