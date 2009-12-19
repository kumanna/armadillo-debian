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


//! \addtogroup unwrap_cube
//! @{



template<typename T1>
class unwrap_cube
  {
  public:
  inline unwrap_cube(const T1& A)
    {
    arma_type_check< is_arma_cube_type<T1>::value == false >::apply();
    }
  };



template<typename eT>
class unwrap_cube< Cube<eT> >
  {
  public:

  inline unwrap_cube(const Cube<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Cube<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap_cube< OpCube<T1, op_type> >
  {
  public:

  inline unwrap_cube(const OpCube<T1, op_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  typedef typename T1::elem_type elem_type;
  const Cube<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_cube_type>
class unwrap_cube< GlueCube<T1, T2, glue_cube_type> >
  {
  public:

  inline unwrap_cube(const GlueCube<T1, T2, glue_cube_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  typedef typename T1::elem_type elem_type;
  const Cube<elem_type> M;
  
  };



template<typename eT>
class unwrap_cube< subview_cube<eT> >
  {
  public:

  inline unwrap_cube(const subview_cube<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Cube<eT> M;
  
  };



//
//
//



template<typename T1>
class unwrap_cube_check
  {
  private:
  template<typename eT> inline unwrap_cube_check(const T1& A, const Cube<eT>& B);
  template<typename eT> inline unwrap_cube_check(const T1& A, const subview_cube<eT>&  B);
  };



template<typename eT>
class unwrap_cube_check< Cube<eT> >
  {
  public:

  inline
  unwrap_cube_check(const Cube<eT>& A, const Cube<eT>& B)
    : M_local( (&A == reinterpret_cast<const Cube<eT>*>(&B)) ? new Cube<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Cube<eT>*>(&B)) ? (*M_local)      : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_cube_check(const Cube<eT>& A, const subview_cube<eT>& B)
    : M_local( (&A == reinterpret_cast<const Cube<eT>*>(&B.m)) ? new Cube<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Cube<eT>*>(&B.m)) ? (*M_local)      : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  ~unwrap_cube_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  // the order below is important
  const Cube<eT>* M_local;
  const Cube<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap_cube_check< OpCube<T1, op_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline
  unwrap_cube_check(const OpCube<T1,op_type>& A, const Cube<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_cube_check()
    {
    arma_extra_debug_sigprint();
    }
  
  const Cube<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_cube_type>
class unwrap_cube_check< GlueCube<T1, T2, glue_cube_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline
  unwrap_cube_check(const GlueCube<T1, T2, glue_cube_type>& A, const Cube<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_cube_check()
    {
    arma_extra_debug_sigprint();
    }
  
  
  const Cube<elem_type> M;
  
  };



template<typename eT>
class unwrap_cube_check< subview_cube<eT> >
  {
  public:

  template<typename T2>
  inline unwrap_cube_check(const subview_cube<eT>& A, const T2& junk)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Cube<eT> M;
  
  };



//! if the given object in not a cube, unwrap it into the given 'out' cube (i.e. do not create another cube)
//! and provide a reference to the 'out' cube.
//! if the given object is a cube, set the size of the 'out' cube to be the same as the given object
//! and provide a reference to the given object.

template<typename T1>
class unwrap_cube_write
  {
  private:
  template<typename eT> inline unwrap_cube_write(Cube<eT>& out, const T1& in);
  };


//template <>
template<typename eT>
class unwrap_cube_write< Cube<eT> >
  {
  public:
  
  inline
  unwrap_cube_write(Cube<eT>& out, const Cube<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    out.copy_size(A);
    }
  
  
  inline
  ~unwrap_cube_write()
    {
    arma_extra_debug_sigprint();
    }
  
  
  const Cube<eT>& M;
  };



template<typename T1, typename op_type>
class unwrap_cube_write< OpCube<T1, op_type> >
  {
  public:
  typedef typename T1::elem_type eT;

  //template<typename eT>
  inline
  unwrap_cube_write(Cube<eT>& out, const OpCube<T1,op_type>& A)
    : M(out)
    {
    arma_extra_debug_sigprint();
    out = A;
    }
  
  
  inline
  ~unwrap_cube_write()
    {
    arma_extra_debug_sigprint();
    }
  
  const Cube<eT>& M;
  };



template<typename T1, typename T2, typename glue_type>
class unwrap_cube_write< GlueCube<T1, T2, glue_type> >
  {
  public:
  typedef typename T1::elem_type eT;
  
  inline
  unwrap_cube_write(Cube<eT>& out, const GlueCube<T1, T2, glue_type>& A)
    : M(out)
    {
    arma_extra_debug_sigprint();
    out = A;
    }
  
  
  inline
  ~unwrap_cube_write()
    {
    arma_extra_debug_sigprint();
    }
  
  
  const Cube<eT>& M;
  };




//template<>
template<typename eT>
class unwrap_cube_write< subview_cube<eT> >
  {
  public:
  
  inline
  unwrap_cube_write(Cube<eT>& out, const subview_cube<eT>& A)
    : M(out)
    {
    arma_extra_debug_sigprint();
    out = A;
    }
  
  
  inline
  ~unwrap_cube_write()
    {
    arma_extra_debug_sigprint();
    }
  
  
  const Cube<eT>& M;
  };



//! @}
