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


//! \addtogroup unwrap
//! @{


template<typename T1>
class unwrap
  {
  public:
  inline unwrap(const T1& A)
    {
    arma_type_check< is_arma_type<T1>::value == false >::apply();
    }
  };



template<>
template<typename eT>
class unwrap< basic_mat<eT> >
  {
  public:

  inline unwrap(const basic_mat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const basic_mat<eT>& M;
  
  };



template<>
template<typename eT>
class unwrap< basic_rowvec<eT> >
  {
  public:

  inline unwrap(const basic_rowvec<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const basic_rowvec<eT>& M;
  
  };



template<>
template<typename eT>
class unwrap< basic_colvec<eT> >
  {
  public:

  inline unwrap(const basic_colvec<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const basic_colvec<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap< op_data<T1, op_type> >
  {
  public:

  inline unwrap(const op_data<T1, op_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  typedef typename T1::elem_type elem_type;
  const basic_mat<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_type>
class unwrap< glue_data<T1, T2, glue_type> >
  {
  public:

  inline unwrap(const glue_data<T1, T2, glue_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  typedef typename T1::elem_type elem_type;
  const basic_mat<elem_type> M;
  
  };


template<>
template<typename eT>
class unwrap< subview<eT> >
  {
  public:

  inline unwrap(const subview<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const basic_mat<eT> M;
  
  };


template<>
template<typename eT>
class unwrap< diagview<eT> >
  {
  public:

  inline unwrap(const diagview<eT> & A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const basic_mat<eT> M;
  
  };



//
//
//


template<typename T1>
class unwrap_to_elem_access
  {
  public:
  inline unwrap_to_elem_access(const T1& A)
    {
    arma_type_check< is_arma_type<T1>::value == false >::apply();
    }
  };



template<>
template<typename eT>
class unwrap_to_elem_access< basic_mat<eT> >
  {
  public:

  inline unwrap_to_elem_access(const basic_mat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col); }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  
  const basic_mat<eT>& M;
  
  };



template<>
template<typename eT>
class unwrap_to_elem_access< basic_rowvec<eT> >
  {
  public:

  inline unwrap_to_elem_access(const basic_rowvec<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col); }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const basic_rowvec<eT>& M;
  
  };



template<>
template<typename eT>
class unwrap_to_elem_access< op_data< basic_rowvec<eT>, op_trans> >
  {
  public:

  // NOTE:
  // currently accessing M.n_rows and M.n_cols will give wrong results
  
  inline unwrap_to_elem_access(const op_data<basic_rowvec<eT>, op_trans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col); }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const basic_rowvec<eT>& M;
  
  };



template<>
template<typename eT>
class unwrap_to_elem_access< basic_colvec<eT> >
  {
  public:

  inline unwrap_to_elem_access(const basic_colvec<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col); }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const basic_colvec<eT>& M;
  
  };



template<>
template<typename eT>
class unwrap_to_elem_access< op_data<basic_colvec<eT>, op_trans> >
  {
  public:

  inline unwrap_to_elem_access(const op_data<basic_colvec<eT>, op_trans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  // NOTE: use of in_row and in_col is swapped (due to transpose operation)
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_col,in_row); }
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_col,in_row); }
  
  const basic_colvec<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap_to_elem_access< op_data<T1, op_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline unwrap_to_elem_access(const op_data<T1, op_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline elem_type  operator[](const u32 i) const { return M[i]; }
  inline elem_type  operator()(const u32 i) const { return M(i); }
  
  inline elem_type  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col); }
  inline elem_type          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const basic_mat<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_type>
class unwrap_to_elem_access< glue_data<T1, T2, glue_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline unwrap_to_elem_access(const glue_data<T1, T2, glue_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline elem_type  operator[](const u32 i) const { return M[i]; }
  inline elem_type  operator()(const u32 i) const { return M(i); }
  
  inline elem_type  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col); }
  inline elem_type          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
    
  const basic_mat<elem_type> M;
  
  };


template<>
template<typename eT>
class unwrap_to_elem_access< subview<eT> >
  {
  public:

  inline unwrap_to_elem_access(const subview<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT operator[](const u32 i) const { return M[i]; }
  inline eT operator()(const u32 i) const { return M(i); }
  
  inline eT operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col); }
  inline eT         at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const subview<eT>& M;
  
  };


//
//
//

template<typename T1>
class unwrap_check
  {
  private:
  template<typename eT> inline unwrap_check(const T1& A, const basic_mat<eT>&    B);
  template<typename eT> inline unwrap_check(const T1& A, const basic_rowvec<eT>& B);
  template<typename eT> inline unwrap_check(const T1& A, const basic_colvec<eT>& B);
  template<typename eT> inline unwrap_check(const T1& A, const subview<eT>&      B);
  template<typename eT> inline unwrap_check(const T1& A, const diagview<eT>&     B);
  };


template <>
template<typename eT>
class unwrap_check< basic_mat<eT> >
  {
  public:

  inline
  unwrap_check(const basic_mat<eT>& A, const basic_mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_mat<eT>*>(&B)) ? new basic_mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_mat<eT>*>(&B)) ? (*M_local)           : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const basic_mat<eT>& A, const basic_rowvec<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_mat<eT>*>(&B)) ? new basic_mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_mat<eT>*>(&B)) ? (*M_local)           : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const basic_mat<eT>& A, const basic_colvec<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_mat<eT>*>(&B)) ? new basic_mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_mat<eT>*>(&B)) ? (*M_local)           : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const basic_mat<eT>& A, const subview<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_mat<eT>*>(&B.m)) ? new basic_mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_mat<eT>*>(&B.m)) ? (*M_local)           : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const basic_mat<eT>& A, const diagview<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_mat<eT>*>(&B.m)) ? new basic_mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_mat<eT>*>(&B.m)) ? (*M_local)           : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  // the order below is important
  const basic_mat<eT>* M_local;
  const basic_mat<eT>& M;
  
  };



template <>
template<typename eT>
class unwrap_check< basic_rowvec<eT> >
  {
  public:

  inline
  unwrap_check(const basic_rowvec<eT>& A, const basic_mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B)) ? new basic_rowvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const basic_rowvec<eT>& A, const basic_rowvec<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B)) ? new basic_rowvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  unwrap_check(const basic_rowvec<eT>& A, const basic_colvec<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B)) ? new basic_rowvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const basic_rowvec<eT>& A, const subview<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B.m)) ? new basic_rowvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B.m)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const basic_rowvec<eT>& A, const diagview<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B.m)) ? new basic_rowvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_rowvec<eT>*>(&B.m)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }

  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      delete M_local;
    }
  
  
  // the order below is important
  const basic_rowvec<eT>* M_local;
  const basic_rowvec<eT>& M;
  
  };




template <>
template<typename eT>
class unwrap_check< basic_colvec<eT> >
  {
  public:

  inline
  unwrap_check(const basic_colvec<eT>& A, const basic_mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B)) ? new basic_colvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const basic_colvec<eT>& A, const basic_rowvec<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B)) ? new basic_colvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const basic_colvec<eT>& A, const basic_colvec<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B)) ? new basic_colvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const basic_colvec<eT>& A, const subview<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B.m)) ? new basic_colvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B.m)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const basic_colvec<eT>& A, const diagview<eT>& B)
    : M_local( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B.m)) ? new basic_colvec<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const basic_colvec<eT>*>(&B.m)) ? (*M_local)              : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      delete M_local;
    }
  
  
  // the order below is important
  const basic_colvec<eT>* M_local;
  const basic_colvec<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap_check< op_data<T1, op_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  //template<typename eT>
  inline
  unwrap_check(const op_data<T1,op_type>& A, const basic_mat<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  //template<typename eT>
  inline
  unwrap_check(const op_data<T1,op_type>& A, const basic_rowvec<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  //template<typename eT>
  inline
  unwrap_check(const op_data<T1,op_type>& A, const basic_colvec<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    }
  
  const basic_mat<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_type>
class unwrap_check< glue_data<T1, T2, glue_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline
  unwrap_check(const glue_data<T1, T2, glue_type>& A, const basic_mat<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const glue_data<T1, T2, glue_type>& A, const basic_rowvec<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const glue_data<T1, T2, glue_type>& A, const basic_colvec<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    }
  
  
  const basic_mat<elem_type> M;
  
  };




template<>
template<typename eT>
class unwrap_check< subview<eT> >
  {
  public:

  template<typename T2>
  inline unwrap_check(const subview<eT>& A, const T2& junk)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const basic_mat<eT> M;
  
  };


template<>
template<typename eT>
class unwrap_check< diagview<eT> >
  {
  public:

  template<typename T2>
  inline unwrap_check(const diagview<eT> & A, const T2& junk)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const basic_mat<eT> M;
  
  };


//! @}
