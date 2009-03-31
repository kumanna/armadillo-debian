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


//! \addtogroup subview
//! @{


//! Class for storing data required to construct or apply operations to a submatrix
//! (i.e. where the submatrix starts and ends as well as a reference/pointer to the original matrix),
template<typename eT>
class subview : public arma_base<eT, subview<eT> >
  {
  public:  const basic_mat<eT>& m;
  protected:     basic_mat<eT>* m_ptr;
  
  public:
  
  typedef eT elem_type;
  typedef typename get_pod_type<elem_type>::pod_type pod_type;

  const u32 aux_row1;
  const u32 aux_col1;
  
  const u32 aux_row2;
  const u32 aux_col2;
  
  const u32 n_rows;
  const u32 n_cols;
  const u32 n_elem;
  
  
  protected:
  
  inline subview(const basic_mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2);
  inline subview(      basic_mat<eT>& in_m, const u32 in_row1, const u32 in_col1, const u32 in_row2,  const u32 in_col2);
  
  
  public:
  
  inline ~subview();
  
  inline void operator+= (const eT val);
  inline void operator-= (const eT val);
  inline void operator*= (const eT val);
  inline void operator/= (const eT val);
  
  // deliberately returning void
  template<typename T1> inline void operator=  (const arma_base<eT,T1>& x);
  template<typename T1> inline void operator+= (const arma_base<eT,T1>& x);
  template<typename T1> inline void operator-= (const arma_base<eT,T1>& x);
  template<typename T1> inline void operator%= (const arma_base<eT,T1>& x);
  template<typename T1> inline void operator/= (const arma_base<eT,T1>& x);
  
  inline void operator=  (const subview& x);
  inline void operator+= (const subview& x);
  inline void operator-= (const subview& x);
  inline void operator%= (const subview& x);
  inline void operator/= (const subview& x);
  
  inline static void extract(basic_mat<eT>& out, const subview& in);
  
  inline static void  plus_inplace(basic_mat<eT>& out, const subview& in);
  inline static void times_inplace(basic_mat<eT>& out, const subview& in);
  inline static void minus_inplace(basic_mat<eT>& out, const subview& in);
  inline static void schur_inplace(basic_mat<eT>& out, const subview& in);
  inline static void   div_inplace(basic_mat<eT>& out, const subview& in);
  
  inline void fill(const eT val);
  inline void zeros();
  
  inline eT& operator[](const u32 i);
  inline eT  operator[](const u32 i) const;
  
  inline eT& operator()(const u32 i);
  inline eT  operator()(const u32 i) const;
  
  inline eT& operator()(const u32 in_row, const u32 in_col);
  inline eT  operator()(const u32 in_row, const u32 in_col) const;
  
  inline eT&         at(const u32 in_row, const u32 in_col);
  inline eT          at(const u32 in_row, const u32 in_col) const;
  
  inline       eT* colptr(const u32 in_col);
  inline const eT* colptr(const u32 in_col) const;
  
  inline bool check_overlap(const subview& x) const;
  
  
  private:
  
  friend class basic_mat<eT>;
  subview();
  };



template<typename eT>
class subview_col : public subview<eT>
  {
  public:
  
  typedef eT elem_type;
  typedef typename get_pod_type<elem_type>::pod_type pod_type;
  
  inline void operator= (const subview<eT>&   x);
  inline void operator= (const subview_col&   x);
  
  template<typename T1>
  inline void operator= (const arma_base<eT,T1>& x);
  
  
  protected:
  
  inline subview_col(const basic_mat<eT>& in_m, const u32 in_col);
  inline subview_col(      basic_mat<eT>& in_m, const u32 in_col);
  
  
  private:
  
  friend class basic_mat<eT>;
  subview_col();
  };
  


template<typename eT>
class subview_row : public subview<eT>
  {
  public:
  
  typedef eT elem_type;
  typedef typename get_pod_type<elem_type>::pod_type pod_type;
  
  inline void operator= (const subview<eT>&   x);
  inline void operator= (const subview_row&   x);
  
  template<typename T1>
  inline void operator= (const arma_base<eT,T1>& x);
  
  
  protected:
  
  inline subview_row(const basic_mat<eT>& in_m, const u32 in_row);
  inline subview_row(      basic_mat<eT>& in_m, const u32 in_row);
  
  
  private:
  
  friend class basic_mat<eT>;
  subview_row();
  };



//! @}
