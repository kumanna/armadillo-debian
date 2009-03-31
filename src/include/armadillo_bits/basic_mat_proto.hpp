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


//! \addtogroup basic_mat
//! @{

//! Dense matrix class

template<typename eT>
class basic_mat : public arma_base<eT, basic_mat<eT> >
  {
  public:
  
  typedef eT elem_type;  //!< the type of elements stored in the matrix
  
  typedef typename get_pod_type<elem_type>::pod_type pod_type;
  //!< if eT is std::complex, pod_type is the underlying type used by std::complex.
  //!< otherwise pod_type is the same as elem_type


  const u32 n_rows;     //!< number of rows in the matrix (read-only)
  const u32 n_cols;     //!< number of columns in the matrix (read-only)
  const u32 n_elem;     //!< number of elements in the matrix (read-only)
        u32 pad;        //!< padding for 64 bit alignment
  const eT* const mem;  //!< pointer to memory used by the matrix (memory is read-only)
  
  protected:
        eT  mem_local[ 16 ];
  
  
  public:
  
  inline     ~basic_mat();
  inline      basic_mat();
  
  inline      basic_mat(const u32 in_rows, const u32 in_cols);
  inline void  set_size(const u32 in_rows, const u32 in_cols);
  
  inline                  basic_mat(const char* text);
  inline const basic_mat& operator=(const char* text);
  inline                  basic_mat(const std::string& text);
  inline const basic_mat& operator=(const std::string& text);
  
  inline basic_mat(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols);
  
  inline const basic_mat&  operator=(const eT val);
  inline const basic_mat& operator+=(const eT val);
  inline const basic_mat& operator-=(const eT val);
  inline const basic_mat& operator*=(const eT val);
  inline const basic_mat& operator/=(const eT val);
  
  inline                   basic_mat(const basic_mat& m);
  inline const basic_mat&  operator=(const basic_mat& m);
  inline const basic_mat& operator+=(const basic_mat& m);
  inline const basic_mat& operator-=(const basic_mat& m);
  inline const basic_mat& operator*=(const basic_mat& m);
  inline const basic_mat& operator%=(const basic_mat& m);
  inline const basic_mat& operator/=(const basic_mat& m);

  template<typename T1, typename T2>
  inline explicit basic_mat(const arma_base<pod_type,T1>& A, const arma_base<pod_type,T2>& B);

  inline                   basic_mat(const subview<eT>& X);
  inline const basic_mat&  operator=(const subview<eT>& X);
  inline const basic_mat& operator+=(const subview<eT>& X);
  inline const basic_mat& operator-=(const subview<eT>& X);
  inline const basic_mat& operator*=(const subview<eT>& X);
  inline const basic_mat& operator%=(const subview<eT>& X);
  inline const basic_mat& operator/=(const subview<eT>& X);
  
  //inline explicit          basic_mat(const diagview<eT>& X);
  inline                   basic_mat(const diagview<eT>& X);
  inline const basic_mat&  operator=(const diagview<eT>& X);
  
  inline       subview_row<eT> row(const u32 row_num);
  inline const subview_row<eT> row(const u32 row_num) const;
  
  inline       subview_col<eT> col(const u32 col_num);
  inline const subview_col<eT> col(const u32 col_num) const;
  
  inline       subview<eT> rows(const u32 in_row1, const u32 in_row2);
  inline const subview<eT> rows(const u32 in_row1, const u32 in_row2) const;
  
  inline       subview<eT> cols(const u32 in_col1, const u32 in_col2);
  inline const subview<eT> cols(const u32 in_col1, const u32 in_col2) const;
  
  inline       subview<eT> submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2);
  inline const subview<eT> submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2) const;

  inline       diagview<eT> diag(const s32 in_id = 0);
  inline const diagview<eT> diag(const s32 in_id = 0) const;
    
  inline void swap_rows(const u32 in_row1, const u32 in_row2);
  inline void swap_cols(const u32 in_col1, const u32 in_col2);
  
  template<typename T1, typename op_type> inline                   basic_mat(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_mat&  operator=(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_mat& operator+=(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_mat& operator-=(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_mat& operator*=(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_mat& operator%=(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_mat& operator/=(const op_data<T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                   basic_mat(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_mat&  operator=(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_mat& operator+=(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_mat& operator-=(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_mat& operator*=(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_mat& operator%=(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_mat& operator/=(const glue_data<T1, T2, glue_type>& X);
  
  
  inline eT& operator[] (const u32 i);
  inline eT  operator[] (const u32 i) const;
  inline eT& operator() (const u32 i);
  inline eT  operator() (const u32 i) const;
  
  inline eT& at         (const u32 in_row, const u32 in_col);
  inline eT  at         (const u32 in_row, const u32 in_col) const;
  inline eT& operator() (const u32 in_row, const u32 in_col);
  inline eT  operator() (const u32 in_row, const u32 in_col) const;

  inline bool is_vec() const;
  inline bool is_square() const;

  inline       eT* colptr(const u32 in_col);
  inline const eT* colptr(const u32 in_col) const;
  
  inline       eT* memptr();
  inline const eT* memptr() const;

  inline void print(const std::string extra_text = "") const;

  inline void fill(const eT val);
  inline void zeros();
  inline void zeros(const u32 in_rows, const u32 in_cols);
  
  inline void reset();
  
  inline void save(const std::string name, const file_type type = arma_binary) const;
  inline void load(const std::string name, const file_type type = auto_detect);
  
  
  protected: 
  
  inline void init(const u32 in_rows, const u32 in_cols);
  inline void init(const std::string& text);
  inline void init(const basic_mat& x);
  };


//! @}
