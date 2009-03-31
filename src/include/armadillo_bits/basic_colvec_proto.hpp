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


//! \addtogroup basic_colvec
//! @{

//! Class for column vectors (matrices with only column)

template<typename eT>
class basic_colvec : public basic_mat<eT>, public arma_base_vec<eT, basic_colvec<eT> >
  {
  public:
  
  typedef eT elem_type;
  typedef typename get_pod_type<elem_type>::pod_type pod_type;
  
  
  inline                     basic_colvec();
  inline explicit            basic_colvec(const u32 n_elem);
  
  inline                     basic_colvec(const char* text);
  inline const basic_colvec&    operator=(const char* text);  // TODO: std::string input
  
  inline                     basic_colvec(const basic_colvec& X);
  inline const basic_colvec&    operator=(const basic_colvec& X);
  
  //inline explicit            basic_colvec(const basic_mat<eT>& X);
  inline                     basic_colvec(const basic_mat<eT>& X);
  inline const basic_colvec&    operator=(const basic_mat<eT>& X);
  inline const basic_colvec&   operator*=(const basic_mat<eT>& X);
  
  inline                     basic_colvec(const eT* aux_mem, const u32 aux_length);
  
  template<typename T1, typename T2>
  inline explicit basic_colvec(const arma_base<pod_type,T1>& A, const arma_base<pod_type,T2>& B);

  inline                     basic_colvec(const subview<eT>& X);
  inline const basic_colvec&    operator=(const subview<eT>& X);
  inline const basic_colvec&   operator*=(const subview<eT>& X);
  
  inline                     basic_colvec(const diagview<eT>& X);
  inline const basic_colvec&    operator=(const diagview<eT>& X);
  inline const basic_colvec&   operator*=(const diagview<eT>& X);
  
  template<typename T1, typename op_type> inline                     basic_colvec(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_colvec&    operator=(const op_data<T1, op_type>& X);
  template<typename T1, typename op_type> inline const basic_colvec&   operator*=(const op_data<T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type> inline                     basic_colvec(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_colvec&    operator=(const glue_data<T1, T2, glue_type>& X);
  template<typename T1, typename T2, typename glue_type> inline const basic_colvec&   operator*=(const glue_data<T1, T2, glue_type>& X);
    
  inline void set_size(const u32 n_elem);
  inline void set_size(const u32 n_rows, const u32 n_cols);
  
  inline void zeros();
  inline void zeros(const u32 n_elem);
  inline void zeros(const u32 n_rows, const u32 n_cols);


  inline void load(const std::string name, const file_type type = auto_detect);
  };


//! @}
