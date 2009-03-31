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


//! \addtogroup fn_conv_to
//! @{


//
//
// conversions between various mat types

template<typename out_eT, typename in_eT>
inline
void
copy_complex_elem(out_eT& out, const in_eT& in)
  {
  //arma_extra_debug_sigprint();
  
  out = out_eT(in);
  }



template<typename out_eT, typename in_T>
inline
void
copy_complex_elem(out_eT& out, const std::complex<in_T>& in)
  {
  //arma_extra_debug_sigprint();
  
  out = out_eT( in.real() );
  }



template<typename out_T, typename in_T>
inline
void
copy_complex_elem(std::complex<out_T>& out, const std::complex<in_T>& in)
  {
  //arma_extra_debug_sigprint();
  
  typedef std::complex<out_T> out_eT;
  out = out_eT(in);
  }



//
// scalar family


template<typename out_eT>
class conv_to
  {
  public:
  
  inline static out_eT from(const basic_mat< out_eT >& in);

  template<typename in_eT>
  inline static out_eT from(const basic_mat< in_eT >& in);

  template<typename in_eT, typename T1>
  inline static out_eT from(const arma_base<in_eT,T1>& in);
  
  };


template<typename out_eT>
inline
out_eT
conv_to<out_eT>::from(const basic_mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_elem != 1), "conv_to<>: matrix doesn't have exactly one element" );
  
  return in.mem[0];
  }



template<typename out_eT>
template<typename in_eT>
inline
out_eT
conv_to<out_eT>::from(const basic_mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< is_supported_elem_type<out_eT>::value == false >::apply();
  
  arma_debug_check( (in.n_elem != 1), "conv_to<>: matrix doesn't have exactly one element" );
  
  out_eT out;
  copy_complex_elem(out, in.mem[0]);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
out_eT
conv_to<out_eT>::from(const arma_base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< is_supported_elem_type<out_eT>::value == false >::apply();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to<out_eT>::from( tmp.M );
  }



//
//
// basic_mat family

template<typename out_eT>
class conv_to< basic_mat<out_eT> >
  {
  public:
  
  template<typename in_eT>
  inline static basic_mat<out_eT> from(const basic_mat< in_eT >& in);
  
  template<typename in_T>
  inline static basic_mat<out_eT> from(const basic_mat< std::complex<in_T> >& in);
  
  template<typename in_eT, typename T1>
  inline static basic_mat<out_eT> from(const arma_base<in_eT,T1>& in);
  
  
  
//   template<typename in_eT>
//   inline static basic_mat<out_eT> from(const std::vector< in_eT >& in);
//   
//   template<typename in_eT>
//   inline static basic_mat<out_eT> from(const std::vector< std::complex<in_eT> >& in);
  
  
  template<typename in_eT>
  inline static basic_mat<out_eT> from(const itpp::Mat< in_eT >& in);
  
  template<typename in_T>
  inline static basic_mat<out_eT> from(const itpp::Mat< std::complex<in_T> >& in);
  };



template<typename out_eT>
template<typename in_eT>
inline
basic_mat<out_eT>
conv_to< basic_mat<out_eT> >::from(const basic_mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
basic_mat<out_eT>
conv_to< basic_mat<out_eT> >::from(const basic_mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  basic_mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
basic_mat<out_eT>
conv_to< basic_mat<out_eT> >::from(const arma_base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to< basic_mat<out_eT> >::from( tmp.M );
  }



template<typename out_eT>
template<typename in_eT>
inline
basic_mat<out_eT>
conv_to< basic_mat<out_eT> >::from(const itpp::Mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  basic_mat<out_eT> out(in.rows(), in.cols());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
basic_mat<out_eT>
conv_to< basic_mat<out_eT> >::from(const itpp::Mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  basic_mat<out_eT> out(in.rows(), in.cols());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



//
//
// basic_rowvec family

template<typename out_eT>
class conv_to< basic_rowvec<out_eT> >
  {
  public:
  
  inline static basic_rowvec<out_eT> from(const basic_mat< out_eT >& in);
  
  template<typename in_eT>
  inline static basic_rowvec<out_eT> from(const basic_mat< in_eT >& in);
  
  template<typename in_T>  
  inline static basic_rowvec<out_eT> from(const basic_mat< std::complex<in_T> >& in);
  
  template<typename in_eT, typename T1>
  inline static basic_rowvec<out_eT> from(const arma_base<in_eT,T1>& in);
  
  
  
  template<typename in_eT>
  inline static basic_rowvec<out_eT> from(const itpp::Vec< in_eT >& in);
  
  template<typename in_T>  
  inline static basic_rowvec<out_eT> from(const itpp::Vec< std::complex<in_T> >& in);
  
  //inline static basic_rowvec<out_eT> from(const basic_colvec< out_eT >& in);
  //template<typename in_eT> inline static basic_rowvec<out_eT> from(const basic_colvec< in_eT >& in);
  //template<typename in_T>  inline static basic_rowvec<out_eT> from(const basic_colvec< std::complex<in_T> >& in);
  
  };



template<typename out_eT>
inline
basic_rowvec<out_eT>
conv_to< basic_rowvec<out_eT> >::from(const basic_mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_rows > 1), "conv_to<>: given matrix has more than one row");
  basic_rowvec<out_eT> out(in.n_cols);
  
  syslib::copy_elem(out.memptr(), in.mem, out.n_elem);
  
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
basic_rowvec<out_eT>
conv_to< basic_rowvec<out_eT> >::from(const basic_mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_rows > 1), "conv_to<>: given matrix has more than one row");
  basic_rowvec<out_eT> out(in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
basic_rowvec<out_eT>
conv_to< basic_rowvec<out_eT> >::from(const basic_mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  arma_debug_check( (in.n_rows > 1), "conv_to<>: given matrix has more than one row");
  basic_rowvec<out_eT> out(in.n_cols);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
basic_rowvec<out_eT>
conv_to< basic_rowvec<out_eT> >::from(const arma_base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();

  const unwrap<T1> tmp(in.get_ref());
  
  return conv_to< basic_rowvec<out_eT> >::from( tmp.M );
  }



template<typename out_eT>
template<typename in_eT>
inline
basic_rowvec<out_eT>
conv_to< basic_rowvec<out_eT> >::from(const itpp::Vec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  basic_rowvec<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
basic_rowvec<out_eT>
conv_to< basic_rowvec<out_eT> >::from(const itpp::Vec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  basic_rowvec<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



// template<typename out_eT>
// inline
// basic_rowvec<out_eT>
// conv_to< basic_rowvec<out_eT> >::from(const basic_colvec<out_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   return trans(in);
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_eT>
// inline
// basic_rowvec<out_eT>
// conv_to< basic_rowvec<out_eT> >::from(const basic_colvec<in_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   basic_rowvec<out_eT> out(in.n_rows);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     out_mem[i] = out_eT( in_mem[i] );
//     }
//   
//   return out;
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_T>
// inline
// basic_rowvec<out_eT>
// conv_to< basic_rowvec<out_eT> >::from(const basic_colvec< std::complex<in_T> >& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   typedef std::complex<in_T> in_eT; 
//   
//   basic_rowvec<out_eT> out(in.n_rows);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     copy_complex_elem(out_mem[i], in_mem[i]);
//     }
//   
//   return out;
//   }



//
//
// basic_colvec family

template<typename out_eT>
class conv_to< basic_colvec<out_eT> >
  {
  public:
  
  inline static basic_colvec<out_eT> from(const basic_mat< out_eT >& in);
  
  template<typename in_eT>
  inline static basic_colvec<out_eT> from(const basic_mat< in_eT >& in);
  
  template<typename in_T>
  inline static basic_colvec<out_eT> from(const basic_mat< std::complex<in_T> >& in);
  
  template<typename in_eT, typename T1>
  inline static basic_colvec<out_eT> from(const arma_base<in_eT,T1>& in);
  
  
  
  template<typename in_eT>
  inline static basic_colvec<out_eT> from(const itpp::Vec< in_eT >& in);
  
  template<typename in_T>
  inline static basic_colvec<out_eT> from(const itpp::Vec< std::complex<in_T> >& in);

//   inline static basic_colvec<out_eT> from(const basic_rowvec< out_eT >& in);
//   template<typename in_eT> inline static basic_colvec<out_eT> from(const basic_rowvec< in_eT >& in);
//   template<typename in_T>  inline static basic_colvec<out_eT> from(const basic_rowvec< std::complex<in_T> >& in);
  
  };



template<typename out_eT>
inline
basic_colvec<out_eT>
conv_to< basic_colvec<out_eT> >::from(const basic_mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_cols > 1), "conv_to<>: given matrix has more than one column");
  basic_colvec<out_eT> out(in.n_rows);
  
  syslib::copy_elem(out.memptr(), in.mem, out.n_elem);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
basic_colvec<out_eT>
conv_to< basic_colvec<out_eT> >::from(const basic_mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (in.n_cols > 1), "conv_to<>: given matrix has more than one column");
  basic_colvec<out_eT> out(in.n_rows);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
basic_colvec<out_eT>
conv_to< basic_colvec<out_eT> >::from(const basic_mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  arma_debug_check( (in.n_cols > 1), "conv_to<>: given matrix has more than one column");
  basic_colvec<out_eT> out(in.n_rows);
  
  const in_eT*  in_mem = in.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
basic_colvec<out_eT>
conv_to< basic_colvec<out_eT> >::from(const arma_base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());
  
  return conv_to< basic_colvec<out_eT> >::from( tmp.M );
  }



template<typename out_eT>
template<typename in_eT>
inline
basic_colvec<out_eT>
conv_to< basic_colvec<out_eT> >::from(const itpp::Vec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  basic_colvec<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
basic_colvec<out_eT>
conv_to< basic_colvec<out_eT> >::from(const itpp::Vec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  basic_colvec<out_eT> out(in.length());
  
  const in_eT*  in_mem = in._data();
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



// template<typename out_eT>
// inline
// basic_colvec<out_eT>
// conv_to< basic_colvec<out_eT> >::from(const basic_rowvec<out_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   return trans(in);
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_eT>
// inline
// basic_colvec<out_eT>
// conv_to< basic_colvec<out_eT> >::from(const basic_rowvec<in_eT>& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   basic_colvec<out_eT> out(in.n_cols);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     out_mem[i] = out_eT( in_mem[i] );
//     }
//   
//   return out;
//   }
// 
// 
// 
// template<typename out_eT>
// template<typename in_T>
// inline
// basic_colvec<out_eT>
// conv_to< basic_colvec<out_eT> >::from(const basic_rowvec< std::complex<in_T> >& in)
//   {
//   arma_extra_debug_sigprint();
//   
//   typedef std::complex<in_T> in_eT; 
//   
//   basic_colvec<out_eT> out(in.n_cols);
//   
//   const in_eT*  in_mem = in.mem;
//   out_eT*      out_mem = out.memptr();
//   
//   for(u32 i=0; i<out.n_elem; ++i)
//     {
//     copy_complex_elem(out_mem[i], in_mem[i]);
//     }
//   
//   return out;
//   }


//
//
// itpp::Mat family

template<typename out_eT>
class conv_to< itpp::Mat<out_eT> >
  {
  public:
  
  inline static itpp::Mat<out_eT> from(const basic_mat< out_eT >& in);
  
  inline static itpp::Mat<out_eT> from(const basic_colvec< out_eT >& in);
  
  inline static itpp::Mat<out_eT> from(const basic_rowvec< out_eT >& in);

  
  
  template<typename in_eT>
  inline static itpp::Mat<out_eT> from(const basic_mat< in_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Mat<out_eT> from(const basic_colvec< in_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Mat<out_eT> from(const basic_rowvec< in_eT >& in);
  
  
  template<typename in_T>
  inline static itpp::Mat<out_eT> from(const basic_mat< std::complex<in_T> >& in);
  
  template<typename in_T>
  inline static itpp::Mat<out_eT> from(const basic_colvec< std::complex<in_T> >& in);
  
  template<typename in_T>
  inline static itpp::Mat<out_eT> from(const basic_rowvec< std::complex<in_T> >& in);
  
  
  
  template<typename in_eT, typename T1>
  inline static itpp::Mat<out_eT> from(const arma_base<in_eT,T1>& in);
  
  };



template<typename out_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Mat<out_eT> out(in.n_rows, in.n_cols);
  
  syslib::copy_elem(out._data(), in.mem, in.n_elem);
  
  return out;
  }



template<typename out_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_colvec<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const basic_mat<out_eT>& >(in) );
  }



template<typename out_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_rowvec<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const basic_mat<out_eT>& >(in) );
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT* in_mem = in.memptr();
  out_eT*     out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_colvec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const basic_mat<in_eT>& >(in) );
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_rowvec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const basic_mat<in_eT>& >(in) );
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  itpp::Mat<out_eT> out(in.n_rows, in.n_cols);
  
  const in_eT* in_mem = in.memptr();
  out_eT*     out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_colvec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const basic_mat< std::complex<in_T> >& >(in) );
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const basic_rowvec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  return conv_to< itpp::Mat<out_eT> >::from( reinterpret_cast<const basic_mat< std::complex<in_T> >& >(in) );
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
itpp::Mat<out_eT>
conv_to< itpp::Mat<out_eT> >::from(const arma_base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to< itpp::Mat<out_eT> >::from( tmp.M );
  }



//
//
// itpp::Vec family

template<typename out_eT>
class conv_to< itpp::Vec<out_eT> >
  {
  public:
  
  inline static itpp::Vec<out_eT> from(const basic_mat< out_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Vec<out_eT> from(const basic_mat< in_eT >& in);
  
  template<typename in_T>  
  inline static itpp::Vec<out_eT> from(const basic_mat< std::complex<in_T> >& in);
  
  
  
  template<typename in_eT>
  inline static itpp::Vec<out_eT> from(const basic_colvec< in_eT >& in);
  
  template<typename in_eT>
  inline static itpp::Vec<out_eT> from(const basic_rowvec< in_eT >& in);
  
  
  
  template<typename in_T>  
  inline static itpp::Vec<out_eT> from(const basic_colvec< std::complex<in_T> >& in);
  
  template<typename in_T>  
  inline static itpp::Vec<out_eT> from(const basic_rowvec< std::complex<in_T> >& in);
  
  
  
  template<typename in_eT, typename T1>
  inline static itpp::Vec<out_eT> from(const arma_base<in_eT,T1>& in);
  };



template<typename out_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const basic_mat<out_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in.n_cols != 1) && (in.n_rows != 1) ), "conv_to<>: given matrix can't be interpreted as a vector");
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  syslib::copy_elem(out._data(), in.mem, in.n_elem);
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const basic_mat<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( ( (in.n_cols != 1) && (in.n_rows != 1) ), "conv_to<>: given matrix can't be interpreted as a vector");
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const basic_mat< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  arma_debug_check( ( (in.n_cols != 1) && (in.n_rows != 1) ), "conv_to<>: given matrix can't be interpreted as a vector");
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const basic_colvec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const basic_colvec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const basic_rowvec<in_eT>& in)
  {
  arma_extra_debug_sigprint();
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    out_mem[i] = out_eT( in_mem[i] );
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_T>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const basic_rowvec< std::complex<in_T> >& in)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<in_T> in_eT; 
  
  itpp::Vec<out_eT> out(in.n_elem);
  
  const in_eT*  in_mem = in.memptr();
  out_eT*      out_mem = out._data();
  
  for(u32 i=0; i<in.n_elem; ++i)
    {
    copy_complex_elem(out_mem[i], in_mem[i]);
    }
  
  return out;
  }



template<typename out_eT>
template<typename in_eT, typename T1>
inline
itpp::Vec<out_eT>
conv_to< itpp::Vec<out_eT> >::from(const arma_base<in_eT,T1>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.get_ref());

  return conv_to< itpp::Vec<out_eT> >::from( tmp.M );
  }


//! @}
