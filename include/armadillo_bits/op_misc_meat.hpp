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



//! \addtogroup op_misc
//! @{


template<typename T1>
inline
void
op_log::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_log>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::log(A_mem[i]);
    }
  }
  
  

template<typename T1>
inline
void
op_log::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_log>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::log(A_mem[i]);
    }
  }
  
  

template<typename T1>
inline
void
op_trunc_log::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trunc_log>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = trunc_log(A_mem[i]);
    }
  }
  
  

template<typename T1>
inline
void
op_trunc_log::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_trunc_log>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = trunc_log(A_mem[i]);
    }
  }
  
  

template<typename T1>
inline
void
op_log10::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_log10>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::log10(A_mem[i]);
    }
  }
  
  

template<typename T1>
inline
void
op_log10::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_log10>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::log10(A_mem[i]);
    }
  }



template<typename T1>
inline
void
op_exp::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_exp>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::exp(A_mem[i]);
    }
  }
  
  

template<typename T1>
inline
void
op_exp::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_exp>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::exp(A_mem[i]);
    }
  }
  
  

template<typename T1>
inline
void
op_trunc_exp::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_trunc_exp>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = trunc_exp(A_mem[i]);
    }
  }



template<typename T1>
inline
void
op_trunc_exp::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_trunc_exp>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = trunc_exp(A_mem[i]);
    }
  }



template<typename T1>
inline
void
op_sqrt::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sqrt>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::sqrt(A_mem[i]);
    }
  }



template<typename T1>
inline
void
op_sqrt::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_sqrt>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::sqrt(A_mem[i]);
    }
  }



template<typename T1>
inline
void
op_square::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_square>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const eT val = A_mem[i];
    out_mem[i] = val*val;
    }
  }



template<typename T1>
inline
void
op_square::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_square>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const eT val = A_mem[i];
    out_mem[i] = val*val;
    }
  }



template<typename T1, typename T2>
inline
T1 
op_pow::internal_pow(const T1 base, const T2 exponent)
  {
  return std::pow(base, exponent);
  }



template<typename T2>
inline
char
op_pow::internal_pow(const char base, const T2 exponent)
  {
  typedef char out_type;
  return out_type( std::pow(double(base), exponent) );
  }



template<typename T2>
inline
unsigned char
op_pow::internal_pow(const unsigned char base, const T2 exponent)
  {
  typedef unsigned char out_type;
  return out_type( std::pow(double(base), exponent) );
  }



template<typename T2>
inline
int
op_pow::internal_pow(const int base, const T2 exponent)
  {
  typedef int out_type;
  return out_type( std::pow(double(base), exponent) );
  }



template<typename T2>
inline
unsigned int
op_pow::internal_pow(const unsigned int base, const T2 exponent)
  {
  typedef unsigned int out_type;
  return out_type( std::pow(double(base), exponent) );
  }



template<typename T1>
inline
void
op_pow::apply(Mat<typename T1::pod_type>& out, const Op<T1,op_pow>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    //out_mem[i] = std::pow(A_mem[i], in.aux);
    out_mem[i] = op_pow::internal_pow(A_mem[i], in.aux);
    }
  
  }



template<typename T1>
inline
void
op_pow::apply(Cube<typename T1::pod_type>& out, const OpCube<T1,op_pow>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(in.m);
  const Cube<eT>& A         = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    //out_mem[i] = std::pow(A_mem[i], in.aux);
    out_mem[i] = op_pow::internal_pow(A_mem[i], in.aux);
    }
  
  }



template<typename T1>
inline
void
op_pow::apply(Mat< std::complex<typename T1::pod_type> >& out, const Op<T1,op_pow>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type  T;
  typedef std::complex<T>       eT;
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  //out.set_size(A.n_rows, A.n_cols);
  out.copy_size(A);
  
  eT* out_ptr = out.memptr();
  
  
  if(in.aux.imag() == T(0))
    {
    const T in_aux_real = in.aux.real();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = std::pow(A.mem[i], in_aux_real);
      }
    }
  else
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = std::pow(A.mem[i], in.aux);
      }
    }
  }



template<typename T1>
inline
void
op_pow::apply(Cube< std::complex<typename T1::pod_type> >& out, const OpCube<T1,op_pow>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::pod_type  T;
  typedef std::complex<T>       eT;
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  const unwrap_cube<T1> tmp(in.m);
  
  const Cube<eT>& A = tmp.M;
  const u32 n_elem  = A.n_elem;
  
  //out.set_size(A.n_rows, A.n_cols, A.n_slices);
  out.copy_size(A);
  
  eT* out_ptr = out.memptr();
  
  
  if(in.aux.imag() == T(0))
    {
    const T in_aux_real = in.aux.real();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = std::pow(A.mem[i], in_aux_real);
      }
    }
  else
    {
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = std::pow(A.mem[i], in.aux);
      }
    }
  }



template<typename T1>
inline
T1 
op_pow_s32::internal_pow(const T1 base, const int exponent)
  {
  return std::pow(base, exponent);
  }



inline
char
op_pow_s32::internal_pow(const char base, const int exponent)
  {
  typedef char out_type;
  return out_type( std::pow(double(base), exponent) );
  }



inline
unsigned char
op_pow_s32::internal_pow(const unsigned char base, const int exponent)
  {
  typedef unsigned char out_type;
  return out_type( std::pow(double(base), exponent) );
  }



inline
int
op_pow_s32::internal_pow(const int base, const int exponent)
  {
  typedef int out_type;
  return out_type( std::pow(double(base), exponent) );
  }



inline
unsigned int
op_pow_s32::internal_pow(const unsigned int base, const int exponent)
  {
  typedef unsigned int out_type;
  return out_type( std::pow(double(base), exponent) );
  }



template<typename T1>
inline
void
op_pow_s32::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_pow_s32>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  const int exponent = (in.aux_u32_b == 0) ? in.aux_u32_a : -in.aux_u32_a;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    //out_mem[i] = std::pow(A_mem[i], exponent);  // causes problems with gcc 4.1/4.2 for base that has an integer type
    out_mem[i] = op_pow_s32::internal_pow(A_mem[i], exponent);
    }
  
  }



template<typename T1>
inline
void
op_pow_s32::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_pow_s32>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;
  
        eT* out_mem = out.memptr();  
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  const int exponent = (in.aux_u32_b == 0) ? in.aux_u32_a : -in.aux_u32_a;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    //out_mem[i] = std::pow(A_mem[i], exponent);  // causes problems with gcc 4.1/4.2 for base that has an integer type
    out_mem[i] = op_pow_s32::internal_pow(A_mem[i], exponent);
    }
  
  }



template<typename T1>
inline
void
op_conj::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_conj>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_write<T1> tmp(out, in.m);
  const Mat<eT>& A     = tmp.M;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::conj(A_mem[i]);
    }
  }



template<typename T1>
inline
void
op_conj::apply(Cube<typename T1::elem_type>& out, const OpCube<T1,op_conj>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_cube_write<T1> tmp(out, in.m);
  const Cube<eT>& A         = tmp.M;

        eT* out_mem = out.memptr();
  const eT* A_mem   = A.memptr();
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = std::conj(A_mem[i]);
    }
  }



//! @}
