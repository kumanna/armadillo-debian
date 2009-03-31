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


//! \addtogroup gemm
//! @{



//! \brief
//! Partial emulation of ATLAS/BLAS gemm(), using caching for speedup.
//! Matrix 'C' is assumed to have been set to the correct size (i.e. taking into account transposes)

template<const bool do_trans_A=false, const bool do_trans_B=false, const bool use_alpha=false, const bool use_beta=false>
class gemm_emul_cache
  {
  public:
  
  template<typename eT>
  inline
  static
  void
  apply
    (
          basic_mat<eT>& C,
    const basic_mat<eT>& A,
    const basic_mat<eT>& B,
    const eT alpha = eT(1),
    const eT beta  = eT(0)
    )
    {
    arma_extra_debug_sigprint();
    
    const u32 A_n_rows = A.n_rows;
    const u32 A_n_cols = A.n_cols;
    
    const u32 B_n_rows = B.n_rows;
    const u32 B_n_cols = B.n_cols;
    
    if( (do_trans_A == false) && (do_trans_B == false) )
      {
      podarray<eT> tmp(A_n_cols);
      eT* A_rowdata = tmp.memptr();
      
      for(u32 row_A=0; row_A < A_n_rows; ++row_A)
        {
        
        for(u32 col_A=0; col_A < A_n_cols; ++col_A)
          {
          A_rowdata[col_A] = A.at(row_A,col_A);
          }
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          const eT* B_coldata = B.colptr(col_B);
          
          eT acc = eT(0);
          for(u32 i=0; i < B_n_rows; ++i)
            {
            acc += A_rowdata[i] * B_coldata[i];
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(row_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(row_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(row_A,col_B) = acc + beta*C.at(row_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(row_A,col_B) = alpha*acc + beta*C.at(row_A,col_B);
            }
          
          }
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) )
      {
      for(u32 col_A=0; col_A < A_n_cols; ++col_A)
        {
        // col_A is interpreted as row_A when storing the results in matrix C
        
        const eT* A_coldata = A.colptr(col_A);
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          const eT* B_coldata = B.colptr(col_B);
          
          eT acc = eT(0);
          for(u32 i=0; i < B_n_rows; ++i)
            {
            acc += A_coldata[i] * B_coldata[i];
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,col_B) = acc + beta*C.at(col_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,col_B) = alpha*acc + beta*C.at(col_A,col_B);
            }
          
          }
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) )
      {
      basic_mat<eT> B_tmp = trans(B);
      gemm_emul_cache<false, false, use_alpha, use_beta>::apply(C, A, B_tmp, alpha, beta);
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) )
      {
      // mat B_tmp = trans(B);
      // dgemm_arma<true, false,  use_alpha, use_beta>::apply(C, A, B_tmp, alpha, beta);
      
      
      // By using the trans(A)*trans(B) = trans(B*A) equivalency,
      // transpose operations are not needed
      
      podarray<eT> tmp(B.n_cols);
      eT* B_rowdata = tmp.memptr();
      
      for(u32 row_B=0; row_B < B_n_rows; ++row_B)
        {
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          B_rowdata[col_B] = B.at(row_B,col_B);
          }
        
        for(u32 col_A=0; col_A < A_n_cols; ++col_A)
          {
          const eT* A_coldata = A.colptr(col_A);
          
          eT acc = eT(0);
          for(u32 i=0; i < A_n_rows; ++i)
            {
            acc += B_rowdata[i] * A_coldata[i];
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,row_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,row_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,row_B) = acc + beta*C.at(col_A,row_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,row_B) = alpha*acc + beta*C.at(col_A,row_B);
            }
          
          }
        }
      
      }
    }
    
  };



//! Partial emulation of ATLAS/BLAS gemm(), non-cached version.
//! Matrix 'C' is assumed to have been set to the correct size (i.e. taking into account transposes)
template<const bool do_trans_A=false, const bool do_trans_B=false, const bool use_alpha=false, const bool use_beta=false>
class gemm_emul_simple
  {
  public:
  
  template<typename eT>
  inline
  static
  void
  apply
    (
          basic_mat<eT>& C,
    const basic_mat<eT>& A,
    const basic_mat<eT>& B,
    const eT alpha = eT(1),
    const eT beta  = eT(0)
    )
    {
    arma_extra_debug_sigprint();
    
    const u32 A_n_rows = A.n_rows;
    const u32 A_n_cols = A.n_cols;
    
    const u32 B_n_rows = B.n_rows;
    const u32 B_n_cols = B.n_cols;
    
    if( (do_trans_A == false) && (do_trans_B == false) )
      {
      for(u32 row_A = 0; row_A < A_n_rows; ++row_A)
        {
        for(u32 col_B = 0; col_B < B_n_cols; ++col_B)
          {
          const eT* B_coldata = B.colptr(col_B);
          
          eT acc = eT(0);
          for(u32 i = 0; i < B_n_rows; ++i)
            {
            acc += A.at(row_A,i) * B_coldata[i];
            }
          
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(row_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(row_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(row_A,col_B) = acc + beta*C.at(row_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(row_A,col_B) = alpha*acc + beta*C.at(row_A,col_B);
            }
          }
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == false) )
      {
      for(u32 col_A=0; col_A < A_n_cols; ++col_A)
        {
        // col_A is interpreted as row_A when storing the results in matrix C
        
        const eT* A_coldata = A.colptr(col_A);
        
        for(u32 col_B=0; col_B < B_n_cols; ++col_B)
          {
          const eT* B_coldata = B.colptr(col_B);
          
          eT acc = eT(0);
          for(u32 i=0; i < B_n_rows; ++i)
            {
            acc += A_coldata[i] * B_coldata[i];
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,col_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,col_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,col_B) = acc + beta*C.at(col_A,col_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,col_B) = alpha*acc + beta*C.at(col_A,col_B);
            }
          
          }
        }
      }
    else
    if( (do_trans_A == false) && (do_trans_B == true) )
      {
      for(u32 row_A = 0; row_A < A_n_rows; ++row_A)
        {
        for(u32 row_B = 0; row_B < B_n_rows; ++row_B)
          {
          eT acc = eT(0);
          for(u32 i = 0; i < B_n_cols; ++i)
            {
            acc += A.at(row_A,i) * B.at(row_B,i);
            }
          
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(row_A,row_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(row_A,row_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(row_A,row_B) = acc + beta*C.at(row_A,row_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(row_A,row_B) = alpha*acc + beta*C.at(row_A,row_B);
            }
          }
        }
      }
    else
    if( (do_trans_A == true) && (do_trans_B == true) )
      {
      for(u32 row_B=0; row_B < B_n_rows; ++row_B)
        {
        
        for(u32 col_A=0; col_A < A_n_cols; ++col_A)
          {
          const eT* A_coldata = A.colptr(col_A);
          
          eT acc = eT(0);
          for(u32 i=0; i < A_n_rows; ++i)
            {
            acc += B.at(row_B,i) * A_coldata[i];
            }
        
          if( (use_alpha == false) && (use_beta == false) )
            {
            C.at(col_A,row_B) = acc;
            }
          else
          if( (use_alpha == true) && (use_beta == false) )
            {
            C.at(col_A,row_B) = alpha * acc;
            }
          else
          if( (use_alpha == false) && (use_beta == true) )
            {
            C.at(col_A,row_B) = acc + beta*C.at(col_A,row_B);
            }
          else
          if( (use_alpha == true) && (use_beta == true) )
            {
            C.at(col_A,row_B) = alpha*acc + beta*C.at(col_A,row_B);
            }
          
          }
        }
      
      }
    }
    
  };



//! \brief
//! Wrapper for ATLAS/BLAS dgemm function, using template arguments to control the arguments passed to dgemm.
//! Matrix 'C' is assumed to have been set to the correct size (i.e. taking into account transposes)

template<const bool do_trans_A=false, const bool do_trans_B=false, const bool use_alpha=false, const bool use_beta=false>
class gemm
  {
  public:
  
  template<typename eT>
  inline
  static
  void
  apply_blas_type( basic_mat<eT>& C, const basic_mat<eT>& A, const basic_mat<eT>& B, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    arma_extra_debug_sigprint();
    
    if( ((A.n_elem <= 64u) && (B.n_elem <= 64u)) )
      {
      gemm_emul_simple<do_trans_A, do_trans_B, use_alpha, use_beta>::apply(C,A,B,alpha,beta);
      }
    else
      {
      #if defined(ARMA_USE_ATLAS)
        {
        arma_extra_debug_print("atlas::cblas_gemm()");
        
        atlas::cblas_gemm<eT>
          (
          atlas::CblasColMajor,
          (do_trans_A) ? atlas::CblasTrans : atlas::CblasNoTrans,
          (do_trans_B) ? atlas::CblasTrans : atlas::CblasNoTrans,
          C.n_rows,
          C.n_cols,
          (do_trans_A) ? A.n_rows : A.n_cols,
          (use_alpha) ? alpha : eT(1),
          A.mem,
          (do_trans_A) ? A.n_rows : C.n_rows,
          B.mem,
          (do_trans_B) ? C.n_cols : ( (do_trans_A) ? A.n_rows : A.n_cols ),
          (use_beta) ? beta : eT(0),
          C.memptr(),
          C.n_rows
          );
        }
      #elif defined(ARMA_USE_BLAS)
        {
        arma_extra_debug_print("blas::gemm_()");
        
        const char trans_A = (do_trans_A) ? 'T' : 'N';
        const char trans_B = (do_trans_B) ? 'T' : 'N';
        
        const int m   = C.n_rows;
        const int n   = C.n_cols;
        const int k   = (do_trans_A) ? A.n_rows : A.n_cols;
        
        const eT local_alpha = (use_alpha) ? alpha : eT(1);
        
        const int lda = (do_trans_A) ? k : m;
        const int ldb = (do_trans_B) ? n : k;
        
        const eT local_beta  = (use_beta) ? beta : eT(0);
        
        arma_extra_debug_print( arma_boost::format("blas::gemm_(): trans_A = %c") % trans_A );
        arma_extra_debug_print( arma_boost::format("blas::gemm_(): trans_B = %c") % trans_B );
        
        blas::gemm_<eT>
          (
          &trans_A,
          &trans_B,
          &m,
          &n,
          &k,
          &local_alpha,
          A.mem,
          &lda,
          B.mem,
          &ldb,
          &local_beta,
          C.memptr(),
          &m
          );
        }
      #else
        {
        gemm_emul_cache<do_trans_A, do_trans_B, use_alpha, use_beta>::apply(C,A,B,alpha,beta);
        }
      #endif
      }
    }
  
  
  
  //! immediate multiplication of matrices A and B, storing the result in C
  template<typename eT>
  inline
  static
  void
  apply( basic_mat<eT>& C, const basic_mat<eT>& A, const basic_mat<eT>& B, const eT alpha = eT(1), const eT beta = eT(0) )
    {
    if( (A.n_elem <= 64u) && (B.n_elem <= 64u) )
      {
      gemm_emul_simple<do_trans_A, do_trans_B, use_alpha, use_beta>::apply(C,A,B,alpha,beta);
      }
    else
      {
      gemm_emul_cache<do_trans_A, do_trans_B, use_alpha, use_beta>::apply(C,A,B,alpha,beta);
      }
    }
  
  
  
  inline
  static
  void
  apply
    (
          basic_mat<float>& C,
    const basic_mat<float>& A,
    const basic_mat<float>& B,
    const float alpha = float(1),
    const float beta  = float(0)
    )
    {
    gemm<do_trans_A, do_trans_B, use_alpha, use_beta>::apply_blas_type(C,A,B,alpha,beta);
    }
  
  
  
  inline
  static
  void
  apply
    (
          basic_mat<double>& C,
    const basic_mat<double>& A,
    const basic_mat<double>& B,
    const double alpha = double(1),
    const double beta  = double(0)
    )
    {
    gemm<do_trans_A, do_trans_B, use_alpha, use_beta>::apply_blas_type(C,A,B,alpha,beta);
    }
  
  
  
  inline
  static
  void
  apply
    (
          basic_mat< std::complex<float> >& C,
    const basic_mat< std::complex<float> >& A,
    const basic_mat< std::complex<float> >& B,
    const std::complex<float> alpha = std::complex<float>(1),
    const std::complex<float> beta  = std::complex<float>(0)
    )
    {
    gemm<do_trans_A, do_trans_B, use_alpha, use_beta>::apply_blas_type(C,A,B,alpha,beta);
    }
  
  
  
  inline
  static
  void
  apply
    (
          basic_mat< std::complex<double> >& C,
    const basic_mat< std::complex<double> >& A,
    const basic_mat< std::complex<double> >& B,
    const std::complex<double> alpha = std::complex<double>(1),
    const std::complex<double> beta  = std::complex<double>(0)
    )
    {
    gemm<do_trans_A, do_trans_B, use_alpha, use_beta>::apply_blas_type(C,A,B,alpha,beta);
    }
  
  };



//! @}
