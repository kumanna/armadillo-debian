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


//! \addtogroup ostream
//! @{

class arma_ostream
  {
  public:
  
  template<typename eT>
  inline static u32 set_flags(std::ostream& o, const basic_mat<eT>& m);
  
  template<typename T>
  inline static u32 set_flags(std::ostream& o, const basic_mat< std::complex<T> >& m);
  
  };



template<typename eT>
inline
u32
arma_ostream::set_flags(std::ostream& o, const basic_mat<eT>& m)
  {
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.fill(' ');
  
  bool use_layout_B = false;
  bool use_layout_C = false;
  
  for(u32 i=0; i<m.n_elem; ++i)
    {
    const eT val = m.mem[i];
    
    if(
      val >= eT(+100) ||
      ( (is_signed<eT>::value == true) && (val <= eT(-100)) ) ||
      ( (is_non_integral<eT>::value == true) && (val > eT(0)) && (val <= eT(+1e-4)) ) ||
      ( (is_non_integral<eT>::value == true) && (is_signed<eT>::value == true) && (val < eT(0)) && (val >= eT(-1e-4)) ) 
      )
      {
      use_layout_C = true;
      break;
      }
      
    if(
      (val >= eT(+10)) || ( (is_signed<eT>::value == true) && (val <= eT(-10)) )
      )
      {
      use_layout_B = true;
      }
    }
  
  u32 cell_width;
  
  if(use_layout_C == true)
    {
    o.setf(ios::scientific);
    o.unsetf(ios::fixed);
    o.precision(4);
    cell_width = 13;
    }
  else
  if(use_layout_B == true)
    {
    o.unsetf(ios::scientific);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 10;
    }
  else
    {
    o.unsetf(ios::scientific);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 9;
    }
 
  return cell_width;
  }



//! "better than nothing" settings for complex numbers
template<typename T>
inline
u32
arma_ostream::set_flags(std::ostream& o, const basic_mat< std::complex<T> >& m)
  {
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.fill(' ');
  
  o.setf(ios::scientific);
  o.unsetf(ios::fixed);
  o.precision(3);
  const u32 cell_width = 27;
  return cell_width;
  }



//! Print a matrix to the specified stream
template<typename eT>
inline
std::ostream&
operator<< (std::ostream& o, const basic_mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  const ios::fmtflags orig_flags = o.flags();
  const u32 cell_width = arma_ostream::set_flags(o, m);
  
  for(u32 row=0; row != m.n_rows; ++row)
    {
    for(u32 col=0; col != m.n_cols; ++col)
      {
      o.width(cell_width);
      o << m.at(row,col);
      }
    
    o << '\n';
    }
  
  o.flush();
  o.flags(orig_flags);
  
  return o;
  }





//! @}
