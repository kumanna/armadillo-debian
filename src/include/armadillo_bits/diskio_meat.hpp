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


//! \addtogroup diskio
//! @{


template<typename eT>
inline
std::string
diskio::gen_txt_header(const basic_mat<eT>& x)
  {
  arma_type_check<diskio::is_supported_type<eT>::value == false>::apply();

  if(is_u8<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IU001");
    }
  else
  if(is_s8<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IS001");
    }
  else
  if(is_u16<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IU002");
    }
  else
  if(is_s16<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IS002");
    }
  else
  if(is_u32<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IU004");
    }
  else
  if(is_s32<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_IS004");
    }
  else
  if(is_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FN004");
    }
  else
  if(is_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FN008");
    }
  else
  if(is_complex_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FC008");
    }
  else
  if(is_complex_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_TXT_FC016");
    }
  else
    {
    return std::string();
    }
  
  }



template<typename eT>
inline
std::string
diskio::gen_bin_header(const basic_mat<eT>& x)
  {
  arma_type_check<diskio::is_supported_type<eT>::value == false>::apply();
  
  if(is_u8<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IU001");
    }
  else
  if(is_s8<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IS001");
    }
  else
  if(is_u16<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IU002");
    }
  else
  if(is_s16<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IS002");
    }
  else
  if(is_u32<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IU004");
    }
  else
  if(is_s32<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_IS004");
    }
  if(is_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FN004");
    }
  else
  if(is_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FN008");
    }
  else
  if(is_complex_float<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FC008");
    }
  else
  if(is_complex_double<eT>::value == true)
    {
    return std::string("ARMA_MAT_BIN_FC016");
    }
  else
    {
    return std::string();
    }
  
  }



inline
char
diskio::conv_to_hex_char(const u8 x)
  {
  char out;
  switch(x)
    {
    case  0: out = '0'; break;
    case  1: out = '1'; break;
    case  2: out = '2'; break;
    case  3: out = '3'; break;
    case  4: out = '4'; break;
    case  5: out = '5'; break;
    case  6: out = '6'; break;
    case  7: out = '7'; break;
    case  8: out = '8'; break;
    case  9: out = '9'; break;
    case 10: out = 'a'; break;
    case 11: out = 'b'; break;
    case 12: out = 'c'; break;
    case 13: out = 'd'; break;
    case 14: out = 'e'; break;
    case 15: out = 'f'; break;
    default: out = '-'; break;
    }

  return out;  
  }



inline
void
diskio::conv_to_hex(char* out, const u8 x)
  {
  const u8 a = x / 16;
  const u8 b = x - 16*a;

  out[0] = conv_to_hex_char(a);
  out[1] = conv_to_hex_char(b);
  }



//! Append a quasi-random string to the given filename.
//! The rand() function is deliberately not used,
//! as rand() has an internal state that changes
//! from call to call. Such states should not be
//! modified in scientific applications, where the
//! results should be reproducable and not affected 
//! by saving data.
inline
std::string
diskio::gen_tmp_name(const std::string& x)
  {
  const std::string* ptr_x = &x;
  const u8* ptr_ptr_x = reinterpret_cast<const u8*>(&ptr_x);
  
  const char* extra = ".tmp_";
  const u32 extra_size = 5;
  
  const u32 tmp_size = 2*sizeof(u8*) + 2*2;
  char tmp[tmp_size];
  
  u32 char_count = 0;
  
  for(u32 i=0; i<sizeof(u8*); ++i)
    {
    conv_to_hex(&tmp[char_count], ptr_ptr_x[i]);
    char_count += 2;
    }
  
  const u32 x_size = x.size();
  u8 sum = 0;
  
  for(u32 i=0; i<x_size; ++i)
    {
    sum += u8(x[i]);
    }
  
  conv_to_hex(&tmp[char_count], sum);
  char_count += 2;
  
  conv_to_hex(&tmp[char_count], u8(x_size));
  
  
  std::string out;
  out.resize(x_size + extra_size + tmp_size);
  
  
  for(u32 i=0; i<x_size; ++i)
    {
    out[i] = x[i];
    }
  
  for(u32 i=0; i<extra_size; ++i)
    {
    out[x_size + i] = extra[i];
    }
  
  for(u32 i=0; i<tmp_size; ++i)
    {
    out[x_size + extra_size + i] = tmp[i];
    }
  
  return out;
  }



//! Safely rename a file.
//! Before renaming, test if we can write to the final file.
//! This should prevent:
//! (i)  overwriting files that have been write protected,
//! (ii) overwriting directories.
inline
void
diskio::safe_rename(const std::string& old_name, const std::string& new_name)
  {
  std::fstream f(new_name.c_str(), std::fstream::out | std::fstream::app);
  f.put(' ');
  
  const bool writing_problem = (f.good() == false);
  f.close();
  
  arma_warn( writing_problem, "trouble writing ", new_name );
  
  if(writing_problem == false)
    {
    std::remove(new_name.c_str());
    
    const int mv_result = std::rename(old_name.c_str(), new_name.c_str());
    arma_warn( (mv_result != 0), "trouble writing ", new_name );
    }
  
  }



//! Save a matrix as raw text (no header, human readable).
//! Non-complex matrices can be loaded in Matlab and Octave.
template<typename eT>
inline
void
diskio::save_raw_ascii(const basic_mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out);
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", tmp_name);
    }
  else
    {
    u32 cell_width;
    
    // TODO: need sane values for complex numbers
    
    if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
      {
      f.setf(ios::scientific);
      f.precision(8);
      cell_width = 16;
      }
    
    for(u32 row=0; row != x.n_rows; ++row)
      {
      for(u32 col=0; col != x.n_cols; ++col)
        {
        f.put(' ');
        
        if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
          {
          f.width(cell_width);
          }
        
        f << x.at(row,col);
        }
        
      f.put('\n');
      }
    
    const bool writing_problem = (f.good() == false);
    
    arma_warn(writing_problem, "trouble writing ", tmp_name );
    
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_ascii(const basic_mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str());
  
  diskio::save_arma_ascii(x, tmp_name, f);
  
  const bool writing_problem = (f.good() == false);
  
  f.flush();
  f.close();
  
  arma_warn( writing_problem, "trouble writing ", tmp_name );
  
  if(writing_problem == false)
    {
    diskio::safe_rename(tmp_name, final_name);
    }
  }



//! Save a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void 
diskio::save_arma_ascii(const basic_mat<eT>& x, const std::string& name, std::ofstream& f)
  {
  arma_extra_debug_sigprint();
  
  if(f.is_open() == false)
    {
    arma_debug_print("unable to write ", name);
    }
  else
    {
    const ios::fmtflags orig_flags = f.flags();
    
    f << diskio::gen_txt_header(x) << '\n';
    f << x.n_rows << ' ' << x.n_cols << '\n';
    
    u32 cell_width;
    
    // TODO: need sane values for complex numbers
    
    if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )
      {
      f.setf(ios::scientific);
      f.precision(8);
      cell_width = 16;
      }
      
    for(u32 row=0; row != x.n_rows; ++row)
      {
      for(u32 col=0; col != x.n_cols; ++col)
        {
        f.put(' ');
        
        if( (is_float<eT>::value == true) || (is_double<eT>::value == true) )        
          {
          f.width(cell_width);
          }
        
        f << x.at(row,col);
        }
      
      f.put('\n');
      }
    
    f.flags(orig_flags);
    }
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_binary(const basic_mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::ofstream f(tmp_name.c_str(), std::fstream::binary);
  
  diskio::save_arma_binary(x, tmp_name, f);
  
  const bool writing_problem = (f.good() == false);
  
  f.flush();
  f.close();
  
  arma_warn( writing_problem, "trouble writing ", tmp_name );
  
  if(writing_problem == false)
    {
    diskio::safe_rename(tmp_name, final_name);
    }
  }



//! Save a matrix in binary format,
//! with a header that stores the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::save_arma_binary(const basic_mat<eT>& x, const std::string& name, std::ofstream& f)
  {
  arma_extra_debug_sigprint();
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", name);
    }
  else
    {
    f << diskio::gen_bin_header(x) << '\n';
    f << x.n_rows << ' ' << x.n_cols << '\n';
    
    f.write(reinterpret_cast<const char*>(x.mem), x.n_elem*sizeof(eT));
    }
  
  }


//
// TODO:
// add functionality to save the image in a normalised format,
// i.e. scaled so that every value falls in the [0,255] range.

//! Save a matrix as a PGM greyscale image
template<typename eT>
inline
void
diskio::save_pgm_binary(const basic_mat<eT>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  
  std::fstream f(tmp_name.c_str(), std::fstream::out | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_print("unable to write ", tmp_name);
    }
  else
    {
    f << "P5" << '\n';
    f << x.n_cols << ' ' << x.n_rows << '\n';
    f << 255 << '\n';
    
    const u32 n_elem = x.n_rows * x.n_cols;
    podarray<u8> tmp(n_elem);
    
    u32 i = 0;
    
    for(u32 row=0; row != x.n_rows; ++row)
      {
      for(u32 col=0; col != x.n_cols; ++col)
        {
        tmp[i] = u8( x(row,col) );  // TODO: add round() ?
        ++i;
        }
      }
    
    f.write(reinterpret_cast<const char*>(tmp.mem), n_elem);
    
    const bool writing_problem = (f.good() == false);
    
    arma_warn(writing_problem, "trouble writing ", tmp_name );
    
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    }
  
  }



//! Save a matrix as a PGM greyscale image
template<typename T>
inline
void
diskio::save_pgm_binary(const basic_mat< std::complex<T> >& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  const uchar_mat tmp = conv_to<uchar_mat>::from(x);
  diskio::save_pgm_binary(tmp,name);
  }



//! Load a matrix as raw text (no header, human readable).
//! Can read matrices saved as text in Matlab and Octave.
//! NOTE: this is much slower than reading a file with a header.
template<typename eT>
inline
void
diskio::load_raw_ascii(basic_mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();

  std::fstream f;
  f.open(name.c_str(), std::fstream::in);
  
  bool load_okay = true;
  
  if(f.is_open() == false)
    {
    arma_print("unable to read ", name );
    load_okay = false;
    }
  else
    {
    //std::fstream::pos_type start = f.tellg();
    
    //
    // work out the size
    
    
    u32 f_n_rows = 0;
    u32 f_n_cols = 0;
    
    bool f_n_cols_found = false;
    
    std::string line_string;
    std::string token;
    
    while( (f.good() == true) && (load_okay == true) )
      {
      std::getline(f, line_string);
      if(line_string.size() == 0)
        break;
      
      std::stringstream line_stream(line_string);
      
      u32 line_n_cols = 0;
      while (line_stream >> token)
        line_n_cols++;
      
      if(f_n_cols_found == false)
        {
        f_n_cols = line_n_cols;
        f_n_cols_found = true;
        }
      else
        {
        if(line_n_cols != f_n_cols)
          {
          arma_print("inconsistent number of columns in ", name );
          load_okay = false;
          }
        }
      
      ++f_n_rows;
      }
      
    if(load_okay == true)
      {
      f.clear();
      f.seekg(0, ios::beg);
      //f.seekg(start);
      
      x.set_size(f_n_rows, f_n_cols);
    
      eT val;
      
      for(u32 row=0; row != x.n_rows; ++row)
        {
        for(u32 col=0; col != x.n_cols; ++col)
          {
          // f >> token;
          // x.at(row,col) = eT( strtod(token.c_str(), 0) );
          
          f >> val;
          x.at(row,col) = val;
          }
        }
      }
    
    if(f.good() == false)
      {
      arma_print("trouble reading ", name );
      load_okay = false; 
      }
    
    f.close();
    }
  
  
  if(load_okay == false)
    {
    x.reset();
    }
  
  }



//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_ascii(basic_mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f(name.c_str());
  diskio::load_arma_ascii(x, name, f);
  f.close();
  }
  


//! Load a matrix in text format (human readable),
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_ascii(basic_mat<eT>& x, const std::string& name, std::ifstream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  if(f.is_open() == false)
    {
    arma_print("unable to read ", name);
    load_okay = false;
    }
  else
    {
    std::string f_header;
    u32 f_n_rows;
    u32 f_n_cols;
    
    f >> f_header;
    f >> f_n_rows;
    f >> f_n_cols;
    
    if(f_header == diskio::gen_txt_header(x))
      {
      x.set_size(f_n_rows, f_n_cols);
      
      for(u32 row=0; row != x.n_rows; ++row)
        {
        for(u32 col=0; col != x.n_cols; ++col)
          {
          f >> x.at(row,col);
          }
        }
      
      if(f.good() == false)
        {
        arma_print("trouble reading ", name);
        load_okay = false;
        }
      }
    else
      {
      arma_print("incorrect header in ", name );
      load_okay = false;
      }
  
    }
  
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Load a matrix in binary format,
//! with a header that indicates the matrix type as well as its dimensions
template<typename eT>
inline
void
diskio::load_arma_binary(basic_mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::ifstream f;
  f.open(name.c_str(), std::fstream::binary);
  diskio::load_arma_binary(x, name, f);
  f.close();
  }



template<typename eT>
inline
void
diskio::load_arma_binary(basic_mat<eT>& x, const std::string& name, std::ifstream& f)
  {
  arma_extra_debug_sigprint();
  
  bool load_okay = true;
  
  if(f.is_open() == false)
    {
    arma_print("unable to read ", name);
    load_okay = false;
    }
  else
    {
    std::string f_header;
    u32 f_n_rows;
    u32 f_n_cols;
    
    f >> f_header;
    f >> f_n_rows;
    f >> f_n_cols;
    
    if(f_header == diskio::gen_bin_header(x))
      {
      //f.seekg(1, ios::cur);  // NOTE: this may not be portable, as on a Windows machine a newline could be two characters
      f.get();
      
      x.set_size(f_n_rows,f_n_cols);
      f.read( reinterpret_cast<char *>(x.memptr()), x.n_elem*sizeof(eT));
      
      if(f.good() == false)
        {
        arma_print("trouble reading ", name);
        load_okay = false;
        }
      }
    else
      {
      arma_print("incorrect header in ", name);
      load_okay = false;
      }
    
    }
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



inline
void
diskio::pnm_skip_comments(std::fstream& f)
  {
  while( isspace(f.peek()) )
    {
    while( isspace(f.peek()) )
      f.get();
  
    if(f.peek() == '#')
      {
      while( (f.peek() != '\r') && (f.peek()!='\n') )
        f.get();
      }
    }
  }



//! Load a PGM greyscale image as a matrix
template<typename eT>
inline
void
diskio::load_pgm_binary(basic_mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  bool load_okay = true;
  
  if(f.is_open() == false)
    {
    arma_print("unable to read ", name );
    load_okay = false;
    }
  else
    {
    std::string f_header;
    f >> f_header;
    
    if(f_header == "P5")
      {
      u32 f_n_rows   = 0;
      u32 f_n_cols   = 0;
      int f_maxval = 0;
    
      diskio::pnm_skip_comments(f);
    
      f >> f_n_cols;
      diskio::pnm_skip_comments(f);
    
      f >> f_n_rows;
      diskio::pnm_skip_comments(f);
    
      f >> f_maxval;
      f.get();
      
      if( (f_maxval > 0) || (f_maxval <= 65535) )
        {
        x.set_size(f_n_rows,f_n_cols);
        
        if(f_maxval <= 255)
          {
          const u32 n_elem = f_n_cols*f_n_rows;
          podarray<u8> tmp(n_elem);
          
          f.read( reinterpret_cast<char*>(tmp.memptr()), n_elem);
          
          u32 i = 0;
          
          //cout << "f_n_cols = " << f_n_cols << endl;
          //cout << "f_n_rows = " << f_n_rows << endl;
          
          
          for(u32 row=0; row != f_n_rows; ++row)
            {
            for(u32 col=0; col != f_n_cols; ++col)
              {
              x.at(row,col) = eT(tmp[i]);
              ++i;
              }
            
            }
          }
        else
          {
          const u32 n_elem = f_n_cols*f_n_rows;
          podarray<u16> tmp(n_elem);
          
          f.read( reinterpret_cast<char *>(tmp.memptr()), n_elem*2);
          
          u32 i = 0;
          
          for(u32 row=0; row != f_n_rows; ++row)
            {
            for(u32 col=0; col != f_n_cols; ++col)
              {
              x.at(row,col) = eT(tmp[i]);
              ++i;
              }
            
            }
          
          }
        
        }
      
      if(f.good() == false)
        {
        arma_print("trouble reading ", name);
        load_okay = false;
        }
      }
    else
      {
      arma_print("unsupported header in ", name);
      load_okay = false;
      }
    
    f.close();
    }
  
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! Load a PGM greyscale image as a matrix
template<typename T>
inline
void
diskio::load_pgm_binary(basic_mat< std::complex<T> >& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  uchar_mat tmp;
  tmp.load(name);
  x = conv_to< basic_mat< std::complex<T> > >::from(tmp);
  }



//! Try to load a matrix by automatically determining its type
template<typename eT>
inline
void
diskio::load_auto_detect(basic_mat<eT>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  static const std::string ARMA_MAT_TXT = "ARMA_MAT_TXT";
  static const std::string ARMA_MAT_BIN = "ARMA_MAT_BIN";
  static const std::string           P5 = "P5";
  
  std::fstream f;
  f.open(name.c_str(), std::fstream::in | std::fstream::binary);
  
  if(f.is_open() == false)
    {
    arma_print("diskio::load_auto_detect(): unable to read ", name );
    x.reset();
    }
  else
    {
    podarray<char> raw_header(1 + ARMA_MAT_TXT.length());
    
    f.read(raw_header.memptr(), ARMA_MAT_TXT.length());
    raw_header[ARMA_MAT_TXT.length()] = '\0';
    
    const std::string header = raw_header.mem;
    
    if(ARMA_MAT_TXT == header.substr(0,ARMA_MAT_TXT.length()))
      {
      load_arma_ascii(x, name);
      }
    else
    if(ARMA_MAT_BIN == header.substr(0,ARMA_MAT_BIN.length()))
      {
      load_arma_binary(x, name);
      }
    else
    if(P5 == header.substr(0,P5.length()))
      {
      load_pgm_binary(x, name);
      }
    else
      {
      load_raw_ascii(x, name);
      }
    
    f.close();
    }
  
  }



template<typename T1>
inline
void
diskio::save_field(const field<T1>& x, const std::string& final_name)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_basic_mat<T1>::value == false>::apply();
  
  const std::string tmp_name = diskio::gen_tmp_name(final_name);
  std::ofstream f( tmp_name.c_str(), std::fstream::binary );
  
  if(f.is_open() == false)
    {
    arma_print("couldn't write ", tmp_name);
    }
  else
    {
    f << "ARMA_FLD_BIN" << '\n';
    f << x.n_rows << '\n';
    f << x.n_cols << '\n';
    
    for(u32 i=0; i<x.n_elem; ++i)
      {
      diskio::save_arma_binary(x[i], tmp_name, f);
      }
    
    const bool writing_problem = (f.good() == false);
    
    arma_warn(writing_problem, "trouble writing ", tmp_name );
    
    f.flush();
    f.close();
    
    if(writing_problem == false)
      {
      diskio::safe_rename(tmp_name, final_name);
      }
    
    }
  
  }



template<typename T1>
inline
void
diskio::load_field(field<T1>& x, const std::string& name)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_basic_mat<T1>::value == false>::apply();
  
  bool load_okay = true;
  
  std::ifstream f( name.c_str() );
  
  if(f.fail())
    {
    arma_print("unable to read ", name );
    load_okay = false;
    }
  else
    {
    std::string f_type;
    u32         f_n_rows;
    u32         f_n_cols;
    
    f >> f_type;
    if(f_type != "ARMA_FLD_BIN")
      {
      arma_print("unsupported field type in ", name);
      load_okay = false;
      }
    else
      {
      f >> f_n_rows;
      f >> f_n_cols;
      
      x.set_size(f_n_rows, f_n_cols);
      
      f.get();      
      
      for(u32 i=0; i<x.n_elem; ++i)
        {
        diskio::load_arma_binary(x[i], name, f);
        
        if(f.good() == false)
          {
          arma_print("trouble reading ", name);
          load_okay = false;
          break;
          }
        }
      }
    }
  
  f.close();
  
  
  if(load_okay == false)
    {
    x.reset();
    }
  }



//! @}

