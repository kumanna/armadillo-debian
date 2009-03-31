Name:           armadillo
Version:        0.6.2
Release:        1
Summary:        Fast C++ matrix library with interfaces to LAPACK and ATLAS

Group:          Development/Libraries
License:        LGPLv3+
URL:            http://arma.sourceforge.net/
Source:         http://download.sourceforge.net/arma/%{name}-%{version}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
BuildRequires:  cmake, gcc-c++, libstdc++-devel, atlas-devel, lapack-devel, blas-devel, boost-devel

%description
Armadillo is a streamlined C++ linear algebra library (matrix maths)
aiming towards a good balance between speed and ease of use.
Integer, floating point and complex numbers are supported,
as well as a subset of trigonometric and statistics functions.
Optional integration with LAPACK and ATLAS libraries is also provided.
A delayed evaluation approach is employed (during compile time)
to combine several operations into one and reduce (or eliminate) 
the need for temporaries. This is accomplished through recursive
templates and template meta-programming.  This library is useful
if C++ has been decided as the language of choice (due to speed
and/or integration capabilities), rather than another language
like Matlab (TM) or Octave.  The library is distributed under a 
license that is useful in both open-source and commercial contexts.


%package devel
Summary:        Development headers and docs for the Armadillo C++ library
Group:          Development/Libraries
Requires:       %{name} = %{version}-%{release}
Requires:       libstdc++-devel, atlas-devel, lapack-devel, blas-devel, boost-devel

%description devel
This package contains files necessary for development using the
Armadillo C++ library. It contains header files, example programs,
user documentation (reference guide), and the technical documentation.


%prep
%setup -q


%build
cmake -DCMAKE_INSTALL_PREFIX:PATH=%{_prefix} -DINCLUDE_INSTALL_DIR:PATH=%{_includedir} -DLIB_INSTALL_DIR:PATH=%{_libdir} .
make VERBOSE=1 %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/%{_docdir}/%{name}-%{version}
cp -r -p CREDITS.txt INSTALLATION.txt LICENSE.txt LICENSE_GPL.txt LICENSE_LGPL.txt README.txt index.html examples docs_user docs_tech $RPM_BUILD_ROOT/%{_docdir}/%{name}-%{version}


%clean
rm -rf $RPM_BUILD_ROOT


%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig


%files
%defattr(-,root,root,-)
%{_libdir}/*.so.*


%files devel
%defattr(-,root,root,-)
%{_libdir}/*.so
%{_includedir}/armadillo
%{_includedir}/armadillo_bits/
%{_includedir}/armadillo_itpp
%{_docdir}/%{name}-%{version}/


%changelog
* Wed Mar 24 2009  Conrad Sanderson  <conradsand at ieee dot org>
- Added explicit dependence on libstdc++-devel 

* Wed Mar 17 2009  Conrad Sanderson  <conradsand at ieee dot org>
- Simplified specification of directories
- Removed library packages specified by "Requires", as library dependencies are detected automatically

* Wed Mar 12 2009  Conrad Sanderson  <conradsand at ieee dot org>
- Modified to generate separate devel package (subsumes previous doc package)
- Removed redundant packages specified by "BuildRequires"
- Added CMake installation prefixes to allow for x86_64

* Wed Feb  4 2009  Conrad Sanderson  <conradsand at ieee dot org>
- Modified to generate separate doc package

* Thu Jan 28 2009  Conrad Sanderson  <conradsand at ieee dot org>
- Added argument to cmake: -DCMAKE_INSTALL_PREFIX=/usr 

* Thu Jan 22 2009  Conrad Sanderson  <conradsand at ieee dot org>
- Initial spec file prepared


  
