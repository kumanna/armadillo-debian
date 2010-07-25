=== Contents ===

1. Installation
   1.0. Preliminaries
   1.1. Manual Installation
   1.2. Installation on Unix-like systems
   1.3. Installation on MS Windows

2. Compiling Programs and Linking
   2.0. Examples
   2.1. Compiling & Linking on Unix-like systems
   2.2. Compiling & Linking on MS Windows

3. Caveats
   3.0. Support for ATLAS
   3.1. Support for ACML and Intel MKL

4. Documentation / Reference Manual

5. Using Armadillo with IT++

6. License

7. Bug Reports

8. Credits



=== 1.0. Installation: Preliminaries ===

While this library has gone through testing, not all possible
cases have been covered yet. As such, its functionality may
not be 100% correct. If you find a bug, either in the library
or the documentation, we are interested in hearing about it.

Armadillo makes extensive use of template metaprogramming,
recursive templates and template based function overloading.
As such, C++ compilers which do not fully implement the C++
standard may not work correctly.

The functionality of Armadillo is partly dependent on other
libraries -- mainly LAPACK and BLAS. Armadillo can work without
LAPACK or BLAS, but its functionality will be reduced.
In particular, basic functionality will be available
(e.g. matrix addition and multiplication), but things like
eigen decomposition or will not be. Matrix multiplication
(mainly for big matrices) may not be as fast.

For manual installation on all systems, see section 1.1.

For installation on Unix-like systems, see section 1.2.
Unix-like systems include:
  - Linux
  - Mac OS X
  - FreeBSD
  - Solaris
  - CygWin

For installation on MS Windows, see section 1.3.



=== 1.1. Manual Installation ===

The manual installation is comprised of 3 steps:

* Step 1:
  Copy the entire "include" folder to a convenient location
  and tell your compiler to use that location for header files
  (in addition to the locations it uses already).
  Alternatively, you can use the "include" folder directly.

* Step 2:
  Modify "include/armadillo_bits/config.hpp" to indicate 
  which libraries are currently available on your system.
  For example, if you have LAPACK and BLAS present, 
  uncomment the following lines:
  
  #define ARMA_USE_LAPACK
  #define ARMA_USE_BLAS

* Step 3:
  If you have LAPACK and/or BLAS present, configure your 
  compiler to link with these libraries. 
  
  You can also link with the the equivalent of LAPACK and BLAS,
  e.g. Intel's MKL or AMD's ACML. Under Mac OS X, link using 
  -framework Accelerate



=== 1.2. Installation on Unix-like systems ===

If you have installed Armadillo using an RPM or DEB package,
you don't need to do anything else. Otherwise read on.

You can use the manual installation process as described in
section 1.1, or the following CMake based automatic installation.

* Step 1:
  If CMake is not already be present on your system, download
  it from http://www.cmake.org

  On major Linux systems (such as Fedora, Ubuntu, Debian, etc),
  cmake is available as a pre-built package, though it may need
  to be explicitly installed (using a tool such as PackageKit,
  yum, rpm, apt, aptitude, etc).

* Step 2:
  If you have BLAS and/or LAPACK, install them before installing
  Armadillo. Under Mac OS X this is not necessary.
  
  On Linux systems it is recommended that the following libraries
  are present: LAPACK, BLAS, ATLAS and Boost. LAPACK and BLAS are
  the most important. If you have ATLAS and Boost, it's also necessary
  to have the corresponding header files installed.

* Step 3a:
  Open a shell (command line), change into the directory that was
  created by unpacking the armadillo archive, and type the following
  commands:

  cmake .
  make 

  The full stop separated from "cmake" by a space is important.
  CMake will figure out what other libraries are currently installed
  and will modify Armadillo's configuration correspondingly.
  CMake will also generate a run-time armadillo library, which is a 
  combined alias for all the relevant libraries present on your system
  (e.g. BLAS, LAPACK and ATLAS).
  
  If you need to re-run cmake, it's a good idea to first delete the
  "CMakeCache.txt" file (not "CMakeLists.txt").

* Step 3b:
  If you have access to administrator privileges (e.g. root),
  first enable the privileges (e.g. through "su" or "sudo")
  and then type the following command:
  
  make install

  If you don't have administrator privileges, 
  type the following command:

  make install DESTDIR=my_usr_dir

  where "my_usr_dir" is for storing C++ headers and library files.
  Make sure your C++ compiler is configured to use the sub-directories
  present within this directory.



=== 1.3. Installation on MS Windows ===

There is currently no automatic installation for Windows.
Please use the manual installation process described in
section 1.1.

Please contact the authors if you'd like to contribute
a CMake based installation solution for Windows.

Pre-compiled BLAS and LAPACK libraries for Windows are
provided in the "examples/libs_win32" folder.
They can be also obtained from:

  - http://www.stanford.edu/~vkl/code/libs.html
  - http://icl.cs.utk.edu/lapack-for-windows/lapack/
  - http://www.fi.muni.cz/~xsvobod2/misc/lapack/



=== 2.0. Compiling Programs and Linking: Examples ===

The "examples" directory contains several quick example programs
that use the Armadillo library. If Armadillo was installed manually
(i.e. according to section 1.1), you will also need to explicitly
link your programs with the libraries that were specified in
"include/armadillo_bits/config.hpp".

"example1.cpp" doesn't need any external libraries.
"example2.cpp" requires the LAPACK library or its equivalent
(e.g. the Accelerate framework on Mac OS X). You may get errors
at compile or run time if LAPACK functions are not available.

NOTE: as Armadillo is a template library, we recommended that
      optimisation is enabled during compilation.



=== 2.1. Compiling & Linking on Unix-like systems ===

Please see "examples/Makefile", which may may need to be configured
for your system. If Armadillo header files were installed in a
non-standard location, you will need to modify "examples/Makefile"
to tell the compiler where they are.

If Armadillo was installed manually and you specified that
LAPACK and BLAS are available, instead of using "-larmadillo",
use the following:
  - under Linux, use "-llapack -lblas"
  - under Mac OS X, use "-framework Accelerate"
  - under the Sun Studio compiler, try "-library=sunperf"

NOTE: on Ubuntu and Debian based systems you may need to add 
      "-lgfortran" to the compiler flags.



=== 2.2. Compiling & Linking on MS Windows ===

As a courtesy, we've provided pre-compiled 32 bit versions of
LAPACK and BLAS for Windows, as well as MSVC project files to
compile example1.cpp and example2.cpp. The project files are
stored in the following folders:
  examples/example1_win32
  examples/example2_win32

If you're not using MSVC, you will need to manually modify 
"include/armadillo_bits/config.hpp" to enable the use of
LAPACK and BLAS. Please see sections 1.1 & 1.3 for more info.

The MSCV project files were tested on Windows XP (32 bit) with
Visual C++ 2008 (Express Edition). You may need to make adaptations
for 64 bit systems, later versions of Windows and/or the compiler.

To preserve our sanity, we (Armadillo developers) don't use Windows
on a regular basis, and as such can't help you with the adaptations.

If you encounter issues with the MS Visual C++ compiler,
the following high-quality compilers are useful alternatives:

  - Intel's C++ compiler
    http://software.intel.com/en-us/intel-compilers/

  - GCC (part MinGW)
    http://www.mingw.org/

  - GCC (part of CygWin)
    http://www.cygwin.com/

If using Intel's C++ compiler, you'll need version 10.0 or better.
If using GCC, you'll need version 4.0 or better.



=== 3.0. Caveats: Support for ATLAS ===

Armadillo can use the ATLAS library for faster versions of
certain LAPACK and BLAS functions. Not all ATLAS functions are
currently used, and as such LAPACK should still be installed.

The minimum recommended version of ATLAS is 3.8.
Old versions (e.g. 3.6) can produce incorrect results
as well as corrupting memory, leading to random crashes.

Users of Ubuntu and Debian based systems should explicitly
check that version 3.6 is not installed. It's better to
remove the old version and use the standard LAPACK library.



=== 3.1. Caveats: Support for ACML and Intel MKL ===

Armadillo can work with AMD Core Math Library and Intel's
Math Kernel Library (MKL), however there are several caveats.

On Linux systems, ACML and MKL are typically installed in a
non-standard location, which can cause problems during linking.

Before installing Armadillo, the system should know where the ACML or MKL
libraries are located (for example, "/opt/intel/mkl/10.2.2.025/lib/em64t/").
This can be achieved by setting the LD_LIBRARY_PATH environment variable,
or, for a more permanent solution, adding the location of the libraries
to "/etc/ld.so.conf". It may also be possible to store a text file 
with the location in the "/etc/ld.so.conf.d" directory.
In the latter two cases you will need to run "ldconfig" afterwards.

The default installations of ACML 4.4.0 and MKL 10.2.2.025 are known 
to have issues with SELinux, which is turned on by default in Fedora
(and possibly RHEL). The problem may manifest itself during run-time,
where the run-time linker reports permission problems.
It is possible to work around the problem by applying an appropriate
SELinux type to all ACML and MKL libraries.

If you have ACML or MKL installed and they are persistently giving
you problems during linking, you can disable the support for them
by editing the "CMakeLists.txt" file, deleting "CMakeCache.txt" and
re-running the CMake based installation.  Specifically, comment out
the lines containing:
  INCLUDE(ARMA_FindMKL)
  INCLUDE(ARMA_FindACMLMP)
  INCLUDE(ARMA_FindACML)



=== 4. Documentation / Reference Manual ===

A reference manual (user documentation) is available at
http://arma.sourceforge.net or in the "docs_user" directory.
Use a web browser to open the "docs_user/index.html" file.

The user documentation explains how to use Armadillo's classes
and functions, with snippets of example code.

The technical documentation (produced with the aid of Doxygen) is
available in the "docs_tech" directory. Use a web browser to open
the "docs_tech/index.html" file. The technical documentation is only 
useful if you want to understand the internals of Armadillo.



=== 5. Using Armadillo with IT++ ===

If you wish to use the IT++ library in conjunction with Armadillo,
use #include "armadillo_itpp" instead of #include "armadillo"
in your code. See also the "examples/example_itpp.cpp" file.



=== 6. License ===

Please see the "LICENSE.txt" file.



=== 7. Bug Reports ===

If you find a bug, either in the library or the documentation,
we are interested in hearing about it. Please send a report to:
  Conrad Sanderson <conradsand at ieee dot org>
or post a message on the discussion board:
  http://sourceforge.net/apps/phpbb/arma/



=== 8. Credits ===

Main developers:
- Conrad Sanderson - http://itee.uq.edu.au/~conrad/
- Ian Cullinan
- Dimitrios Bouzas

Contributors:
- Eric R. Anderson
- Benoît Bayol
- Salim Bcoin
- Justin Bedo
- Darius Braziunas
- Ted Campbell
- Chris Davey
- Dirk Eddelbuettel
- Romain Francois
- Charles Gretton
- Edmund Highcock
- Kshitij Kulshreshtha
- Oka Kurniawan
- David Lawrence
- Carlos Mendes
- Artem Novikov
- Martin Orlob
- Ken Panici
- Adam Piątyszek
- Jayden Platell
- Vikas Reddy
- Ola Rinta-Koski
- Laurianne Sitbon
- Paul Torfs
- Yong Kang Wong
