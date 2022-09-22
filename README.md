0. PREREQUISITES
==================

    * GNU compilers for C, C++ and Fortran

    * A working MPI installation

    * Python 3.x.x to build PETSc

1. INSTALLATION
==================
  
    * Run "python install.py"                                   
    
    * Compile any joe.cpp file as "joecxx {CFLAGS} joe.cpp {LIBS}"


2. SOME POINTERS
==================

    * The turbulence models have not been included in the compilation process for the libraries
      as a lot of times changes need to be made in there. Some turbulence model files have been
      given as guidelines for users to write their own turbulence model classes in joe.cpp.
