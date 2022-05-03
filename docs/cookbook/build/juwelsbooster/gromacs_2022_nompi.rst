::
    
  #!/bin/bash
  
  # Download the source code
  wget https://ftp.gromacs.org/gromacs/gromacs-2022.tar.gz
  
  # extract
  tar xvf gromacs-2022.tar.gz
  
  # move into the source directory
  cd gromacs-2022
  
  # make the build directory
  mkdir build
  cd build
  
  
  # load the required modules
  module load GCC/11.2.0
  module load CUDA/11.5
  module load CMake/3.21.1


  module list

  cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DGMX_BUILD_OWN_FFTW=yes -DGMX_GPU=CUDA -DGMX_MPI=OFF -DGMX_OPENMP=ON -DCMAKE_INSTALL_PREFIX=$(pwd)/.. -DGMX_HWLOC=ON
  
  make -j20
  
  make install