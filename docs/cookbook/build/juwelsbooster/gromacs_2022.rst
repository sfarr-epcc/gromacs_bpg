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
  module load Stages/2022 
  module load GCC/11.2.0
  module load CUDA/11.5
  module load OpenMPI/4.1.2
  module load CMake/3.21.1
  module load FFTW/3.3.10
  
  cmake .. -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DGMX_FFT_LIBRARY=fftw3 -DGMX_GPU=CUDA -DGMX_MPI=ON -DGMX_OPENMP=ON -DCMAKE_INSTALL_PREFIX=$(pwd)/.. -DGMX_HWLOC=ON
  
  make -j20
  
  make install