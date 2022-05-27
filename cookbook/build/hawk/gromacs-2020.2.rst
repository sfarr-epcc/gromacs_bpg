::

  module load fftw

  cd gromacs-2020.2
  mkdir build
  cd build
  
  cmake .. -DCMAKE_INSTALL_PREFIX=${HOME}/gromacs/2020.2 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
           -DGMX_MPI=on -DGMX_SIMD=AVX2_256 -DGMX_GPU=off -DGMX_BUILD_SHARED_EXE=off -DBUILD_SHARED_LIBS=off \
	   -DGMX_FFT_LIBRARY=fftw3 -DCMAKE_PREFIX_PATH=${FFTW_ROOT}

  make -j 64
  make install
  
