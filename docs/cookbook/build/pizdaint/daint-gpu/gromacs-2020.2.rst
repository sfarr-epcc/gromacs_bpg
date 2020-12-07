::

   module load daint-gpu
   module swap PrgEnv-cray PrgEnv-gnu
   module load cray-fftw
   module load craype-accel-nvidia60
   
   cd gromacs-2020.2
   mkdir build
   cd build
   
   cmake .. -DCMAKE_INSTALL_PREFIX=${HOME}/gromacs/2020.2 -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
            -DGMX_MPI=on -DGMX_GPU=on -DGMX_SIMD=AVX2_256 -DGMX_FFT_LIBRARY=fftw3 -DGMX_HWLOC=on
   
   make -j 12
   make install
   
