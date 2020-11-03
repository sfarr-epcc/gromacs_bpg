====================
Performance Cookbook
====================

This part of the GROMACS best practice guide assumes your simulations are prepared appropriately and focuses on providing guidance related to running GROMACS simulations, i.e. executing ``mdrun``.


-----------------------------------------------------------
Top tips for running ``mdrun`` and getting good performance
-----------------------------------------------------------
This section provides top tips for how best to run your simulation with ``mdrun`` so as to make good use of available hardware and/or to obtain results in the shortest time possible, be it on a laptop, a multi-GPU desktop workstation, a departmental cluster, or a large supercomputer.


*Provide key highlights based on the more detailed "getting good performance from mdrun" manual page. Link out the manual page for more information, as discussed the BPG can act in part as a guide (a bridge) to the details in the manual, as MarkA did in 2016 with dissemination of then up-to-date guidance on performance*


Single compute node, desktop workstation, or laptop
---------------------------------------------------
Can run thread-mpi (``gmx mdrun``), no need for MPI version (``gmx_mpi mdrun``) - may even be slower


++++++++
CPU only
++++++++
- experiment with thread-MPI ranks x OpenMP threads per rank by varying ``-ntmpi N`` and ``-ntomp M`` such that ``M x N =`` total number of cores
  - consider threads within socket
- hyperthreading (SMT) - may get some benefit for thousands of particles per core
  - adjust so ``M x N =`` number of virtual cores (hardware threads) rather than physical cores

+++++++++
CPU + GPU
+++++++++
- Default one thread-mpi rank per GPU (default) likely optimal
  - use ``-ntomp`` to use remaining CPU cores (similar considerations for hyperthreading/SMT as for CPU only nodes)
  - if there are more processor sockets than GPUs, probably want one rank per socket and share access to the same GPU using ``-gputasks`` mapping

GPU offload options:

- Default behaviour is to offload short-range non-bonded (``-nb gpu``), PME (``-pme gpu``), and bonded (``-bonded gpu``)
- PME offload: one rank only, and limitations (list)

  - If multiple GPUs, can use ``-gputasks`` to dedicate one GPU for PME offload rank, PP calculations on other GPUs

    - this would happen anyway if total number of ranks = number of GPUs on node
    
- Try offload options taking into account `performance considerations for GPU tasks <http://manual.gromacs.org/current/user-guide/mdrun-performance.html#performance-considerations-for-gpu-tasks>`_ given your hardware 
- FFTs offloaded if PME offloaded, but can avoid with ``-pmefft cpu`` - consider for older GPU with newer CPU
- Constraint calculations and coordinate updates default to CPU, but can be offloaded with ``-update gpu``
    


Multiple networked compute nodes in a cluster or supercomputer
--------------------------------------------------------------
Use ``gmx_mpi`` not ``gmx``, launch with parallel application launcher which determines number of ranks (not ``-ntmpi``)

- Dedicated PME ranks
- ``gmx tune_pme`` (without ``dlb``) - link out to manual, don't repeat explanation here
- Advice from MarkA still valid?: "Best when the two groups of {PP, PME} MPI ranks have sizes that are composite numbers with lots of common factors, e.g. {48,16} > {40,24} >> {42,22}"

  
++++++++
CPU only
++++++++
- Experiment with number of dedicated PME ranks (using ``gmx tune_pme``), OpenMP threading for PP ranks (can choose different threading for PME rank if deemed useful)

+++++++++
CPU + GPU
+++++++++
  - In addition to GPU offloading considerations already relevant for single node, investigate which is faster:

    - tuned multi-rank (multicore) CPU PME computation
    - single-rank offloaded PME calculation on one GPU (in entire simulation)






----------------------------------------------------------
Acting on performance-related warnings found in ``md.log``
----------------------------------------------------------
This section provides guidance on identifying, understanding, and acting on performance-related warnings and suggestions issued by ``mdrun`` that you may encounter in ``md.log``. These may suggest more efficient ways to launch ``mdrun``, or spot if the GROMACS installation you are using has not been built to give good performance on the hardware and suggest how to improve this. 


*Provide practical guidance, and link out to Installation Guide in manual (again, acting as bridge / guide to the manual)*


- SIMD warning in ``md.log``, e.g.: ::

   Highest SIMD level requested by all nodes in run: AVX_512
   SIMD instructions selected at compile time:       AVX2_256
   This program was compiled for different hardware than you are running on,
   which could influence performance. This build might have been configured on a
   login node with only a single AVX-512 FMA unit (in which case AVX2 is faster),
   while the node you are running on has dual AVX-512 FMA units.

  - proviso regarding narrower SIMD width (AVX2_256) potentially better (than AVX_512) for Intel Skylake & Cascade Lake
   
- any other warnings? (e.g. CUDA related such as wrong or suboptimal sm_arch, OpenCL-related, ...?)



  
------------------------------------------------------
Building and running GROMACS on PRACE/EuroHPC machines
------------------------------------------------------
This section provides concrete recipes showing how best to build GROMACS for good performance on some of the largest EU-based supercomputers available to EU researchers through PRACE (and in future also EuroHPC), and guidance on how best to run GROMACS on these machines. 

General guidance for each machine is complemented by an analysis illustrating the effect of key runtime execution options on ``mdrun`` performance for a range of benchmark simulations, comparison with which should help you determine how best to run your own simulations on these machines. 


Benchmark simulations
---------------------
Brief description of the benchmarks 





Piz Daint (CSCS, Switzerland)
------------------------------------------------
++++++++
Hardware
++++++++
Focus on XC50 GPU partitition

Each node:

- Processor: 1 x 12-core Intel Xeon E5-2690 v3 @ 2.60GHz (one socket)
- Memory: 64GB RAM
- GPU: 1 x NVidia P100
- CRAY Aries interconnect

++++++++
Software
++++++++
- Cray MPICH MPI library
- Cray-optimised FFTW
- Cray-libsci provides BLAS & LAPACK
- craype-accel-nvidia60 targets the correct SM architecture to compile for the P100 GPU
  
+++++
Build
+++++
Build instructions Piz Daint's XC50 GPU partition: 

.... include:: cookbook/build/pizdaint/daint-gpu/gromacs-2020.2.rst
   
::

   module load daint-gpu
   module swap PrgEnv-cray PrgEnv-gnu
   module load cray-fftw
   module load craype-accel-nvidia60
   
   cd gromacs-2020.2
   mkdir build
   cd build
   
   cmake .. -DCMAKE_INSTALL_PREFIX=${HOME}/gromacs/2020.2 -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DGMX_MPI=on -DGMX_GPU=on -DGMX_SIMD=AVX2_256 -DGMX_FFT_LIBRARY=fftw3 -DGMX_HWLOC=on
   
   make -j 12
   make install

   
+++
Run
+++
Example job script to run ``mdrun`` on Piz Daint's XC50 GPU partition for 1 hour on 4 nodes, with 1 MPI rank per node and 12 OpenMP threads per rank, without hyperthreading

::

   #!/bin/bash -l
   #SBATCH --job-name=benchmark
   #SBATCH --time=01:00:00
   #SBATCH --nodes=1
   #SBATCH --ntasks-per-node=1
   #SBATCH --cpus-per-task=12
   #SBATCH --ntasks-per-core=1     # 1 = no hyperthreading, 2 = with hyperthreading
   #SBATCH --hint=nomultithread    # nomultithread = no hyperthreading, multithread = hyperthreading
   #SBATCH --partition=normal
   #SBATCH --constraint=gpu
   
   module swap PrgEnv-cray PrgEnv-gnu
   module load daint-gpu
   module load cray-fftw 
   module load craype-accel-nvidia60
   
   export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
   export CRAY_CUDA_MPS=1
   
   export PATH=${HOME}/gromacs/2020.2/bin:$PATH
   
   srun gmx_mpi mdrun -s benchmark.tpr -ntomp ${OMP_NUM_THREADS}

   
+++++++++++
Performance
+++++++++++
Keeping in mind best practice to associate one rank (one spatial domain) with each GPU means we run ``mdrun`` with one MPI rank per node since there is only GPU on each node. We choose 12 OpenMP threads (``-ntomp 12``) per MPI rank in order to use all CPU cores, or 24 threads per rank if hyperthreading is enabled. Since there is only one processor (one socket) on Piz Daint GPU compute nodes, we do not expect to suffer a performance penalty by having a single rank spanning across all 12 cores using multithreading. Using default GPU offloaded settings means PME is offloaded to GPU, which can currently only happen on a single rank, corresponding to running ``mdrun`` with ``-npme 1``. 

[Figure showing scaling of benchmark performance for these options, i.e. with single PME rank offloaded to GPU]

Disabling PME GPU offloading by choosing ``-pme cpu`` means PME computations are done on CPU cores. Consider using ``tune_pme`` to determine optimal number of PME ranks.  

[Figure showing scaling of benchmark performance with increasing number of PME ranks starting with 1 (=1 node)] 




HAWK (HLRS, Germany)
--------------------


++++++++
Hardware
++++++++
Each node:

- Processors: 2 x 64-core AMD EPYC 7742 @2.25 GHz (two sockets)
- Memory: 256GB RAM
- Interconnect: InfiniBand HDR200

++++++++
Software
++++++++
Software stack on system (available preinstalled as modules):

- HPE MPT (not compared performance with OpenMPI)
- OpenMPI
- FFTW (claims Zen2 architecture-specific build)

  
+++++
Build
+++++

::

  module load fftw

  cd gromacs-2020.2
  mkdir build
  cd build
  
  cmake .. -DCMAKE_INSTALL_PREFIX=${HOME}/gromacs/2020.2 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DGMX_MPI=on -DGMX_SIMD=AVX2_256 -DGMX_GPU=off -DGMX_BUILD_SHARED_EXE=off -DBUILD_SHARED_LIBS=off -DGMX_FFT_LIBRARY=fftw3 -DCMAKE_PREFIX_PATH=${FFTW_ROOT}

  make -j 64
  make -j 64 check
  make install
  


+++
Run
+++



+++++++++++
Performance
+++++++++++








---------------------------------------------------------------------
GROMACS Reference Benchmarks Performance on PRACE/EuroHPC machines
---------------------------------------------------------------------
This section provides a reference set of benchmark simulation performance results representative of good obtainable performance on PRACE/EuroHPC machines with GROMACS built and run according to best practice outlined in this guide.

These results are intended as a convenient reference to help researchers estimate compute time requirements for their proposed research in preparation for applying to HPC resource allocation calls.










