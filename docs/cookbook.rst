====================
Performance Cookbook
====================

The performance cookbook part of the GROMACS best practice guide
assumes your simulations are prepared appropriately and provides
concrete guidance on how best to run GROMACS simulations,
i.e. execute ``mdrun``, so as to make good use of available hardware
and obtain results in the shortest time possible, be it on a laptop, a
multi-GPU desktop workstation, a departmental cluster, or a large
supercomputer.

This complements and provides a bridge into navigating the detailed
information provided in the "getting good performance from mdrun" page
in the GROMACS manual:

http://manual.gromacs.org/current/user-guide/mdrun-performance.html

As well as general guidance applicable to whatever machine you may be
running GROMACS on, this guide provides concrete examples of
performance on PRACE (and in future EuroHPC) machines.

-------------------------------------------------------------------
General guidance for running ``mdrun`` and getting good performance
-------------------------------------------------------------------


Single compute node, desktop workstation, or laptop
---------------------------------------------------

When using a laptop or desktop workstation, we typically run the
thread-mpi version of mdrun (``gmx mdrun``) rather than the fully
MPI-enabled version (``gmx_mpi mdrun``) which is essential for running
GROMACS on multi-node shared-memory systems such as HPC machines.


++++++++
CPU only
++++++++

On a single-node CPU-only machine one should experiment with
different combinations of numbers of thread-MPI ranks (controlled by
varying by varying ``-ntmpi N`` and numbers of OpenMP threads per rank
(controlled by varying ``-ntomp M``), such that ``M x N`` is equal to
the total number of CPU cores available on the node. 

Enabling usage of all logical cores (SMT - simultaneous multithreads),
also known as hyperthreads for Intel processors, may yield some
benefit. In this case ``-ntmpi N`` and ``-ntomp M`` should be chosen
such that ``M x N`` equals the number of virtual, also known as
logical cores (hardware threads) rather than physical cores.

For details see http://manual.gromacs.org/current/user-guide/mdrun-performance.html

+++++++++
CPU + GPU
+++++++++

- Default one thread-mpi rank per GPU (default) likely optimal
 
   * use ``-ntomp`` to use remaining CPU cores (similar considerations
  for hyperthreading/SMT as for CPU only nodes)
  
   * if there are more processor sockets than GPUs, probably want one
    rank per socket and share access to the same GPU using
    ``-gputasks`` mapping

GPU offload options:

- Default behaviour is to offload short-range non-bonded (``-nb
  gpu``), PME (``-pme gpu``), and bonded (``-bonded gpu``)
- PME offload: one rank only, and limitations (list)
   * If multiple GPUs, can use ``-gputasks`` to dedicate one GPU for
     PME offload rank, PP calculations on other GPUs. This would
     happen anyway if total number of ranks = number of GPUs on node
    
- Try offload options taking into account `performance considerations
  for GPU tasks
  <http://manual.gromacs.org/current/user-guide/mdrun-performance.html#performance-considerations-for-gpu-tasks>`_
  given your hardware
- FFTs offloaded if PME offloaded, but can avoid with ``-pmefft
  cpu`` - consider for older GPU with newer CPU
- Constraint calculations and coordinate updates default to CPU, but
  can be offloaded with ``-update gpu``
    


Multiple networked compute nodes in a cluster or supercomputer
--------------------------------------------------------------

Use ``gmx_mpi`` not ``gmx``, launch with parallel application launcher
(``mpirun``, ``mpiexec``, ``srun``, or ``aprun``) which determines
number of ranks (not ``-ntmpi``)

For these larger runs, GROMACS can benefit from having MPI ranks dedicated
to PME ranks and spawns a certain fraction of these based on internal settings.

You should examine the md.log file produced and check if there is any
note suggesting to use more or fewer PME ranks in future. To determine
this in a systematic way you can make use of the ``gmx tune_pme`` tool
as described in

http://manual.gromacs.org/current/onlinehelp/gmx-tune_pme.html#gmx-tune-pme


 
++++++++
CPU only
++++++++

Start with one rank per core. Once scaling of performance with
increasing numbers of cores gets bad (expected for around ~200
particles per core), try introducing multiple OpenMP threads per rank
as this may give some benefit.

Experiment with number of dedicated PME ranks (using ``gmx tune_pme``)

OpenMP threading for PME rank can be different than for PP ranks (see
    ``-ntomp_pme`` option). 


    
+++++++++
CPU + GPU
+++++++++
In addition to GPU offloading considerations already relevant for single node, investigate which is faster:

- pme-tuned tuned multi-rank (multicore) CPU-based PME computation

- single-rank offloaded PME calculation on one GPU in entire simulation



      




      
----------------------------------------------------------
Getting good GROMACS performance on PRACE/EuroHPC machines
----------------------------------------------------------

This section provides concrete recipes showing how best to build
GROMACS for good performance on some of the largest EU-based
supercomputers available to EU researchers through PRACE (and in
future also EuroHPC), and guidance on how best to run GROMACS on these
machines.

General guidance for each machine is complemented by an analysis
illustrating the effect of key runtime execution options on ``mdrun``
performance for a range of benchmark simulations. Comparison between
this and your own simulations should help you determine how best to
run your own simulations on these machines.


Benchmark simulations
---------------------

A brief description is provided below of the benchmarks used to
illustrate how to obtain good performance on PRACE/EuroHPC machines
for a range of system sizes and types .

Benchmarks prefixed with "bench" are available from the Dept. of
Theoretical and Computational Biophysics at the MPI for biophysical
Chemistry, Göttingen: https://www.mpibpc.mpg.de/grubmueller/bench

Benchmarks suffixed with "_HBS" are available from the UK's HECBioSim
consortium of computational biomolecular researchers:
https://www.hecbiosim.ac.uk/benchmarks



- **20k_HBS**:
   * 3NIR Crambin
   * Total number of atoms: 19,605
   * Protein atoms: 642
   * Water atoms: 18,963

- **benchMEM**:
   * Protein in membrane, surrounded by water
   * Total number of atoms: 82k
     

- **465k_HBS**:
   * hEGFR Dimer of 1IVO and 1NQL
   * Total number of atoms: 465,399
   * Protein atoms: 21,749
   * Lipid atoms: 134,268
   * Water atoms: 309,087
   * Ions: 295
    
- **benchRIB**:
   * Ribosome in water
   * Total number of atoms: 2M

- **benchPEP**:
   * Peptides in water
   * Total number of atoms: 12M

  


HAWK (HLRS, Germany)
--------------------
https://www.hlrs.de/systems/hpe-apollo-hawk/

HAWK is currently (November 2020) listed as number 16 on the Top500,
the 6th largest European HPC system, and is accessible through PRACE
access mechanisms.

**Hardware**

Each node has:

- Processors: 2 x 64-core AMD EPYC 7742 @2.25 GHz (two sockets)
- Memory: 256GB RAM
- Interconnect: InfiniBand HDR200

**Software**

Software stack on system (available preinstalled as modules):

- HPE MPT
- OpenMPI
- FFTW (Zen2 architecture specific build)
  
**Build**

A multi-node MPI-enabled version of GROMACS with good performance on HAWK can be built with the following options:
::

  module load fftw

  cd gromacs-2020.2
  mkdir build
  cd build
  
  cmake .. -DCMAKE_INSTALL_PREFIX=${HOME}/gromacs/2020.2 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DGMX_MPI=on -DGMX_SIMD=AVX2_256 -DGMX_GPU=off -DGMX_BUILD_SHARED_EXE=off -DBUILD_SHARED_LIBS=off -DGMX_FFT_LIBRARY=fftw3 -DCMAKE_PREFIX_PATH=${FFTW_ROOT}

  make -j 64
  make install
  


**Run**

Example job script to run ``mdrun`` on Piz Daint's XC50 GPU partition for 1 hour on 4 nodes, with 1 MPI rank per node and 12 OpenMP threads per rank, without hyperthreading

.. include:: cookbook/run/hawk/jobscript.rst



**MPI vs OpenMP hybrid parallelism and SIMT multithreading**

Examining performance of the chosen benchmarks on HAWK we observe a
very significant effect on performance of the choice of different MPI
vs OpenMP hybrid decompositions. 

Simultaneous multithreading (SMT) gives some performance boost


.. csv-table:: 20k_HBS benchmark (20k atoms) on HAWK with ``-dlb no`` and ``-tunepme no``
   :file: cookbook/results/hawk/20k_HBS_dlbNO_tunepmeNO.csv
   :header-rows: 2


.. csv-table:: benchMEM benchmark (81k atoms) on HAWK with ``-dlb no`` and ``-tunepme no``
   :file: cookbook/results/hawk/benchMEM_dlbNO_tunepmeNO.csv
   :header-rows: 2

.. csv-table:: 465k_HBS benchmark (461k atoms) on HAWK with ``-dlb no`` and ``-tunepme no``
   :file: cookbook/results/hawk/465k_HBS_dlbNO_tunepmeNO.csv
   :header-rows: 2

.. csv-table:: benchRIB benchmark (2M atoms) on HAWK with ``-dlb no`` and ``-tunepme no``
   :file: cookbook/results/hawk/benchRIB_dlbNO_tunepmeNO.csv
   :header-rows: 2

.. csv-table:: benchPEP benchmark (12M atoms) on HAWK with ``-dlb no`` and ``-tunepme no``
   :file: cookbook/results/hawk/benchPEP_dlbNO_tunepmeNO.csv
   :header-rows: 2



**PME ranks**

	 







  

Piz Daint (CSCS, Switzerland)
-----------------------------

**Hardware**

Focus on XC50 GPU partitition

Each node:

- Processor: 1 x 12-core Intel Xeon E5-2690 v3 @ 2.60GHz (one socket)
- Memory: 64GB RAM
- GPU: 1 x NVidia P100
- CRAY Aries interconnect


**Software**

- Cray MPICH MPI library
- Cray-optimised FFTW
- Cray-libsci provides BLAS & LAPACK
- craype-accel-nvidia60 targets the correct SM architecture to compile for the P100 GPU
  

**Build**

Build instructions Piz Daint's XC50 GPU partition: 

.. include:: cookbook/build/pizdaint/daint-gpu/gromacs-2020.2.rst
   
   
**Run**

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


Keeping in mind best practice to associate one rank (one spatial domain) with each GPU means we run ``mdrun`` with one MPI rank per node since there is only GPU on each node. We choose 12 OpenMP threads (``-ntomp 12``) per MPI rank in order to use all CPU cores, or 24 threads per rank if hyperthreading is enabled. Since there is only one processor (one socket) on Piz Daint GPU compute nodes, we do not expect to suffer a performance penalty by having a single rank spanning across all 12 cores using multithreading. Using default GPU offloaded settings means PME is offloaded to GPU, which can currently only happen on a single rank, corresponding to running ``mdrun`` with ``-npme 1``. 

[Figure showing scaling of benchmark performance for these options, i.e. with single PME rank offloaded to GPU]

Disabling PME GPU offloading by choosing ``-pme cpu`` means PME computations are done on CPU cores. Consider using ``tune_pme`` to determine optimal number of PME ranks.  

[Figure showing scaling of benchmark performance with increasing number of PME ranks starting with 1 (=1 node)] 






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

- Error (abort):

  ::

   Feature not implemented:
   PME GPU does not support PME decomposition

  - Need to pass ``-npme 1`` as an option to ``mdrun``.

- Error (abort):

  ::

   Inconsistency in user input:
   Update task on the GPU was required,
   but the following condition(s) were not satisfied:
   Domain decomposition without GPU halo exchange is not supported.
   With separate PME rank(s), PME must use direct communication.
   

  - Is because of "Update and constraints on a GPU is currently not supported with domain decomposition, free-energy, virtual sites, Ewald surface correction, replica exchange, constraint pulling, orientation restraints and computational electrophysiology."

- Error (abort):

  ::

     128 OpenMP threads were requested. Since the non-bonded force buffer reduction
     is prohibitively slow with more than 64 threads, we do not allow this. Use 64
     or less OpenMP threads.

- Performance warning:

  ::
  
     NOTE: PME load balancing increased the non-bonded workload by more than 50%.
     For better performance, use (more) PME ranks (mdrun -npme),
     or if you are beyond the scaling limit, use fewer total ranks (or nodes).

     - Use ``tune_pme``
       
- Error (abort):

  ::
  
     Update task on the GPU was required, but the following condition(s) were not satisfied:
     Virtual sites are not supported.
     Non-connecting constraints are not supported
     The number of coupled constraints is higher than supported in the CUDA LINCS

- Fatal error:

  ::
     
     There is no domain decomposition for 11 ranks that is compatible with the
     given box and a minimum cell size of 0.79375 nm
     Change the number of ranks or mdrun option -rdd or -dds
     Look in the log file for details on the domain decomposition
     
     
- any other warnings? (e.g. CUDA related such as wrong or suboptimal sm_arch, OpenCL-related, ...?)






------------------------------------------------------------------
GROMACS Reference Benchmarks Performance on PRACE/EuroHPC machines
------------------------------------------------------------------

This section provides a reference set of benchmark simulation
performance results representative of good obtainable performance on
PRACE/EuroHPC machines with GROMACS built and run according to best
practice outlined in this guide.

These results are intended as a convenient reference to help
researchers estimate compute time requirements for their proposed
research in preparation for applying to HPC resource allocation calls.










