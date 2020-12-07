====================
Performance Cookbook
====================

The performance cookbook part of the GROMACS best practice guide
assumes your simulations are prepared appropriately and provides
concrete guidance on how best to run GROMACS simulations, i.e. execute
``mdrun``, so as to make good use of available hardware and obtain
results in the shortest time possible, be it on a laptop, a multi-GPU
desktop workstation, a departmental cluster, and especially on large
supercomputers. This complements and provides a bridge into navigating
the detailed information provided in the "getting good performance
from mdrun" page in the GROMACS manual:

http://manual.gromacs.org/current/user-guide/mdrun-performance.html

GROMACS can generally be launched without specifying anything other
than the essential simulation parameters, as it has built-in
heuristics that enable it to detect the underlying hardware and use
accumulated insights about good performance embedded in the code to
make usually reasonable choices given the available number and types
of CPU cores and/or GPUs. By default GROMACS also adapts dynamically
during execution to improve performance. However for any given
simulation it is possible that better than default choices
exist. Understanding how to control these by explicitly specifying
parallel execution options and how best to approach obtaining optimal
use of available hardware can make a significant difference to
throughput and hence scientific results achieved over a given
timespan, as well as to financial (cost) and environmental (energy
usage) efficiency. For using high-performance computing resources
there is additionally clear value in knowing what scale of resources
(number of cores / nodes / GPUs) are efficient to use, or how to go
about finding this out. 

As well as general guidance applicable to whatever machine you may be
running GROMACS on, the performance cookbook provides concrete
examples showing how to obtain good performance on a number of
specific PRACE (and in future EuroHPC) machines, both as an
illustration of the application of general best practice process for
obtaining good performance, and to promote efficient usage of the
named machines. The cookbook also provides a reference set of (near)
optimal benchmark performance results obtained on these machines using
best practice in order to aid estimation of required compute time
allocations for researchers requesting such time. 


-------------------------------------------------------------------
General guidance for running ``mdrun`` and getting good performance
-------------------------------------------------------------------

Single node
-----------

When using GROMACS on a single node, that is to say a machine such as
a laptop, desktop workstation, or server, where all processor cores or
GPUs have access to a single shared memory, we typically run the
thread-mpi version of mdrun (``gmx mdrun``).

++++++++
CPU only
++++++++

On a single-node CPU-only machine one can experiment with different
combinations of numbers of thread-MPI ranks (controlled by varying
``-ntmpi N``) and numbers of OpenMP threads per rank (controlled by
varying ``-ntomp M`` and setting the environment variable
``OMP_NUM_THREADS=M``), such that ``M x N`` is equal to the total
number of CPU cores available on the node. Many processors support
simultaneous multithreading (SMT), known as hyperthreading for Intel
processors, whereby each physical core can run multiple threads or
processes. Enabling multithreading may boost performance. If enabled,
``-ntmpi N`` and ``-ntomp M`` should be chosen such that ``M x N``
equals the number of logical cores identified by the operating
systems, which is equal to the number of physical cores times by the
number of multithreads per physical core.

For details, see `Running mdrun within a single node
<http://manual.gromacs.org/current/user-guide/mdrun-performance.html#running-mdrun-within-a-single-node>`_
and `Process(-or) level parallelization via OpenMP <http://manual.gromacs.org/current/user-guide/mdrun-performance.html#process-or-level-parallelization-via-openmp>`_



+++++++++
CPU + GPU
+++++++++

- Default behaviour one thread-mpi rank per GPU (default) likely optimal

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
    

For details, see `Running mdrun within a single node
<http://manual.gromacs.org/current/user-guide/mdrun-performance.html#running-mdrun-within-a-single-node>`_,
`Node level parallelization via GPU offloading and thread-MPI
<http://manual.gromacs.org/current/user-guide/mdrun-performance.html#node-level-parallelization-via-gpu-offloading-and-thread-mpi>`_,
and `Running mdrun with GPUs <http://manual.gromacs.org/current/user-guide/mdrun-performance.html#running-mdrun-with-gpus>`_
  

Multiple networked compute nodes in a cluster or supercomputer
--------------------------------------------------------------

In order to run GROMACS on a multi-node distributed-memory machine
such as a supercomputer we need to run the MPI-enabled version,
``gmx_mpi mdrun`` or simply ``mdrun_mpi``, rather than ``gmx mdrun``
or simply ``mdrun``. We must launch In addition and launch with parallel application
launcher (``mpirun``, ``mpiexec``, ``srun``, or ``aprun``) which
determines number of ranks (not ``-ntmpi``)

::

   mpirun -np N gmx_mpi mdrun <options> 

For these larger runs, GROMACS can benefit from having MPI ranks dedicated
to PME ranks and spawns a certain fraction of these based on internal settings.

You should examine the md.log file produced and check if there is any
note suggesting to use a larger or smaller number of PME ranks.

A performance-optimal number of PME ranks for a given total number of
ranks can be determined in a systematic way using the ``gmx tune_pme``
tool as described in `gmx tune_pme <http://manual.gromacs.org/current/onlinehelp/gmx-tune_pme.html#gmx-tune-pme>`_


 
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

This section provides guidance and concrete recipes showing how to
build and run GROMACS for good performance on some of the largest
EU-based supercomputers available to EU researchers through PRACE (and
in future also EuroHPC).

General guidance for each machine is complemented by an analysis
illustrating the effect of key runtime execution choices on ``mdrun``
performance for a range of benchmark simulations. Comparison between
these and your own simulations should help you determine how best to
run your own simulations on these machines and how to obtain good
performance. 


Benchmarks
----------

A brief description is provided below of the benchmarks used to
illustrate how to obtain good performance on PRACE/EuroHPC machines
for a range of system sizes and types.

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

- Processors: 2 x 64-core AMD EPYC 7742 @2.25
- Memory: 256GB RAM
- Interconnect: InfiniBand HDR200

**Software**

Relevant software stack on system (available to all users via environment modules):

- HPE MPT and OpenMPI MPI libraries
- FFTW (Zen2 architecture-specific build)
  
**Build**

A multinode-capable MPI-enabled version of GROMACS with good performance on HAWK can be built as follows:

.. include:: build/hawk/gromacs-2020.2.rst


**Run**

The example job script below shows how to run GROMACS on HAWK for 1 hour on 8 nodes with 128 MPI ranks per node and 2 OpenMP threads per rank. Each physical core on the two 64-core AMD EPYC processors on HAWK supports two simultaneous multithreads (SMTs) - two logical cores - providing a total of 256 usable logical cores per node. The example launches a total of 256 threads across 128 ranks, implying that we intend to make use of all logical cores. The ``omplace -ht compact`` option should be used when running GROMACS using simultaneous multithreading as it ensures similarly numbered MPI ranks as well as OpenMP threads belonging to the same MPI rank are executed as close together in the processor's and indeed the node's memory hierarchy as possible - in the example, both OpenMP threads run on the same physical core. Not using the compact omplace option was found to be more likely to lead to lower performance when using SMT. 

.. include:: run/hawk/jobscript.rst

	     
+++++++++++++++++++++++++++++++++++++++++++++++++++++++
MPI vs OpenMP hybrid parallelism and SMT multithreading
+++++++++++++++++++++++++++++++++++++++++++++++++++++++

In order to better understand how GROMACS utilises the hardware and
make appropriate choices we can investigate the effect on benchmark
performance of the choice of hybrid MPI \& OpenMP execution and use of
SMT for a given number of HAWK nodes. Doing this systematically is
facilitated in the first instance by disabling dynamic load balancing
(``-dlb no``) and PME tuning (``-tunepme no``), which also allows us
to illustrate the strength of load imbalance in different MPI/OpenMP
execution scenarios. The figures below show benchmark performance on
HAWK for different combinations of MPI ranks and OpenMP threads and
with and without use of multithreading (SMT). Eror bars (offset upward
from measurements) indicate the hypothetical performance if the
combined load imbalance between spatial domains and between PP and PME
ranks reported by GROMACS were absent.

It is clear there is a very significant effect on performance of the
choice of MPI/OpenMP hybrid decomposition. As a general rule on HAWK
it is likely that performance is best for the smallest number of
OpenMP threads per MPI rank, though it may be possible to achieve
better scaling on larger number of nodes with a modest number of
threads per rank. Using threading may allow use of a convenient number
of total nodes by avoiding the fatal error of choosing a number of PP
ranks that have too large a prime factor as largest divisor, causing
GROMACS to abort with a fatal error. 



.. list-table:: 
   :align: center
	   
   * - .. figure:: results/hawk/20k_HBS_1smt_resethway_dlbNO_tunepmeNO-nsperday.svg

          20k_HBS, SMT off (1 thread per core)
     
     - .. figure:: results/hawk/20k_HBS_2smt_resethway_dlbNO_tunepmeNO-nsperday.svg

          20k_HBS, SMT on (2 threads per core)
	  
   * - .. figure:: results/hawk/benchMEM_1smt_resethway_dlbNO_tunepmeNO-nsperday.svg

          benchMEM, SMT off (1 thread per core)
	   
     - .. figure:: results/hawk/benchMEM_2smt_resethway_dlbNO_tunepmeNO-nsperday.svg

          benchMEM, SMT on (2 threads per core)

   * - .. figure:: results/hawk/465k_HBS_1smt_dlbNO_tunepmeNO-nsperday.svg

          465k_HBS, SMT off (1 thread per core)
	   
     - .. figure:: results/hawk/465k_HBS_2smt_dlbNO_tunepmeNO-nsperday.svg

          465k_HBS, SMT on (2 threads per core)

   * - .. figure:: results/hawk/benchRIB_1smt_dlbNO_tunepmeNO-nsperday.svg

          benchRIB, SMT off (1 thread per core)

     - .. figure:: results/hawk/benchRIB_2smt_dlbNO_tunepmeNO-nsperday.svg

          benchRIB, SMT on (2 threads per core)

   * - .. figure:: results/hawk/benchPEP_1smt_dlbNO_tunepmeNO-nsperday.svg

          benchPEP, SMT off (1 thread per core)

     - .. figure:: results/hawk/benchPEP_2smt_dlbNO_tunepmeNO-nsperday.svg

          benchPEP, SMT on (2 threads per core)
	   






**Tuning the number of PME ranks**

Unless disabled with ``-tunepme no`` GROMACS by default attempts to
reduce the load imbalance between PP and PME ranks. However it does
not change the number of dedicated PP and PME ranks. These are fixed,
and unless chosen by specifying the number of PME ranks using
``-npme`` are determined according to various heuristics based on what
is expected be optimal based on the underlying algorithm and what has
often been found to give good performance.

As `mentioned in the manual
<http://manual.gromacs.org/current/user-guide/mdrun-performance.html#running-mdrun-on-more-than-one-node>`_
the domain decomposition (DD) load balancing functionality, which is
enabled by default and disabled with ``-dlb no``, is important for
achieving good performance for spatially heterogeneous
systems. However PME and DD load balancing can interfere with each
other. To improve on the unbalanced performance shown in above
benchmark figures therefore, a systematic approach can be taken by
separately tuning the number of PME for a given total number of ranks
using the ``gmx tune_pme`` tool as described in `the section on
tune_pme
<http://manual.gromacs.org/current/onlinehelp/gmx-tune_pme.html#gmx-tune-pme>`_
in the manual.

On HAWK, we could tune the number of PME ranks to improve the
performance of, for example, the benchRIB benchmark running on 16
nodes with 32 MPI ranks per node (512 ranks in total) and 4 OpenMP
threads per rank, which we saw in the above results is already a good
choice considering unbalanced performance and which has scope to
improve significantly through reduction of observed load imbalance.
The following script will allow us to run ``tune_pme`` on HAWK to do
this:

.. include:: run/hawk/tunepme.rst





  

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

Relevant software stack on system (available to all users via environment modules):

- Cray MPICH MPI library
- Cray-optimised FFTW
- Cray-libsci provides BLAS & LAPACK
- craype-accel-nvidia60 targets the correct SM architecture to compile for the P100 GPU
  

**Build**

Build instructions Piz Daint's XC50 GPU partition: 

.. include:: build/pizdaint/daint-gpu/gromacs-2020.2.rst
   
   
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


**Offload scenarios**
   
Keeping in mind best practice to associate one rank (one spatial domain) with each GPU means we run ``mdrun`` with one MPI rank per node since there is only GPU on each node. We choose 12 OpenMP threads (``-ntomp 12``) per MPI rank in order to use all CPU cores, or 24 threads per rank if hyperthreading is enabled. Since there is only one processor (one socket) on Piz Daint GPU compute nodes, we do not expect to suffer a performance penalty by having a single rank spanning across all 12 cores using multithreading. Using default GPU offloaded settings means PME is offloaded to GPU, which can currently only happen on a single rank, corresponding to running ``mdrun`` with ``-npme 1``. 

[Figure showing scaling of benchmark performance for these options, i.e. with single PME rank offloaded to GPU]

Disabling PME GPU offloading by choosing ``-pme cpu`` means PME computations are done on CPU cores. Consider using ``tune_pme`` to determine optimal number of PME ranks.  



.. list-table:: 
   :align: center
	   
   * - .. figure:: results/pizdaint/20k_HBS_1smt_dlbNO_tunepmeNO-nsperday.svg

          20k_HBS, SMT off (1 thread per core)
     
     - .. figure:: results/pizdaint/20k_HBS_2smt_dlbNO_tunepmeNO-nsperday.svg

          20k_HBS, SMT on (2 threads per core)
	  
   * - .. figure:: results/pizdaint/benchMEM_1smt_dlbNO_tunepmeNO-nsperday.svg

          benchMEM, SMT off (1 thread per core)
	   
     - .. figure:: results/pizdaint/benchMEM_2smt_dlbNO_tunepmeNO-nsperday.svg

          benchMEM, SMT on (2 threads per core)

   * - .. figure:: results/pizdaint/465k_HBS_1smt_dlbNO_tunepmeNO-nsperday.svg

          465k_HBS, SMT off (1 thread per core)
	   
     - .. figure:: results/pizdaint/465k_HBS_2smt_dlbNO_tunepmeNO-nsperday.svg

          465k_HBS, SMT on (2 threads per core)

   * - .. figure:: results/pizdaint/benchRIB_1smt_dlbNO_tunepmeNO-nsperday.svg

          benchRIB, SMT off (1 thread per core)

     - .. figure:: results/pizdaint/benchRIB_2smt_dlbNO_tunepmeNO-nsperday.svg

          benchRIB, SMT on (2 threads per core)

   * - .. figure:: results/pizdaint/benchPEP_1smt_dlbNO_tunepmeNO-nsperday.svg

          benchPEP, SMT off (1 thread per core)

     - .. figure:: results/pizdaint/benchPEP_2smt_dlbNO_tunepmeNO-nsperday.svg

          benchPEP, SMT on (2 threads per core)
	   






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
     

-  Fatal error:

   ::

      The number of ranks selected for particle-particle work (383) contains a large
      prime factor 383. In most cases this will lead to bad performance. Choose a
      number with smaller prime factors or set the decomposition (option -dd)
      manually.






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










