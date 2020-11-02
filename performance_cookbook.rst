====================
Performance Cookbook
====================

This part of the GROMACS best practice guide assumes your simulations are prepared appropriately and focuses on providing guidance related to running GROMACS simulations, i.e. executing ``mdrun``.


---------------------------------------------------------------------------
How to run ``mdrun`` and get good performance
---------------------------------------------------------------------------
This section provides a step-by-step strategy for determining how best to run your simulation with ``mdrun`` so as to make good use of available hardware and/or to obtain results in the shortest time possible, be it on a laptop, a multi-GPU desktop workstation, a departmental cluster, or a large supercomputer.


*provide a structured step-by-step approach for users to follow by identifying ordered priorities based on information in the "getting good performance from mdrun" manual page. Link out the manual page for more information, as discussed the BPG can act in part as a guide (a bridge) to the details in the manual*


Single compute node, desktop workstation, or laptop
---------------------------------------------------


++++++++
CPU only
++++++++
- pinning 



+++++++++
CPU + GPU
+++++++++



Multiple networked compute nodes in a cluster or supercomputer
--------------------------------------------------------------

++++++++
CPU only
++++++++


+++++++++
CPU + GPU
+++++++++







----------------------------------------------------------
Acting on performance-related warnings found in ``md.log``
----------------------------------------------------------
This section provides guidance on identifying, understanding, and acting on performance-related warnings and suggestions issued by ``mdrun`` that you may encounter in ``md.log``. These may suggest more efficient ways to launch ``mdrun``, or spot if the GROMACS installation you are using has not been built to give good performance on the hardware and suggest how to improve this. 


*Provide practical guidance on how to e build is Link out to Installation Guide in manual (again, acting as bridge / guide to the manual)*


- SIMD warning in ``md.log``, e.g.: ::

   Highest SIMD level requested by all nodes in run: AVX_512
   SIMD instructions selected at compile time:       AVX2_256
   This program was compiled for different hardware than you are running on,
   which could influence performance. This build might have been configured on a
   login node with only a single AVX-512 FMA unit (in which case AVX2 is faster),
   while the node you are running on has dual AVX-512 FMA units.

  - proviso regarding narrower SIMD width (AVX2_256) potentially better (than AVX_512) for Intel Skylake & Cascade Lake
   
- any other warnings? (e.g. CUDA related such as wrong or suboptimal sm_arch, OpenCL-related, ...?)



  
---------------------------------------------
Building GROMACS on PRACE/EuroHPC machines
---------------------------------------------
This section provides concrete recipes showing how best to build GROMACS for good performance on some of the largest EU-based supercomputers available to EU researchers through PRACE (and in future also EuroHPC).


Piz Daint (CSCS, Switzerland)
------------------------------------------------

Build instructions Piz Daint's XC50 GPU partition: 

.. include:: build/pizdaint/daint-gpu/gromacs-2020.2.rst
   






HAWK (HLRS, Germany)
--------------------



Marconi M100 (CINECA, Italy)
----------------------------




---------------------------------------------
Running GROMACS on PRACE/EuroHPC machines 
---------------------------------------------
This section provides guidance on how best to run GROMACS on a number of specific PRACE/EuroHPC machines.

eneral guidance for each machine is complemented by an analysis illustrating the effect of key runtime execution options on ``mdrun`` performance for a range of benchmark simulations, comparison with which should help you determine how best to run your own simulations on these machines. 

In each case, assume we are using an optimal GROMACS build as described in this guide. 



Piz Daint (CSCS, Switzerland)
-----------------------------
XC50 GPU partition


HAWK (HLRS, Germany)
--------------------



Marconi M100 (CINECA, Italy)
----------------------------




---------------------------------------------------------------------
GROMACS Reference Benchmarks Performance on PRACE/EuroHPC machines
---------------------------------------------------------------------
This section provides a reference set of benchmark simulation performance results representative of good obtainable performance on PRACE/EuroHPC machines with GROMACS built and run according to best practice outlined in this guide.

These results are intended as a convenient quick reference to help researchers estimate compute time requirements for their proposed research in preparation for applying to HPC resource allocation calls.










