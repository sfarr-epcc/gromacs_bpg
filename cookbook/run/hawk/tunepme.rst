::

   #!/bin/bash -l
   #PBS -N tunepme
   #PBS -l select=16:mpiprocs=32:ompthreads=4
   #PBS -l walltime=01:00:00

   # Change to the directory that the job was submitted from
   cd $PBS_O_WORKDIR

   export OMP_NUM_THREADS=4
   
   export PATH=$HOME/gromacs/2020.2/mpi/AVX2_256/bin:$PATH
   
   gmx_mpi tune_pme -np 512 -mdrun "omplace -nt ${OMP_NUM_THREADS} gmx_mpi mdrun" -s benchRIB.tpr -ntomp $OMP_NUM_THREADS -dlb no 




