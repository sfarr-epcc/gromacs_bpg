::

   #!/bin/bash -l
   #PBS -N benchMEM
   #PBS -l select=2:mpiprocs=128:ompthreads=2
   #PBS -l walltime=00:10:00

   # Change to the directory that the job was submitted from
   cd $PBS_O_WORKDIR
   
   export OMP_NUM_THREADS=2

   mpirun -np 256 omplace -nt 2 -$HOME/gromacs/2020.2/mpi/AVX2_256/bin/gmx_mpi mdrun -s ../benchMEM.tpr -ntomp 2

