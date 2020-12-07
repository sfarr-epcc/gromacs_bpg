::

   #!/bin/bash -l
   #PBS -N benchRIB
   #PBS -l select=8:mpiprocs=128:ompthreads=2
   #PBS -l walltime=01:00:00             
   
   # Change to the directory that the job was submitted from
   cd $PBS_O_WORKDIR

   export OMP_NUM_THREADS=1

   mpirun -np 1024 omplace -nt $OMP_NUM_THREADS -ht compact gmx_mpi mdrun -s benchRIB.tpr -ntomp $OMP_NUM_THREADS






