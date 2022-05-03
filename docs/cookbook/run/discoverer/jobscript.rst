::
    
  #!/bin/bash
  #SBATCH --partition=cn
  #SBATCH --job-name=gmx_mdrun
  #SBATCH --time=24:00:00
  #SBATCH --nodes           2   # 2 Nodes
  #SBATCH --ntasks-per-node 128 # 1 MPI task per physical core
  #SBATCH --cpus-per-task   2   # 2 Open MPI threads per physical core (SMT in use)
  
  module purge

  # Load the most recent OpenMPI+GCC build
  module load gromacs/2021/latest-intel-nogpu-openmpi-gcc

  
  export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
  
  
  # standard run command for mdrun
  mpirun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -s intput.tpr
  
  
  # command for consistent benchmarking performance
  #mpirun gmx_mpi mdrun -ntomp ${SLURM_CPUS_PER_TASK} -v -s input.tpr -resethway -dlb yes -notunepme -noconfout