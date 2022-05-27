::
    
    #!/bin/bash
    #SBATCH --account=<ACCOUNT>
    #SBATCH --nodes=1
    #SBATCH --ntasks=4
    #SBATCH --hint=nomultithread # turn off SMT
    #SBATCH --ntasks-per-node=4  # 1 MPI task per GPU
    #SBATCH --cpus-per-task=12   
    #SBATCH --time=01:00:00
    #SBATCH --partition=booster
    #SBATCH --gres=gpu:4        # Use all 4 GPUs


    # load the required modules
    module load GCC/11.2.0
    module load CUDA/11.5
    module load OpenMPI/4.1.2
    module load FFTW/3.3.10


    # Allow GROMACS to see all 4 GPUs
    export CUDA_VISIBLE_DEVICES=0,1,2,3

    # Enable direct GPU to GPU communications
    export GMX_ENABLE_DIRECT_GPU_COMM=true

    # activate user install of GROMACS 2022
    source <user install location>/gromacs-2022/bin/GMXRC

    srun gmx_mpi mdrun -s INPUT.tpr -ntomp ${SLURM_CPUS_PER_TASK} -nb gpu -pme gpu -bonded gpu -npme 1