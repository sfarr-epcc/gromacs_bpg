.. code:: bash

    #!/bin/bash
    #SBATCH --account=<ACCOUNT>
    #SBATCH --nodes=1
    #SBATCH --ntasks=8           # number of simulations
    #SBATCH --hint=nomultithread
    #SBATCH --ntasks-per-node=8  # number of simulations
    #SBATCH --cpus-per-task=6    # threads per simulation
    #SBATCH --time=01:00:00
    #SBATCH --partition=booster
    #SBATCH --gres=gpu:4         # use all 4 GPUs
    #SBATCH --cuda-mps           # turn on CUDA MPS to enable multiple processes per GPU




    # load the required modules
    module load GCC/11.2.0
    module load CUDA/11.5

    module list

    # Allow GROMACS to see all 4 GPUs
    export CUDA_VISIBLE_DEVICES=0,1,2,3


    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}


    # make sure we are using a non-mpi version of gromacs
    source <user install location>/gromacs-2022_NOMPI/bin/GMXRC

    # input file
    TPR=input.tpr


    EXE="gmx mdrun"  # Non MPI mdrun

    # Set thread mpi to 1 with -ntmpi option
    SETTINGS="-s ${TPR} -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK} -nb gpu -pme gpu -bonded gpu"


    # Options below for 4,8,16,24, or 48 simulations per node.
    # Uncomment the desired "folder" option and the corresponding "gputasks" option,
    # make sure the rest are commented out and ensure the --ntasks,
    # --ntasks-per-node, and --cpus-per-task are set to the corresponding values.

    
    # Currently the 8 simulations option is being used.



    # 4 simulations: (make sure ntasks and ntasks-per-node=4 and cpus-per-task=12)
    #folders=(1 2 3 4)
    #gputasks=(11 00 33 22)

    # 8 simulations: (make sure ntasks and ntasks-per-node=8 and cpus-per-task=6)
    folders=(1 2 3 4 5 6 7 8)
    gputasks=(11 11 00 00 33 33 22 22)

    # 16 simulations: (make sure ntasks and ntasks-per-node=16 and cpus-per-task=3)
    #folders=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
    #gputasks=(11 11 11 11 00 00 00 00 33 33 33 33 22 22 22 22)

    # 24 simulations: (make sure ntasks and ntasks-per-node=24 and cpus-per-task=2)
    #folders=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
    #gputasks=(11 11 11 11 11 11 00 00 00 00 00 00 33 33 33 33 33 33 22 22 22 22 22 22)

    #48 simulations: (make sure ntasks and ntasks-per-node=48 and cpus-per-task=1)
    #folders=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47)
    #gputasks=(11 11 11 11 11 11 11 11 11 11 11 11 00 00 00 00 00 00 00 00 00 00 00 00  33 33 33 33 33 33 33 33 33 33 33 33 22 22 22 22 22 22 22 22 22 22 22 22 )





    # make cpu bindings to ensure cores are not shared between the simulations
    # this works for any of the simulation number options, does not need to be changed
    bindings=($( for (( i=0 ; i<48 ; i+=${SLURM_CPUS_PER_TASK} )); do echo $i-$((${i} + ${SLURM_CPUS_PER_TASK} - 1 )); done ))

    # This works in conjunction with the gputasks settings which map the process to GPU close to the
    # specified GPU.


    # lauch all simulations in the specified folders in parallel
    for i in ${!folders[@]}
    do
        cd ${folders[$i]}
        pwd
        echo ${bindings[$i]}

        # use numactl to bind the simulations to specific cores.
        numactl --physcpubind=${bindings[$i]} ${EXE} ${SETTINGS} -gputasks ${gputasks[$i]} > out.txt 2>&1 &
        cd ..
    done

    wait

