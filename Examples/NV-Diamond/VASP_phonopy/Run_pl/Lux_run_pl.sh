#!/bin/bash
#SBATCH --job-name=nscf
#SBATCH --output=qe.%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00
#SBATCH --partition=gpuq
#SBATCH --dependency=afterany:269811
#SBATCH --account=gpuq
##SBATCH --exclude=krnode15,krnode11
#SBATCH --mail-type=ALL
#SBATCH --mail-user=szhan213@ucsc.edu
#################### Kairay ################################################################
# module add intel/17.0.5.239 impi/2017
# export OMP_NUM_THREADS=1
# NORES=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))
# MPICMD="mpirun -genv I_MPI_FABRICS shm:ofa -n $SLURM_NTASKS"
#PWDIR="/export/data/share/wufeng/programs-intel2017.5/qe-6.1-scal/bin"
#PWDIR="/export/data/share/jxu/qe-codes/qe-6.4.1/bin"
############################################################################################


#################################### Lux ############################################################################################################
 hostname
 echo "Running program on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS total tasks, with each node getting $SLURM_NTASKS_PER_NODE running on cores."
 module load intel/impi
 MPICMD="mpirun -n $SLURM_NTASKS"
#PWDIR=/data/users/jxu153/codes/qe/qe-6.1.0/bin
#PWDIR=/data/groups/ping/kli103/programs/qe-6.6/bin
#PWDIR=/data/groups/ping/kli103/programs/qe-7.2/bin
#PWDIR=/data/groups/ping/jxu153/codes/qe/qe-7.1/bin
PWDIR=/data/users/jxu153/codes/qe/qe-6.1.0/bin

#####################################################################################################################################################

ENV='/home/szhan213/miniconda3/envs/west/bin'

echo "Start:"; date
#PL_CODE='/export/data/share/szhan213/Programs/PL-code/code' #This is the location pl code
PL_CODE="/data/groups/ping/szhan213/Programs/PL-lineshape-code_dev/code"
$ENV/python3 -u $PL_CODE/main.py | tee pl.out

echo "Done:"; date
