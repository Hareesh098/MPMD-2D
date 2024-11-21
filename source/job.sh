#!/bin/bash
# declare a name for this job to be sample_job
#PBS -N Test
# request the parallel queue for this job
##PBS -q parallel
# request a total of 28 processors for this job (1 nodes and 24 processors per node)
#PBS -l nodes=1:ppn=8
# Request to run for specified hrs:min:secs
#PBS -l walltime=500:00:00
# combine PBS standard output and error files
#PBS -j oe
# mail is sent to you when the job starts and when it terminates or aborts
#PBS -m bea
# specify your email address
#PBS -M charan.harish@gmail.com
#change to the directory where you submitted the job
cd $PBS_O_WORKDIR
PREFIX=G2000-P-L-Test
#include the full path to the name of your MPI program
mpiexec.hydra -n 8 /home/hcharan/MPMD-v2.0/source/main $PREFIX
mkdir ../output/$PREFIX
mv ../output/$PREFIX.* ../output/$PREFIX
exit 0


