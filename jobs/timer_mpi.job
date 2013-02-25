#!/bin/bash

# Request 128 cpus. ncpus must be a multiple of 16.
#PBS -l ncpus=128

# Limit to 20 minutes of walltime.
#PBS -l walltime=20:00           

# Merge stdout and stderr into one output file.
#PBS -j oe

# Run in the batch queue.
#PBS -q batch

# Use the name timer.job.
#PBS -N timer_mpi.job

# Set this to your working directory.
execdir=$HOME/asst3/prog2_wsp_mpi

# Load mpi and python 2.7.
source /usr/share/modules/init/bash
module load openmpi/1.6/gnu
module load python

# Move to my $SCRATCH directory.
cd $SCRATCH

# Copy executable and necessary files to $SCRATCH.
cp $execdir/wsp wsp

# Copy over the scoreboard and quit early if we fail to do so.
cp $execdir/scoreboard_token .
if [ $? -ne 0 ]
then
 exit -1
fi

# The script depends on mkinput.py to generate input.
cp $execdir/mkinput.py .

# The script depends on the existance of an input directory.
mkdir -p input

# The script depends on the existance of a Makefile with a run target.
echo "run: " > Makefile

# Copy over the timer script, modifying a few variables:
#   set num_threads to the actual number of CPUS used
#   set num_iterations to 1, since we assume we have full access to the CPU
cat $execdir/timer.bash | sed 's/num_threads=.*/num_threads=$PBS_NCPUS/' \
                        | sed 's/num_iterations=.*/num_iterations=1/'    \
                        > timer.bash

# Run the script.
bash timer.bash
