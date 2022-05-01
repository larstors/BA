#!bin/bash

#$ -cwd
#$ -N crumble
#$ -q softcm.7.day
#$ -l h_vmem=1.75G

./rungrid.pl $SGE_TASK_ID go
