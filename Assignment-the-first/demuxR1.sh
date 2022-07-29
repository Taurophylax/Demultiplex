#!/bin/bash
#SBATCH --partition=bgmp       ### Partition (like a queue in PBS)
#SBATCH --time=0-12:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp         ### Account used for job submission
#SBATCH --mem=12G

python part1.py -s 101 -o "plotR1.png" -i "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"