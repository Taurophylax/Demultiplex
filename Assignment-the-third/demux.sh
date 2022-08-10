#!/bin/bash
#SBATCH --partition=bgmp       ### Partition (like a queue in PBS)
#SBATCH --time=0-18:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp         ### Account used for job submission
#SBATCH --mem=18G

#python demux.py -i R1.fq.gz R2.fq.gz R3.fq.gz R4.fq.gz indexes.txt 
python demux.py -i "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" "/projects/bgmp/shared/2017_sequencing/indexes.txt"