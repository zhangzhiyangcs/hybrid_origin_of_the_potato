#!/bin/bash
#SBATCH -J wrf
#SBATCH --comment=WRF
#SBATCH -n 38
#SBATCH -N 1
#SBATCH -p xhacnormalb
snakemake -k  -s Snakefile --latency-wait 100  -j 38  --stats snakejob.stats >&2 2>> snakejob.log
