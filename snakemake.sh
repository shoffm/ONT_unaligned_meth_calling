#! /bin/bash
# this file is snakemake.sh
module load snakemake samtools || exit 1
module load guppy/5.0.7


snakemake -s snakefile_g --latency-wait 60 --cores 8 -j 3354 --local-cores 4 --cluster "sbatch --ntasks 1 --mem 20g --cpus-per-task 1 --time 70:00:00"

# submission command
# sbatch --cpus-per-task=2 --mem=8g snakemake.sh
