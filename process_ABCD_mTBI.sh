#!/bin/sh
#SBATCH --mail-user=guidoguberman@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

#SBATCH --nodes=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=7-0

export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)

srun nextflow -c ~/scratch/data/singularity.conf run ~/scripts/tractoflow/main.nf --bids ~/scratch/data/ABCD_mTBI_data_complete --output_dir ~/scratch/data/ABCD_mTBI_results_complete_int --wm_seeding false --seeding npv --nbr_seeds 24 --fodf_shells "0 3000" --dti_shells "0 500 1000" -with-singularity ~/scripts/singularity_tractoflow.img -with-report ~/scratch/data/report.html -w ~/scratch/data/ABCD_mTBI_work_complete_int -with-mpi -resume
