#!/bin/sh
#SBATCH --mail-user=guidoguberman@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

#SBATCH --nodes=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=05:00:00

export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)


srun nextflow -c ~/scratch/data/singularity.conf run ~/scripts/connectoflow/main.nf --input ~/scratch/data/connectoflow_tests --template ~/scratch/data/templates/MNI/MNI152_T1_1mm_brain.nii --apply_t1_labels_transfo false --labels_list ~/scratch/data/templates/DKT/dkt.txt -w ~/scratch/data/connectoflow_tests_work --output_dir ~/scratch/data/connectoflow_tests_results --run_commit false --run_afd_rd true --length_weighting true --use_similarity_metric false -with-singularity ~/scripts/connectoflow/connectoflow15102020.img -with-mpi -resume

#srun nextflow -c ~/scratch/data/singularity.conf run ~/scripts/tractoflow/main.nf --bids ~/scratch/data/ABCD_TP3_data --output_dir ~/scratch/data/ABCD_TP3_results --wm_seeding true --seeding npv --nbr_seeds 12 --fodf_shells "0 3000" --dti_shells "0 500 1000" --b0_thr_extract_b0 "0.1" -with-singularity ~/scripts/singularity_tractoflow.img -with-report ~/scratch/data/report.html -w ~/scratch/data/ABCD_TP3_work -with-mpi -resume
