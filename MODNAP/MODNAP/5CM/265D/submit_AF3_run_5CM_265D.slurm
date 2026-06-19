#!/bin/bash

#SBATCH --job-name=af3
#SBATCH --account=gts-amcshan3-paid
#SBATCH --time=72:00:00
#SBATCH -o slurm_run_%j_af3_output.txt
#SBATCH -e slurm_run_%j_af3_log.txt
#SBATCH --mem=80GB
#SBATCH -q inferno
#SBATCH --gres=gpu:A100:1

cd $SLURM_SUBMIT_DIR

input='./af_input'
output='./af_output'

echo $input
echo $output

apptainer exec --bind "$input":/root/af_input --bind "$output":/root/af_output --bind /storage/cedar/cedar0/cedarp-amcshan3-0/alphafold3_param:/root/models --bind /storage/cedar/cedar0/cedarp-amcshan3-0/alphafold3_database:/root/public_databases --nv /storage/cedar/cedar0/cedarp-amcshan3-0/alphafold3/image.sif python /storage/cedar/cedar0/cedarp-amcshan3-0/alphafold3/run_alphafold.py --json_path=/root/af_input/5CM_265D.json --model_dir=/root/models --db_dir=/root/public_databases --output_dir=/root/af_output
