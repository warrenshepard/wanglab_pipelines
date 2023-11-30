#!/bin/bash
# This script runs the HiC-Pro pipeline (full)
# Requires that all fastq files are in individual folders
# Requires a config file called config.txt in the directory the script is run in
# NOTE!!!! This will not run the second step automatically!! 
# following completion of this step, go to terminal and change to HiCPro_Output and run: sbatch HiCPro_step2_HiC-Pro.sh


folder=$(pwd)  # Saves folder as a variable
mkdir -p PBS
mkdir -p log

cat >${folder}/PBS/HiCPro.pbs <<EOF
#!/bin/bash

#SBATCH --job-name=HiC-Pro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH -o ${folder}/log/HiCPro_%j.txt -e ${folder}/log/HiCPro_%j.err.txt

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
cd ${folder}
/dartfs-hpc/rc/lab/W/WangX/Hi-ChIP_master/HiC-Pro-master/bin/HiC-Pro \
	-i ${folder}/fastq \
	-o ${folder}/HiCPro_Output \
	-c ${folder}/config.txt -p

cd HiCPro_Output
sbatch HiCPro_step1_HiC-Pro.sh

EOF
cd ${folder}
sbatch PBS/HiCPro.pbs
