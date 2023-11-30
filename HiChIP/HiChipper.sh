#!/bin/bash
# This script runs the HiChipper pipeline
# Requires HiCPro.sh to be run first (including running the second step manually)
# Requires a yaml file specifying experimental parameters in a yaml/ folder in the directory executing this script.

folder=$(pwd)  # Saves folder as a variable
mkdir -p PBS
mkdir -p log

cat >${folder}/PBS/HiChipper.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=HiChipper
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/HiChipper_%j.txt -e ${folder}/log/HiChipper_%j.err.txt

#------- End of header -------#

cd ${folder}
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
conda activate hichipper_RW
hichipper --out HiChipper yaml/HiChipper.yaml --make-ucsc
EOF
sbatch ${folder}/PBS/HiChipper.pbs