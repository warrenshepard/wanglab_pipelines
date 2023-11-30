#!/bin/bash
# This script runs the diffloop pipeline

folder=$(pwd)  # Saves folder as a variable
mkdir -p PBS
mkdir -p log

cat >${folder}/PBS/Diffloop.pbs <<EOF
#!/bin/bash
#SBATCH --job-name=Diffloop
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/Diffloop_%j.txt -e ${folder}/log/Diffloop_%j.err.txt

#------- End of header -------#

cd ${folder}
source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
conda activate hichipper_RW
Rscript /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/HiChIP/diffLoop.R

EOF
sbatch ${folder}/PBS/Diffloop.pbs