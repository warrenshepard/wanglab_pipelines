#!/bin/bash

# USAGE: 
# This pipeline takes in one positional argument:
# 	$1 - target folder

# This script looks at all the bed files located in the user's desired folder,
# and runs HOMER's findMotifs program offshore on the discovery HPC. 

# Requires that HOMER is available in you PATH variable

mkdir -p PBS
mkdir -p log

folder=$(pwd)

# If a name is not provided
if [ -z "$1" ]; then 
  echo ERROR: TARGET FOLDER WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in one positional argument:
  echo 	\$1 - target folder
  exit 1
fi

mkdir -p ${folder}/homer_bgd/$1/

cat >${folder}/PBS/motif_bgd'.pbs' <<EOF
#!/bin/bash -l

#SBATCH --job-name=$1_motifs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=16:00:00
#SBATCH -o ${folder}/log/motifs_%j.txt -e ${folder}/log/motifs_%j.err.txt

#---------------------End of header---------------------#

# Change to job working directory
cd ${folder}

findMotifsGenome.pl $1/intersects/ARID1A_over_ARID1B.broadPeak hg38 homer_bgd/$1/ -size given -bg $1/CR_HCT_WT_ARID1B_broad_peaks.broadPeak

EOF

	sbatch ${folder}/PBS/motif_bgd.pbs
