#!/bin/bash

# The log fold change will be calculated as log(samp/control), or log(samp) - log(control)

# Script that creates a diffBind comparison, and completes the standard analysis

# This pipeline requires four positional arguments:
# 	$1 - The DiffBind category on which to make the comparison
# 	$2 - First category for comparison (sample group)
# 	$3 - Second category for comparison (control group)
# 	$4 - (Optional) name of the object from diffBind_createObject.sh to read in (default is 'dbObj.rds')

# NOTE: This assumes that there are only two reps. If there are more than two reps, then you must manually change that parameter within the script

# Requires the ATACQC pipeline (will rename in the future)

compare=$1
samp=$2
control=$3
normType=$4
dbObj=$(basename "$5" ".csv")

caps=$(echo "$compare" | tr '[:lower:]' '[:upper:]')

#  All the parameters are not provided
if [ -z "$3" ]; then 
  echo ERROR: NOT ALL PARAMETERS HAVE BEEN SPECIFIED
  echo USAGE:
  echo This pipeline takes in a minimum of three positional argument:
  echo '$1 - The DiffBind category on which to make the comparison'
  echo '$2 - First category for comparison (sample group)'
  echo '$3 - Second category for comparison (control group)'
  echo "\$4 - (Optional) name of the object from diffBind_createObject.sh to read in (default is 'dbObj.rds')"
  exit 1
fi

off=FALSE

if [ "$normType" == Default ]; then
    norm=DBA_NORM_DEFAULT
  elif [ "$normType" == TMM ]; then
    norm=DBA_NORM_TMM
  elif [ "$normType" == RLE ]; then
    norm=DBA_NORM_RLE
  elif [ "$normType" == Loess ]; then
    norm=DBA_NORM_OFFSETS
    off=TRUE
else
  echo ERROR: PROPER NORMALIZATION WAS NOT SPECIFIED
  echo USAGE:
  echo This pipeline takes in three positional arguments:
  echo 	"\$1 - target folder"
  echo 	"\$2 - Normalization Strategy (Default, TMM, RLE, Loess)"
  echo 	"\$3 - (Optional) The name of the generated DiffBind object (default is dbObj)"
  exit 1
fi

# Set up necessary files
mkdir -p PBS
mkdir -p log
subfolder="diffBind/${dbObj}/${samp}_over_${control}_${normType}/"
mkdir -p ${subfolder}

folder=$(pwd)

cat >${folder}/PBS/CorrPlot'.R' <<EOF
# Load libraries
library(DiffBind)
library(dplyr)
library(ggplot2)

rawSamples <- read.csv('$5')
  dbObj <- dba(sampleSheet=rawSamples)
  dbObj <- dba.blacklist(dbObj, blacklist = DBA_BLACKLIST_GRCH38)
  dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel=FALSE)  # Does not run in parallel due to file size
  dbObj <- dba.normalize(dbObj, method = DBA_ALL_METHODS, normalize=${norm}, offsets=${off})

  corr <- plot(dbObj)

  pdf("diffBind/${dbObj}_${normType}_multicorr.pdf")
  corr
  dev.off()

  saveRDS(corr,file="Corr.RDS")
}
EOF

cat >${folder}/PBS/CorrPlot.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=${normType}_${samp}_over_${control}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${samp}_over_${control}_%j.txt -e ${folder}/log/${samp}_over_${control}_%j.err.txt

#------- END OF HEADER -------#
source activate ATACQC

cd ${folder}

Rscript ./PBS/CorrPlot'.R'

EOF

sbatch ${folder}/PBS/CorrPlot.sbatch 