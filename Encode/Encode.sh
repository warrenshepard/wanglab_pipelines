#!/bin/bash

# Creates the json file for encode atac pipeline
mkdir -p json
mkdir -p Encode

folder=$(pwd)  # Save current folder as a variable

# CHANGE FOR FILE EXTENSION
bioRep=_R1_S
suffix1=_R1_001.fastq.gz

# Loop over all files, create the alignment script for each
for file in fastq/*${bioRep}*${suffix1}; do
	base=$(basename "$file" ${suffix1})
	smallBase=${base%${bioRep}*}
	echo ${smallBase}
	Num1="0"
	for rep in fastq/*${smallBase}*${suffix1}; do
	largeBase=$(basename "$rep" ${suffix1})
	number=$(echo "$largeBase" | sed 's/.*_S\([0-9]\+\).*/\1/')
	if [ $Num1 == "0" ]; then
	Num1=$number
	else
	Num2=$number
	fi
	done
	file1="../fastq/${smallBase}_R1_S${Num1}"
	file2="../fastq/${smallBase}_R2_S${Num2}"

	INPUT_JSON=${folder}/json/$smallBase'.json'
	cat >${INPUT_JSON} <<EOF
	{
	"atac.title" : "${smallBase}",
    "atac.description" : "Encode standard output for ${smallBase}.",

    "atac.pipeline_type" : "atac",
    "atac.align_only" : false,
    "atac.true_rep_only" : false,

    "atac.genome_tsv" : "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv",

    "atac.paired_end" : true,

    "atac.fastqs_rep1_R1" : [ "${file1}${suffix1}" ],
    "atac.fastqs_rep1_R2" : [ "${file1}${suffix1/R1/R2}" ],
    "atac.fastqs_rep2_R1" : [ "${file2}${suffix1}" ],
    "atac.fastqs_rep2_R2" : [ "${file2}${suffix1/R1/R2}" ],

    "atac.auto_detect_adapter" : true,

    "atac.multimapping" : 4
	}
EOF

cd Encode
caper hpc submit /dartfs-hpc/rc/lab/W/WangX/ADCM/atac-seq-pipeline/atac.wdl -i "${INPUT_JSON}" --singularity --leader-job-name ${smallBase}
cd ${folder}
done
