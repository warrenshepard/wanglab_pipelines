#!/bin/bash

# This script takes two positional arguments
#   $1 - (REQUIRED) your desired association rule OPTIONS: basalPlusExt | twoClosest | oneClosest
#   $2 - (REQUIRED) your desired mode OPTIONS: default | diffBind
#   $2 - (REQUIRED) a metafile specifying the .bed inputs you'd like to use

# your metafile should be formatted as defined below. Examples/templates will be located in Warren/GREAT/templates

# If using default mode, each line of your metafile should contain the path to a .bed file
# if using diffBind mode, each line should contain the paths to the folder containing the upregulated, downregulated, and unchanged .bed files

# If you have a lot of input, expect this script to take a while... GREAT limits calls to once every 60 seconds

# help command
if [ "$1" == "--help" ]; then
    echo "$(tput setaf 3)USAGE:"
    echo "This script takes three positional arguments:"
    echo "  \$1 (REQUIRED) your desired association rule OPTIONS: basalPlusExt | twoClosest | oneClosest"
    echo "  \$2 (REQUIRED) your desired mode OPTIONS: default | diffBind"
    echo "  \$3 (REQUIRED) a metafile specifying the .bed inputs you'd like to use"
    echo "examples/templates for metafiles are located in Warren/GREAT/templates$(tput sgr0)"
    exit 1
fi

# error message
if [ -z "$3" ]; then
    echo "$(tput setaf 1)ERROR: One or more required positional aguments not supplied"
    echo "This script takes three positional arguments: "
    echo "  \$1 (REQUIRED) your desired association rule OPTIONS: basalPlusExt | twoClosest | oneClosest"
    echo "  \$2 (REQUIRED) your desired mode OPTIONS: default | diffBind"
    echo "  \$3 (REQUIRED) a metafile specifying the .bed inputs you'd like to use"
    echo "examples/templates for metafiles are located in Warren/GREAT/templates$(tput sgr0)"
    exit 1
fi

# checking/assigning parameters
if [ "$1" == "basalPlusExt" ] || [ "$1" == "twoClosest" ] || [ "$1" == "oneClosest" ]; then
    ass_rule=$1
else
    echo "$(tput setaf 1)ERORR: That is not a valid association rule. Please check your spelling$(tput sgr0)"
fi

if [ "$2" == "default" ] || [ "$2" == "diffBind" ]; then
    mode=$2
else
    echo "$(tput setaf 1)ERORR: That is not a valid association rule. Please check your spelling$(tput sgr0)"
fi

meta=$3
folder=$(pwd)

# setup directories
mkdir -p log
mkdir -p PBS
mkdir -p GREAT

# read the metafile inputs
declare -a inputs
while IFS= read -r path; do
  inputs+=("$path")
done < "${meta}"

# list of fixed beds
declare -a fixedBeds

# if the mode is default
if [ "$mode" == "default" ]; then

    # for each file provided
    for input in "${inputs[@]}"; do
        echo "$input"
        base=$(basename "$input")
        base="${base%%.*}"  # removes everything after first *
        mkdir -p "GREAT/${base}"
        mkdir -p "GREAT/${base}/${ass_rule}"

        # fix the format of the .bed files
        fixedBed="GREAT/${base}/${ass_rule}/TEMP_${base}_fixed.bed"
        awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $1 "/" $2 "/" $3}' $input > $fixedBed
        fixedBeds+=("$fixedBed")
    done

    # fix format of fixedBeds so that it can be used in R
    fixedBedsString=$(printf "\"%s\"," "${fixedBeds[@]}")
    fixedBedsString="c(${fixedBedsString%,})"

    # write the R script
    cat >"${folder}/PBS/${mode}_${ass_rule}_GREAT.R" <<EOF
library(rGREAT)

files <- $fixedBedsString

for (file in files) {
    base <- dirname(file)

    bed <- read.delim(file)

    # try catch statement for job submission. If null move to next file
    job <- tryCatch(
        {
            submitGreatJob(bed, species = "hg38")
        },
        error = function(e) {
            cat("Error in submitting job for", file, ":", e$message, "\n")
            next
        }
    )
    if (is.null(job)) {
        next
    }

    job_info <- capture.output(job)
    writeLines(job_info, paste0(base, "/job_info.txt"))

    pdf(file = paste0(base, "/region_gene_associations.pdf"), width = 9, height = 6)
    plotRegionGeneAssociations(job)
    dev.off()

    pdf(file = paste0(base, "/volcano_plot_GO_BP.pdf"))
    plotVolcano(job, ontology = "GO Biological Process")
    dev.off()

    pdf(file = paste0(base, "/volcano_plot_GO_MF.pdf"))
    plotVolcano(job, ontology = "GO Molecular Function")
    dev.off()

    pdf(file = paste0(base, "/volcano_plot_GO_CC.pdf"))
    plotVolcano(job, ontology = "GO Cellular Component")
    dev.off()

    associations <- getRegionGeneAssociations(job)
    associations_df <- as.data.frame(associations)
    associations_df\$annotated_genes <- sapply(associations_df\$annotated_genes, paste, collapse=";")
    associations_df\$dist_to_TSS <- sapply(associations_df\$dist_to_TSS, paste, collapse=";")
    write.csv(associations_df, file = paste0(base, "/region_gene_associations.csv"), row.names = FALSE)

    Sys.sleep(60) # because GREAT limits job submissions to once every 60 seconds
}
EOF
fi

# if the mode is diffBind
if [ "$mode" == "diffBind" ]; then
    mkdir -p "GREAT/diffBind"
    states=("upregulated" "downregulated" "unchanged")

    # for each folder containing upregulated/downregulated/unchanged .bed files
    for input in "${inputs[@]}"; do
        echo "$input"
        base=$(basename "$input")
        mkdir -p "GREAT/diffBind/${base}"

        # for each state
        for state in "${states[@]}"; do
            mkdir -p "GREAT/diffBind/${base}/${state}"
            mkdir -p "GREAT/diffBind/${base}/${state}/${ass_rule}"
            
            # fix the format of the .bed files
            input_bed="${input}/${state}.bed"
            fixedBed="GREAT/diffBind/${base}/${state}/${ass_rule}/TEMP_${base}_fixed.bed"
            awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $1 "/" $2 "/" $3}' $input_bed > $fixedBed
            fixedBeds+=("$fixedBed")
        done
    done

    # fix format of fixedBeds so that it can be used in R
    fixedBedsString=$(printf "\"%s\"," "${fixedBeds[@]}")
    fixedBedsString="c(${fixedBedsString%,})"

    # write the R script
    cat >"${folder}/PBS/${mode}_${ass_rule}_GREAT.R" <<EOF

# library(rGREAT)

files <- $fixedBedsString

for (file in files) {
    path <- dirname(file)

    bed <- read.delim(file)
    job <- submitGreatJob(bed, species = "hg38")

    job_info <- capture.output(job)
    writeLines(job_info, paste0(path, "/job_info.txt"))

    pdf("${folder}/GREAT/${ass_rule}/region_gene_associations.pdf", width = 9, height = 6)
    plotRegionGeneAssociations(job)
    dev.off()


}
EOF
fi

# actually run the R script
cat >"${folder}/PBS/${mode}_${ass_rule}_GREAT.pbs" <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${mode}_${ass_rule}_GREAT  # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=32GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=10:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${mode}_${ass_rule}_GREAT.%j.out
#SBATCH --error=${folder}/log/${mode}_${ass_rule}_GREAT.%j.err
################################
# Enter your code to run below #
################################
source activate GREAT2

cd ${folder}

Rscript ${folder}/PBS/${mode}_${ass_rule}_GREAT.R

EOF
sbatch "${folder}/PBS/${mode}_${ass_rule}_GREAT.pbs"