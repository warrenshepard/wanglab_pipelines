#!/bin/bash

suffix=.sorted.filtered.bam

# Setup: Create needed folders
mkdir -p PBS
mkdir -p log
mkdir -p heatmaps
mkdir -p heatmaps_normalized
mkdir -p bigwig
mkdir -p macs2/narrowPeaks
mkdir -p macs2/broadPeaks

folder=$(pwd)
count=$(find ./aligned -mindepth 1 -type f -name "*${suffix}" -printf x | wc -c)  # Finds total number of files matching extension.

echo "${count} peakcalling to be done..."

# Meta file to know when qc will begin (empty file)
cat >${folder}/'meta.txt' <<EOF
EOF

cd aligned
for file in *_IgG*${suffix}; do
    igGbase=$(basename "$file" "${suffix}")
    groupprefix=${igGbase%%_IgG*}
	echo "${igGbase} completed!" >> ${folder}/'meta.txt'
    for f in $folder/aligned/${groupprefix}*${suffix}; do
        base=$(basename "$f" "${suffix}")
        if [ "$base" != "$igGbase" ]; then
            cat >${folder}/PBS/${base}_macs2'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_macs2  # Name of the job

# Number of compute nodes
#SBATCH --nodes=1

# Number of cores
#SBATCH --cpus-per-task=1

#Number of memory
#SBATCH --mem-per-cpu=32GB

# Number of cores, in this case one
#SBATCH --ntasks-per-node=1

# Walltime (job duration)
#SBATCH --time=24:00:00

# Name of the output files to be created. If not specified the outputs will be joined
#SBATCH --output=${folder}/log/${base}_macs2.%j.out
#SBATCH --error=${folder}/log/${base}_macs2.%j.err
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate peakcalling
cd ${folder}

# Narrow Peak call
macs2 callpeak \
-t $folder/aligned/${base}${suffix} \
-c $folder/aligned/${igGbase}${suffix} \
-f BAMPE \
-n ${base}_narrow \
--outdir $folder/macs2

# Repeat with Broad peak
macs2 callpeak \
-t $folder/aligned/${base}${suffix} \
-c $folder/aligned/${igGbase}${suffix} \
-f BAMPE \
--broad \
-n ${base}_broad \
--outdir $folder/macs2 

mv macs2/${base}_narrow_peaks.narrowPeak macs2/narrowPeaks/${base}_narrow_peaks.narrowPeak
mv macs2/${base}_broad_peaks.broadPeak macs2/broadPeaks/${base}_broad_peaks.broadPeak

conda activate alignment
cd ${folder}

bamCompare \
    -b1 aligned/${base}.sorted.filtered.bam \
	-b2 aligned/${igGbase}.sorted.filtered.bam \
	--scaleFactorsMethod readCount \
	--operation subtract \
	--binSize 10 \
	--extendReads \
	--ignoreDuplicates \
	-of bigwig \
	-o bigwig/${base}.bigwig

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R macs2/narrowPeaks/${base}_narrow_peaks.narrowPeak \
	-S bigwig/${base}.bigwig  \
	-o $folder/heatmaps_normalized/${base}_narrow_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/heatmaps_normalized/${base}_narrow_center.gz \
	--zMin 0 \
	-out $folder/heatmaps_normalized/${base}_narrow_center.pdf --colorList 'white,darkred'

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R macs2/broadPeaks/${base}_broad_peaks.broadPeak \
	-S bigwig/${base}.bigwig  \
	-o $folder/heatmaps_normalized/${base}_broad_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/heatmaps_normalized/${base}_broad_center.gz \
	--zMin 0 \
	-out $folder/heatmaps_normalized/${base}_broad_center.pdf --colorList 'white,darkred'

rm $folder/heatmaps_normalized/${base}_narrow_center.gz
rm $folder/heatmaps_normalized/${base}_broad_center.gz

echo "${base} completed!" >> ${folder}/'meta.txt'

currLine=\$(wc -l < ${folder}/meta.txt)
if ((\$currLine == $count)); then
    source activate base
    sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/CR/merge_macs2_cr.sh
    rm ${folder}/meta.txt
fi
EOF
        sbatch ${folder}/PBS/${base}_macs2'.pbs'
        fi
    done
done




