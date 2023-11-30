#!/bin/bash

suffix=.sorted.filtered.bam

# Setup: Create needed folders
mkdir -p PBS
mkdir -p log
mkdir -p mergedBam
mkdir -p mergedMacs2/narrowPeaks
mkdir -p mergedMacs2/broadPeaks
mkdir -p mergedBigwig

folder=$(pwd)
count=$(find ./aligned -mindepth 1 -type f -name "*_R1${suffix}" -printf x | wc -c)  # Finds total number of files matching extension.

# Meta file to know when qc will begin (empty file)
cat >${folder}/'mergeMeta.txt' <<EOF
EOF

cd aligned
for file in *_IgG*${suffix}; do
    igGbase=$(basename "$file" "${suffix}")
    groupprefix=${igGbase%%_IgG*}
    if [ "$igGbase" = "${groupprefix}_IgG_R1" ]; then
		echo "${igGbase} completed!" >> ${folder}/'mergeMeta.txt'
	fi
    for f in $folder/aligned/${groupprefix}*_R1${suffix}; do
        base=$(basename "$f" "_R1${suffix}")
		echo $base
        if [ "$base" != "$igGbase" ]; then
            cat >${folder}/PBS/${base}_merge_macs2'.pbs' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${base}_merge  # Name of the job

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
#SBATCH --output=${folder}/log/${base}_mergedMacs2.%j.out
#SBATCH --error=${folder}/log/${base}_mergeMacs2.%j.err
################################
# Enter your code to run below #
################################

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate alignment
cd ${folder}

samtools merge -f mergedBam/${base}_merged.bam aligned/${base}_R1${suffix} aligned/${base}_R2${suffix} 
samtools index mergedBam/${base}_merged.bam

samtools sort mergedBam/${base}_merged.bam -o mergedBam/${base}_sortedMerged.bam
samtools index mergedBam/${base}_sortedMerged.bam

rm mergedBam/${base}_merged.bam
rm mergedBam/${base}_merged.bam.bai

conda activate peakcalling

# Narrow Peak call
macs2 callpeak \
-t $folder/mergedBam/${base}_sortedMerged.bam \
-c $folder/aligned/${igGbase}${suffix} \
-f BAMPE \
-n ${base}_narrow \
--outdir $folder/mergedMacs2

# Repeat with Broad peak
macs2 callpeak \
-t $folder/mergedBam/${base}_sortedMerged.bam \
-c $folder/aligned/${igGbase}${suffix} \
-f BAMPE \
--broad \
-n ${base}_broad \
--outdir $folder/mergedMacs2

mv mergedMacs2/${base}_narrow_peaks.narrowPeak mergedMacs2/narrowPeaks/${base}_narrow_peaks.narrowPeak
mv mergedMacs2/${base}_broad_peaks.broadPeak mergedMacs2/broadPeaks/${base}_broad_peaks.broadPeak

conda activate alignment
bamCompare \
    -b1 mergedBam/${base}_sortedMerged.bam \
	-b2 aligned/${igGbase}.sorted.filtered.bam \
	--scaleFactorsMethod readCount \
	--operation subtract \
	--binSize 10 \
	--extendReads \
	--ignoreDuplicates \
	-of bigwig \
	-o mergedBigwig/${base}_merged.bigwig



cd ${folder}

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R mergedMacs2/narrowPeaks/${base}_narrow_peaks.narrowPeak \
	-S mergedBigwig/${base}_merged.bigwig  \
	-o $folder/mergedHeatmaps/${base}_narrow_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/mergedHeatmaps/${base}_narrow_center.gz \
	--zMin 0 \
	-out $folder/mergedHeatmaps/${base}_narrow_center.pdf --colorList 'white,darkred'

computeMatrix reference-point \
	--referencePoint center \
	-b 5000 -a 5000 \
	-R mergedMacs2/broadPeaks/${base}_broad_peaks.broadPeak \
	-S mergedBigwig/${base}_merged.bigwig  \
	-o $folder/mergedHeatmaps/${base}_broad_center.gz --missingDataAsZero -p max

plotHeatmap -m $folder/mergedHeatmaps/${base}_broad_center.gz \
	--zMin 0 \
	-out $folder/mergedHeatmaps/${base}_broad_center.pdf --colorList 'white,darkred'

rm $folder/mergedHeatmaps/${base}_narrow_center.gz
rm $folder/mergedHeatmaps/${base}_broad_center.gz

echo "${base} completed!" >> ${folder}/'mergeMeta.txt'

currLine=\$(wc -l < ${folder}/mergeMeta.txt)
if ((\$currLine == $count)); then
    source activate base
	sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/homerMacs2Motif.sh mergedMacs2/narrowPeaks
	sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/homerMacs2Motif.sh mergedMacs2/broadPeaks
	sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/ChIPseekerMacs2.sh mergedMacs2/narrowPeaks
	sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/ChIPseekerMacs2.sh mergedMacs2/broadPeaks
    rm ${folder}/mergeMeta.txt
fi
EOF
        sbatch ${folder}/PBS/${base}_merge_macs2'.pbs'
        fi
    done
done




