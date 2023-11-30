library(diffloop)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(viridis)


# Setup -------------------------------------------------------------------
# Directory from output of HiChipper (note I moved all mango files to a subfolder titled Mango)
ref <-""
mango_dir <- paste0("HiChipper",ref,"/Mango/")
# Output directory for figures
dir.create(paste0("DiffLoop",ref))
fig_path <- paste0("DiffLoop",ref,"/Figures/")
dir.create(fig_path)

# Figure saving function
SaveFigure <-
  function(plots, name, type = "png", width, height, res) {
    if (type == "png") {
      png(
        paste0(fig_path, name, ".", type),
        width = width,
        height = height,
        units = "in",
        res = 200
      )
    } else {
      pdf(paste0(fig_path, name, ".", type),
          width = width,
          height = height)
    }
    print(plots)
    dev.off()
  }

# Read in mango files and make a loops object
loops <- loopsMake.mango(mango_dir)
groups <- c("dS5_48H_K27Ac","dS5_48H_K27Ac","dS5_UT_K27Ac","dS5_UT_K27Ac")
loops <- updateLDGroups(loops, groups)


# Quality Control ---------------------------------------------------------
# Remove loops that merged on import
loops <- subsetLoops(loops, loops@rowData$loopWidth >= 5000)

# Filter loops w/ FDR > 0.01
loops_fdr <- mangoCorrection(loops, FDR=0.01)
dim(loops) # loops prior to FDR filtering
dim(loops_fdr) # loops after FDR filtering

# Filter loops inconsistent between reps 
km_filt <- loops_fdr[,c(1,2,3,4)]
cm <- km_filt@counts
k_dis <- ((cm[,1]>=5 &cm[,2]==0)|(cm[,2]>=5&cm[,1]==0))
m_dis <- ((cm[,3]>=5 &cm[,4]==0)|(cm[,4]>=5&cm[,3]==0))
qc_filt <- subsetLoops(km_filt, !(k_dis | m_dis))
dim(qc_filt) # loops after final filtering step

# Loop Metrics
loopMetrics(qc_filt)

# Plot PETs 
p <- loopDistancePlot(qc_filt)
SaveFigure(p, "QC_PETs_by_Distance", width = 6, height = 6)

# Plot PCA
p <-pcaPlot(qc_filt) + geom_text_repel(aes(label=groups)) +
  ggtitle("PC Plot with Size Factor Correction") +
  theme(legend.position="none")
SaveFigure(p, "QC_PCA", width = 6, height = 6)

# Differential Analysis ---------------------------------------------------
# Can use Voom or EdgeR, using EdgeR as that is what Ben did, but could be something to look into in the future
km_res <- quickAssoc(qc_filt)

# Annotation
CR_dir <- "../../CR/2023_09_WT/mergedMacs2/broadPeaks"
k27Ac_bed <- paste0(CR_dir, "/", "CR_HCT_WT_K27Ac_broad_peaks.broadPeak")
k4Me1_bed <- paste0(CR_dir, "/", "CR_HCT_WT_K4Me1_broad_peaks.broadPeak")
h3k27ac <- rmchr(padGRanges(bedToGRanges(k27Ac_bed), pad = 1000))
h3k4me1 <- rmchr(padGRanges(bedToGRanges(k4Me1_bed), pad = 1000))
# Defining enhancers in this case by K4me1 as well for primed enhancers, but should look at this with just k27ac as well.
# enhancer <- union(h3k27ac,h3k4me1)
enhancer <- h3k27ac

# Getting TSS data for hg38 - THIS IS SUPER ROUNDABOUT AND PROBABLY WRONG -- had to change to the exact format diffloop wants which was built for hg19
# Had to change the chr identifier which is most likely the biggest problem.
TSS <- read.table("../../../Genomes_and_extra/TSSsites_ARID1A/hg38_TSS_UCSC_paddedUnique.bed",header=FALSE)
TSS$V4 <- as.character(TSS$V4)
split_strings <- strsplit(TSS$V4, "_")
TSS$V1 <- as.character(TSS$V1)
split_chr <- strsplit(TSS$V1, "chr")

# Extract the last element from each split and assign it to the new column
TSS$Gene_ID <- sapply(split_strings, function(x) x[length(x)])
TSS$chr <- sapply(split_chr, function(x) x[length(x)])

promoter <- GRanges(
  seqnames = TSS$chr,
  ranges = IRanges(start = TSS$V2, end = TSS$V3),
  gene = TSS$Gene_ID
)

km_anno <- annotateLoops(km_res, enhancer = enhancer, promoter = promoter)

# Tables
outTable <- paste0("DiffLoop",ref,"/Tables")
dir.create(outTable)

# Filter by FDR
filter <- (km_anno@rowData$FDR < 0.01)
res_filter <- subsetLoops(km_anno, filter)
df <- summary(res_filter)

# Write all loops above FDR
full <- unique(df)
sdf <- rbind(head(df[order(df$logFC),] , 5), head(df[order(df$logFC, decreasing = TRUE),] , 5))
write.table(full[,c(1:11,14,17,18,20)], file = paste0(outTable,"/AllDiff.csv"), sep = ",", quote = FALSE, row.names = FALSE)
write.table(sdf[,c(1:11,14,17,18,20)], file = paste0(outTable,"/TopDiff.csv"), sep = ",", quote = FALSE, row.names = FALSE)

# Keep EP loops
km_res.ep <- keepEPloops(km_res, enhancer, promoter)
df2 <- summary(km_res.ep)

# Filter to e-p loops with FDR < 0.01
filter <- (km_res.ep@rowData$FDR < 0.01)
res_filter <- subsetLoops(km_res.ep, filter)
df2 <- summary(res_filter)

# Strongest Associations
full_ep <- unique(df2)
sdf <- unique(rbind(head(df2[order(df2$logFC),] , 5), head(df2[order(df2$logFC, decreasing = TRUE),] , 5)))
write.table(full_ep[,c(1:11,14,17,18,20)], file = paste0(outTable,"/AllDiff_EP.csv"), sep = ",", quote = FALSE, row.names = FALSE)
write.table(sdf[,c(1:11,14,17,18,20)], file = paste0(outTable,"/TopDiff_EP.csv"), sep = ",", quote = FALSE, row.names = FALSE)


# Loop through and create plots for each strong change in EP looping
for (row in 1:nrow(sdf)){
  # Extract values
  chr=sdf[row,"chr_1"]
  start=sdf[row,"start_1"]-10000
  end=sdf[row,"end_2"]+10000
  tss=sdf[row,"gene.tss"]
  
  # Create plot
  range <- GRanges(seqnames=c(chr),ranges=IRanges(start=c(start),end=c(end)))
  out <- paste0(fig_path,tss,".pdf")
  pdf(out)
  loopPlot(km_anno, range)
  dev.off()
}

