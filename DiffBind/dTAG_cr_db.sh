#!/bin/bash

# The log fold change will be calculated as log(samp/control), or log(samp) - log(control)

# Script that creates a diffBind comparison, and completes the standard analysis

# This pipeline requires four positional arguments:
# 	$1 - The CSV for the diffbind object

# NOTE: This assumes that there are only two reps. If there are more than two reps, then you must manually change that parameter within the script

# Requires the ATACQC pipeline (will rename in the future)

compare="Condition"
caps=$(echo "$compare" | tr '[:lower:]' '[:upper:]')
control="UT"
normType="default"
norm="DBA_NORM_DEFAULT"
dbObj=$(basename "$1" ".csv")
array=("4H" "24H" "72H" "30D") # Comparisons for dTAG

if [ -z "$1" ]; then 
  echo ERROR: Select csv file 
  exit 1
fi

for samp in ${array[@]}; do
off=FALSE

# Set up necessary files
mkdir -p PBS
mkdir -p log
subfolder="diffBind/${dbObj}/${samp}_over_${control}_${normType}/"
mkdir -p ${subfolder}

folder=$(pwd)

cat >${folder}/PBS/${samp}_over_${control}_${normType}_$dbObj'.R' <<EOF
# Load libraries
library(DiffBind)
library(dplyr)
library(ggplot2)
library(ggtext)

rawSamples <- read.csv('$1')

samples <- rawSamples[rawSamples\$$compare == "$samp" | rawSamples\$$compare == "$control", ]

if (!file.exists("diffBind/${dbObj}_${normType}_PCA.pdf") | !file.exists("diffBind/${dbObj}_${normType}_multicorr.pdf")) {
  dbObj <- dba(sampleSheet=rawSamples)
  dbObj <- dba.blacklist(dbObj, blacklist = DBA_BLACKLIST_GRCH38)
  dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel=FALSE)  # Does not run in parallel due to file size
  dbObj <- dba.normalize(dbObj, method = DBA_DESEQ2, normalize=${norm}, offsets=${off})
  
  corr <- plot(dbObj)
  pdf("diffBind/${dbObj}_${normType}_multicorr.pdf")
  plot(dbObj)
  dev.off()
  saveRDS(corr,file="diffBind/multicorr.RDS")
  dba.save(dbObj,file=paste0("${dbObj}","_all_","${normType}"), dir="diffBind")

  pdf("diffBind/${dbObj}_${normType}_PCA.pdf")
    dba.plotPCA(dbObj,DBA_${caps},label=DBA_ID)
  dev.off()
}

dbObj <- dba(sampleSheet=samples)
dbObj <- dba.blacklist(dbObj, blacklist = DBA_BLACKLIST_GRCH38)

dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE, bParallel=FALSE)  # Does not run in parallel due to file size
dbObj <- dba.normalize(dbObj, method = DBA_DESEQ2, normalize=${norm}, offsets=${off})

dbObj <- dba.contrast(dbObj, contrast=c("${compare}","${samp}","${control}"), design = "~${compare}", minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)

dba.save(dbObj,file="${samp}_over_${control}_${normType}", dir="diffBind/${dbObj}/${samp}_over_${control}_${normType}/")

res_all <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
res_all

out <- as.data.frame(res_all) 
write.table(out, file="${subfolder}/${samp}_over_${control}.txt", sep="\t", quote=F, row.names=F)

res_deseq <- as.data.frame(dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1))

# function for making volcanos
ggvolc <- function(data1,
                   size_var = NULL,  # Default value set to NULL
                   fc = 0.05,
                   not_sig_color = "grey82",
                   down_reg_color = "#00798c",
                   up_reg_color = "#d1495b") {
  
  # Validate input
  if(!is.data.frame(data1)) stop("data1 must be a data frame")
  
  # Calculate the size aesthetic outside ggplot
  if(is.null(size_var)) {
    data1\$size_aes <- 1  # Default size if size_var is NULL
    size_aes_range <- c(3, 3)
  } else if (size_var == "Fold") {
    data1\$size_aes <- abs(-log10(data1\$Fold))
    size_aes_range <- c(0, 6)
  } else {
    data1\$size_aes <- abs(data1[[size_var]])
    size_aes_range <- c(min(abs(data1[[size_var]])), max(abs(data1[[size_var]])))
  }
  
  dat1 <- data1 %>%
    dplyr::mutate(threshold = factor(case_when(
      Fold > 0 & FDR <= fc ~ "s_upregulated",
      Fold < 0 & FDR <= fc ~ "s_downregulated",
      TRUE ~ "not_significant"
    ), levels = c("not_significant", "s_downregulated", "s_upregulated")))
  
  dat1.2 <- dat1
  
  color_mapping <- c("s_downregulated" = down_reg_color,
                     "not_significant" = not_sig_color,
                     "s_upregulated" = up_reg_color)
  
  p <- ggplot2::ggplot(dplyr::arrange(dat1.2, threshold)) +
    ggplot2::geom_point(aes(x = Fold, y = log10FDR, color = threshold, size = size_aes),
                        shape = 16, alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "DESeq2 pVolcano",
                  x = "Fold",
                  y = "-log10(FDR)") +
    ggplot2::scale_color_manual(
      values = color_mapping,
      name = "Genes",
      breaks = c("s_downregulated", "not_significant", "s_upregulated"),
      labels = c("Downregulated", "non-significant", "Upregulated")
    )  +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5, alpha=1)))
  
  if (is.null(size_var)) {
    p <- p + scale_size_continuous(guide = "none")  # No legend for size when size_var is NULL
  } else {
    size_legend_name <- ifelse(size_var == "Fold", "Fold", "-log10(Fold)")
    p <- p + scale_size_continuous(name = size_legend_name,
                                   range = size_aes_range) +
      guides(size = guide_legend(override.aes = list(shape = 21, fill = NA)))
  }
  
  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size=14),
      axis.text.x = ggplot2::element_text(margin = margin(t = 2.5, r =0, b = 0, l = 0)),
      axis.text.y = ggplot2::element_text(margin = margin(t = 0, r =2.5, b = 0, l = 0)),
      axis.ticks.length.x = grid::unit(0.25,"cm"),
      axis.ticks.length.y = grid::unit(0.25,"cm"),
      axis.ticks = ggplot2::element_line(color = "#333333", size = .5),
      axis.title = ggplot2::element_text(size=15),
      panel.grid.major = ggplot2::element_line(color=NA),
      panel.border = ggplot2::element_rect(size = 1, color="#333333"),
      legend.title = ggplot2::element_text(hjust=0.5, size=12),
      legend.text = ggplot2::element_text(size=10),
      plot.title = ggtext::element_markdown(color = "#333333", size = 18, face = "bold", margin = margin(0,0,0.5,0, unit = "cm"), hjust=0.5),
      plot.subtitle = ggtext::element_markdown(color = "grey30", size = 12, lineheight = 1.35, hjust=0.5),
      plot.caption = ggtext::element_markdown(color = "grey30", size = 10, lineheight = 1.35, hjust=0.5)
    )
  
  return(p)
}

res_deseq <- res_deseq %>%
  mutate(log10FDR = -log10(FDR))

volcano <- ggvolc(res_deseq)
ggsave(filename = "${subfolder}/DESeq2_test_volcano.pdf", plot = volcano, width = 12, height = 10)


# Create a column to denote which genes are upregulated, downregulated, or not significantly changed
res_deseq <- res_deseq %>% 
	mutate( threshold_KD = case_when(Fold > 0 & FDR <= 0.05 ~ "Upregulated",
	Fold < 0 & FDR <= 0.05 ~ "Downregulated",
	TRUE ~ "Unchanged"))

## Create the volcano plot w/o any labels
p2 <- ggplot(res_deseq, aes(Fold, -log10(FDR))) +
		geom_point(aes(color = threshold_KD), size = 2/5) +
		xlab(expression("log"[2]*"FC")) +  
		ylab(expression("-log"[10]*"FDR")) +
		scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
		guides(colour = guide_legend(override.aes = list(size=1.5))) +
		theme_bw() + theme(legend.position = "none")

pdf("${subfolder}/DESeq2_pVolcano.pdf")
	print(p2)
dev.off()


unchanged <- as.data.frame(res_all) %>%
  filter(FDR > 0.05) %>% 
  select(seqnames, start, end)

write.table(unchanged, file="${subfolder}/unchanged_long.bed", sep="\t", quote=F, row.names=F)

allDB <- out %>% 
  filter((Fold < 0 & FDR < 0.05) | (Fold > 0 & FDR < 0.05)) %>%
  select(seqnames, start, end, Fold, p.value, FDR)

allDBbed <- allDB %>% select(seqnames, start, end)

write.table(allDBbed, file="${subfolder}/allDB.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(allDB, file="${subfolder}/allDB.txt", sep="\t", quote=F, row.names=F, col.names=F)

KO_enrich <- out %>% 
  filter(Fold < 0 & FDR < 0.05) %>% 
  select(seqnames, start, end)
  
# Write to bed file
write.table(KO_enrich, file="${subfolder}/downregulated.bed", sep="\t", quote=F, row.names=F, col.names=F)

WT_enrich <- out %>% 
  filter(Fold > 0 & FDR < 0.05) %>% 
  select(seqnames, start, end)
  
# Write to bed file
write.table(WT_enrich, file="${subfolder}/upregulated.bed", sep="\t", quote=F, row.names=F, col.names=F)

# Note: These have been moved to the bottom because if one group is zero, then the program automatically terminates
png(filename = "${subfolder}/venn_DBvsNonDB.png", height = 1080, width = 1080)
dba.plotVenn(dbObj, contrast=1, main="${samp} over ${control}", bDB=TRUE, bNotDB=TRUE, method = DBA_DESEQ2)
dev.off()

png(filename = "${subfolder}/venn.png", height = 1080, width = 1080)
dba.plotVenn(dbObj, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE, method = DBA_DESEQ2)
dev.off()
EOF

cat >${folder}/PBS/${samp}_over_${control}_${normType}_$dbObj.sbatch <<EOF
#!/bin/bash -l
#SBATCH --job-name=${samp}_over_${control}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7:00:00
#SBATCH -o ${folder}/log/${samp}_over_${control}_${dbObj}_%j.txt -e ${folder}/log/${samp}_over_${control}_${dbObj}_%j.err.txt

#------- END OF HEADER -------#
source activate ATACQC

cd ${folder}

Rscript ./PBS/${samp}_over_${control}_${normType}_$dbObj'.R'

tail -n +2 "${subfolder}/unchanged_long.bed" > ${subfolder}/unchanged.bed
rm ${subfolder}/unchanged_long.bed

sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/pipes/db_homer.sh ${subfolder}
sh /dartfs-hpc/rc/lab/W/WangX/Nicholas/backend/downstreamDiffBind/diffBind_ChIPseeker.sh ${subfolder}
EOF

sbatch ${folder}/PBS/${samp}_over_${control}_${normType}_$dbObj.sbatch 
done