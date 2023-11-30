#!/bin/bash

# Plots the EChO results from EChO.sh. Should automatically trigger

# Script automatically generated using the "slouch" command on Thu Aug 18 16:21:35 EDT 2022

folder=$(pwd)  # Stores current folder as a variable

# Set up required files
mkdir -p PBS
mkdir -p log
mkdir -p merged_EChO_results/analysis

# For each sample group (determined by IgG control)
for file in aligned/*_IgG*.sorted.filtered.bam; do
    groupPrefix=${file%%_IgG*.sorted.filtered.bam}  # 
	groupPrefix=${groupPrefix#"aligned/"}
    echo ${groupPrefix}
    cat >${folder}/PBS/${groupPrefix}_EChO'.R' <<EOF
library(ggplot2)
library(dplyr)
library(tidyr)

getwd()
df <- data.frame(matrix(ncol = 2, nrow = 0))
EOF
    # For each sample in that group
    for f in ./merged_EChO_results/${groupPrefix}*.EChO.bed; do
        echo $f
        base=$(basename "$f" "_foci.EChO.bed")
        base=${base#"./merged_EchO_results/CR_O422T"}
        echo ${base}
            # Add that file to the overall dataframe
            cat >>${folder}/PBS/${groupPrefix}_EChO'.R' <<EOF

df2 <- read.table("$f", sep = "\t")
df2 <- data.frame(sampleGroup="${base}", fragSize = df2\$V4)
str(df2)

df <- rbind(df, df2)
EOF
    done

    # Create the graph
    cat >>${folder}/PBS/${groupPrefix}_EChO'.R' <<EOF

number_groups <- length(unique(df\$sampleGroup))
width <- (2.5 * number_groups) - 2.5
if(number_groups == 2) {
    width <- 5
}
df\$sampleGroup <- as.factor(df\$sampleGroup)
df\$sampleGroup <- factor(df\$sampleGroup, levels = rev(levels(df\$sampleGroup)))
print(levels(df\$sampleGroup))

summary <- df %>%
  group_by(sampleGroup) %>%
  summarise(
    below = sum(fragSize <= 150),
    above = sum(fragSize > 150),
    percent_above = (above / n()) * 100,
    percent_below = (below / n()) * 100
  )

long_summary <- summary %>%
  gather(key = "Metric", value = "Value", -sampleGroup)
long_summary <- long_summary[long_summary\$Metric %in% c("percent_above", "percent_below"),]

pointSize <- 16
lineWidth <- 1 / 2.835

bar <- ggplot(long_summary, aes(x = Metric, y = Value, fill = sampleGroup)) + 
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.2 * length(unique(df\$sampleGroup)), size = 1.2) + 
  labs(y = "Percent", x = "Sample Group") + 
  scale_x_discrete(labels = c("percent_above" = "Above Threshold", 
                              "percent_below" = "Below Threshold")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + 
  theme(
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(size = pointSize * 0.8, colour = "black"),
    axis.title  = element_text(size = pointSize * 0.8, colour = "black"),
    axis.text.x  = element_text(
      size = pointSize * 0.6,
      colour = "black",
      hjust = .5
    ),
    axis.text.y  = element_text(size = pointSize * 0.6, colour = "black"),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize * 0.6, colour = "black"),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.2, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    strip.placement = "outside"
  ) 
ggsave("merged_EChO_results/analysis/${groupPrefix}_bar.pdf", plot = bar, width = 3.25, height = 6)

#png("merged_EChO_results/analysis/${groupPrefix}.png", width = 1080, height = 720)
    plot <- ggplot(df, aes(x=sampleGroup, y=fragSize, fill = sampleGroup)) + 
    geom_hline(yintercept = 147, linetype="dashed", size = 1.2, color = "red") +
    geom_violin() + 
    ggtitle("${groupPrefix} Fragment Size Distributions") + 
    geom_boxplot(width = 0.2, fill="white", outlier.shape = NA, alpha = .7) + 
    scale_y_continuous(limits = c(0, 800)) +
    theme(
    axis.ticks = element_blank(),
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(size = pointSize * 0.8, colour = "black"),
    axis.title  = element_text(size = pointSize * 0.8, colour = "black"),
    axis.text.x  = element_text(
      size = pointSize * 0.6,
      colour = "black",
      angle = 0,
      hjust = .5
    ),
    # axis.text.y  = element_text(size = pointSize * 0.6, colour = "black"),
    axis.text.y  = element_blank(),
    # axis.line.y = element_blank(),
    # panel.grid.major = element_line(color = "lightgrey"),
    # panel.grid.minor = element_line(color = "lightgrey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.background = element_blank(), 
    # plot.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize * 0.75, colour = "black"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.4, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    strip.text = element_text(
      face = "italic",
      size = pointSize * 0.6,
      hjust = 0
    ),
    strip.placement = "outside"
  ) 
ggsave("merged_EChO_results/analysis/${groupPrefix}.pdf", plot = plot, width = width, height = 6)

plot2 <- ggplot(df, aes(x=sampleGroup, y=fragSize, fill = sampleGroup)) + 
    geom_violin() + 
    ggtitle("${groupPrefix} Fragment Size Distributions") + 
    geom_boxplot(width = 0.2, fill="white", outlier.shape = NA, alpha = .7) + 
    geom_hline(aes(yintercept=147)) +
    stat_summary(
    fun = median, 
    geom = "text", 
    aes(label = round(..y.., 2)),  # round to 2 decimal places
    vjust = -.5,
    size = 3
    ) +
    theme(
    axis.ticks = element_blank(),
    text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(size = pointSize * 0.8, colour = "black"),
    axis.title  = element_text(size = pointSize * 0.8, colour = "black"),
    axis.text.x  = element_text(
      size = pointSize * 0.6,
      colour = "black",
      angle = 0,
      hjust = .5
    ),
    # axis.text.y  = element_text(size = pointSize * 0.6, colour = "black"),
    axis.text.y  = element_blank(),
    # axis.line.y = element_blank(),
    # panel.grid.major = element_line(color = "lightgrey"),
    # panel.grid.minor = element_line(color = "lightgrey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.background = element_blank(), 
    # plot.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = pointSize * 0.75, colour = "black"),
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.4, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    strip.text = element_text(
      face = "italic",
      size = pointSize * 0.6,
      hjust = 0
    ),
    strip.placement = "outside"
  )
ggsave("merged_EChO_results/analysis/${groupPrefix}_nums.pdf", plot = plot2, width = width, height = 6)

# dev.off()
EOF

    # Actually run the scripts
    cat >${folder}/PBS/${groupPrefix}_plotEChO'.sbatch' <<EOF
#!/bin/bash -l
# Name of the job
#SBATCH --job-name=${groupPrefix}_plotEChO  # Name of the job

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
#SBATCH --output=${folder}/log/${groupPrefix}_plotEChO.%j.out
#SBATCH --error=${folder}/log/${groupPrefix}_plotEChO.%j.err
################################
# Enter your code to run below #
################################
cd ${folder}

source /dartfs-hpc/rc/lab/W/WangX/sharedconda/miniconda/etc/profile.d/conda.sh
source activate ATACQC

Rscript ${folder}/PBS/${groupPrefix}_EChO'.R'

# rm merged_EChO_results/*${base}*frags*
EOF
    sbatch ${folder}/PBS/${groupPrefix}_plotEChO'.sbatch'
done