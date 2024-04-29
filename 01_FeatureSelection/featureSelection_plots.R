# Load necessary libraries
library(tidyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)

# Read the data
df <- read.csv(file="Boxplot_AUC_featureSelection.csv", header=T, fill=T)
df <- as.data.frame(df)

# Reshape the data to long format
df_long <- reshape2::melt(df, id.vars = "V1")

# Define a custom theme for the plots
custom_theme <- theme_minimal() +
  theme(
    axis.text = element_text(colour = "black", size = 11, family = 'Arial'),
    axis.title = element_text(family = "Arial", size = 11),
    legend.text = element_text(size = 11, family = "Arial"),
  )

# Plotting Boxplot
p1 <- ggplot(df_long, aes(x = V1, y = value)) +
  geom_boxplot() + custom_theme +
  labs(title = NULL, x = "Feature Extraction Methods", y = "AUC") +
  theme_minimal() + coord_flip() + 
  theme(axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y = element_text(color = "black", size = 10))

# Save the plot as png
png("FeatureSelection_Fig1_SupplementaryMaterial1.png", width = 1200, height = 850, res = 300)
p1
dev.off()

# Read the AUC data for dataset 1
AUC1 <- read.csv("AUCs_step1.csv", header=T)
AUC1$frag_t = factor(AUC1$frag, levels=c('500 bp','1000 bp','3000 bp','5000 bp','10000 bp'))

# Plotting AUC for dataset 1
p2 <- ggplot(AUC1, aes(x=type, y=mean, group=frag, color=frag_t)) +  
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0.05)) + 
  geom_line() + geom_point() + theme_linedraw() + facet_grid(frag_t ~ .)  + ggtitle("A)") + theme(strip.text.y.left = element_text(angle = 0))

# Adding labels and themes to the plot
p3 <- p2 + labs(x=NULL, y = "mean AUC", fill = NULL, color = NULL) +  
  theme(axis.text = element_text(colour = "black", size=11,family='Arial')) 

p5 <- p3 + theme(strip.text = element_text(colour = 'black')) + 
  theme(strip.background =element_rect(fill="white")) + 
  theme(axis.title = element_text(family = "Arial", size = 11), 
        axis.text = element_text(family = "Arial"), 
        axis.text.x = element_text(family = "Arial"), 
        legend.text = element_text(size = 11, family = "Arial"))+
  scale_color_manual(values=c("500 bp"="#E41A1C","1000 bp"="#984EA3","3000 bp"="#4DAF4A","5000 bp"="#084594","10000 bp"="#BF5B17"))

# Adjusting fonts and legend position
db1 <- ggpar(p5, font.x = c(10,"black"),font.y = c(10,"black"),font.y.tickslab = c(4,"black"), legend = "top")

# Save the plot as png
png("AUCS_dataset1.png", width = 1200, height = 800, res = 200)
db1
dev.off()

# Read the AUC data for dataset 2
AUC1 <- read.csv("AUCs_step2.csv", header=T)
AUC1$frag_t = factor(AUC1$frag, levels=c('500 bp','1000 bp','3000 bp','5000 bp','10000 bp'))

# Plotting AUC for dataset 2
p2 <- ggplot(AUC1, aes(x=type, y=mean, group=frag, color=frag_t)) +  
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0.05)) + 
  geom_line() + geom_point() + theme_linedraw() + facet_grid(frag_t ~ .)  + ggtitle("B)") + theme(strip.text.y.left = element_text(angle = 0))

# Adding labels and themes to the plot
p3 <- p2 + labs(x=NULL, y = "mean AUC", fill = NULL, color = NULL) +  
  theme(axis.text = element_text(colour = "black", size=11,family='Arial')) 

p5 <- p3 + theme(strip.text = element_text(colour = 'black')) + 
  theme(strip.background =element_rect(fill="white")) + 
  theme(axis.title = element_text(family = "Arial", size = 11), 
        axis.text = element_text(family = "Arial"), 
        axis.text.x = element_text(family = "Arial"), 
        legend.text = element_text(size = 11, family = "Arial"))+
  scale_color_manual(values=c("500 bp"="#E41A1C","1000 bp"="#984EA3","3000 bp"="#4DAF4A","5000 bp"="#084594","10000 bp"="#BF5B17"))

# Adjusting fonts and removing legend
db2 <- ggpar(p5, font.x = c(10,"black"),font.y = c(10,"black"),font.y.tickslab = c(4,"black"), legend = "none")

# Save the plot as png
png("AUCS_dataset2.png", width = 1200, height = 800, res = 200)
db2
dev.off()

# Save combined plots as png
png("SupplementaryMaterial1_AUCs.png", width = 1300, height = 1600, res = 200)
grid.arrange(db1, db2, ncol = 1, heights = c(1.1,1))  # Adjust widths as needed
dev.off()
