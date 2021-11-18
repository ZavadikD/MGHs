library(ape)
library(adegenet)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(cowplot)
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/")

##################SNAP output analysis - dn/ds (inverted ratio-INV)#########################

###put the directory with SNAP raw output files below. Output in .tab format
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/SNAP/")


###make a list of genes and initiate empty dataframe
genes <- c("hop1", "hop2", "mnd1", "dmc1", "mer3", "msh4")

snap <-data.frame("genes" <- genes)

for(i in genes){
  x_sexual <- read.csv(paste0(i,"_sexual.tab"), sep = "\t")
  x_asexual <- read.csv(paste0(i,"_asexual.tab"), sep = "\t")
  asexual <- x_asexual$dn/x_asexual$ds
  sexual <- x_sexual$dn/x_sexual$ds
  print(wilcox.test(sexual, asexual)[3])
  snap[which(snap$X.genes.....genes == i),2] <- wilcox.test(sexual, asexual)[3]
  snap[which(snap$X.genes.....genes == i),3] <- mean(na.omit(sexual))
  snap[which(snap$X.genes.....genes == i),4] <- mean(na.omit(asexual))
  
  c1_asexual <- read.csv(paste0("c1_asexual_",i,".tab"), sep = "\t")
  c1asexual <- c1_asexual$dn/c1_asexual$ds
  print(wilcox.test(sexual, c1asexual)[3])
  snap[which(snap$X.genes.....genes == i),5] <- wilcox.test(sexual, c1asexual)[3]
  snap[which(snap$X.genes.....genes == i),6] <- mean(na.omit(c1asexual))
  
  c2_asexual <- read.csv(paste0("c2_asexual_",i,".tab"), sep = "\t")
  c2asexual <- c2_asexual$dn/c2_asexual$ds
  print(wilcox.test(sexual, c2asexual)[3])
  snap[which(snap$X.genes.....genes == i),7] <- wilcox.test(sexual, c2asexual)[3]
  snap[which(snap$X.genes.....genes == i),8] <- mean(na.omit(c2asexual))
  
  sumasexual <- mean(c(na.omit(c1_asexual$dn/c1_asexual$ds), na.omit(c2_asexual$dn/c2_asexual$ds)))
  print(wilcox.test(sexual, sumasexual)[3])
  snap[which(snap$X.genes.....genes == i),9] <- wilcox.test(sexual, sumasexual)[3]
  snap[which(snap$X.genes.....genes == i),10] <- mean(na.omit(sumasexual))
}
colnames(snap) <- c("gene", "p-value_vs_all_asexual","dn.ds_mean_sexual","dn.ds_mean_all_asexual","p-value_vs_Core_asexual","dn.ds_mean_Core_asexual","p-value_vs_Prasino_asexual", "dn.ds_mean_Prasino_asexual", "p-value_vs_Core+Prasino_asexual", 
                    "dn.ds_mean_Core+Prasino_asexual")
snap$dn.ds_mean_Prasino_asexual[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"
snap$`p-value_vs_Prasino_asexual`[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"

###Here you go. A table with a summary stats on SNAP output
print(snap)
###save wherever you wish
write.csv(snap, "snapdndsINV_values&stats.csv")

##################
#now, the table requires a bit of reshaping - for more convenient plotting
##################
data <- read.csv("reshaped_snapdndsINV_values&stats.csv", sep = ",")
str(data)


data1 <- within(data,{
  p_value_factor <- NA
  p_value_factor[p.value <= 0.05] <- "Significant"
  p_value_factor[p.value >= 0.05] <- "Not significant" 
  p_value_factor[p.value == "No"] <- "Control"
  p_value_factor[p.value == "Control"] <- "Control"})
#data1

SNAPdnds <- ggplot(data=data1, aes(x=gene, y=as.numeric(dn.ds_mean), group=X, color = X, shape = p_value_factor)) +
  geom_point(size=5)+
  #+ yscale("log2", .format = TRUE)+ 
  #geom_line(position = position_dodge(width=0.5))+
  #geom_pointrange(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)),size=0.8,position = position_dodge(width=0.5))+
  #geom_errorbar(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)), width=0.2, size=0.8, position = position_dodge(width=0.5))+
  xlab("Gene") + ylab("dn/ds ratio (SNAP algorithm)") +
  labs( color = "Group", shape = "Difference sexual \nvs asexual")+ scale_color_brewer(palette = "Paired")+theme_bw(base_size = 15)+
  scale_shape_manual(values=c(7:10))

ggsave("SNAPdnds_results.png", SNAPdnds,width = 15, height = 5)

##################SNAP output analysis - pn/ps (inverted ratio-INV)#########################

###put the directory with SNAP raw output files below. Output in .tab format
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/SNAP/")


###make a list of genes and initiate empty dataframe
genes <- c("hop1", "hop2", "mnd1", "dmc1", "mer3", "msh4")

snap <-data.frame("genes" <- genes)

for(i in genes){
  x_sexual <- read.csv(paste0(i,"_sexual.tab"), sep = "\t")
  x_asexual <- read.csv(paste0(i,"_asexual.tab"), sep = "\t")
  asexual <- x_asexual$pn/x_asexual$ps
  sexual <- x_sexual$pn/x_sexual$ps
  print(wilcox.test(sexual, asexual)[3])
  snap[which(snap$X.genes.....genes == i),2] <- wilcox.test(sexual, asexual)[3]
  snap[which(snap$X.genes.....genes == i),3] <- mean(na.omit(sexual))
  snap[which(snap$X.genes.....genes == i),4] <- mean(na.omit(asexual))
  
  c1_asexual <- read.csv(paste0("c1_asexual_",i,".tab"), sep = "\t")
  c1asexual <- c1_asexual$pn/c1_asexual$ps
  print(wilcox.test(sexual, c1asexual)[3])
  snap[which(snap$X.genes.....genes == i),5] <- wilcox.test(sexual, c1asexual)[3]
  snap[which(snap$X.genes.....genes == i),6] <- mean(na.omit(c1asexual))
  
  c2_asexual <- read.csv(paste0("c2_asexual_",i,".tab"), sep = "\t")
  c2asexual <- c2_asexual$pn/c2_asexual$ps
  print(wilcox.test(sexual, c2asexual)[3])
  snap[which(snap$X.genes.....genes == i),7] <- wilcox.test(sexual, c2asexual)[3]
  snap[which(snap$X.genes.....genes == i),8] <- mean(na.omit(c2asexual))
  
  sumasexual <- mean(c(na.omit(c1_asexual$pn/c1_asexual$ps), na.omit(c2_asexual$pn/c2_asexual$ps)))
  print(wilcox.test(sexual, sumasexual)[3])
  snap[which(snap$X.genes.....genes == i),9] <- wilcox.test(sexual, sumasexual)[3]
  snap[which(snap$X.genes.....genes == i),10] <- mean(na.omit(sumasexual))
}
colnames(snap) <- c("gene", "p-value_vs_all_asexual","pn.ps_mean_sexual","pn.ps_mean_all_asexual","p-value_vs_Core_asexual","pn.ps_mean_Core_asexual","p-value_vs_Prasino_asexual", "pn.ps_mean_Prasino_asexual", "p-value_vs_Core+Prasino_asexual", 
                    "pn.ps_mean_Core+Prasino_asexual")
snap$pn.ps_mean_Prasino_asexual[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"
snap$`p-value_vs_Prasino_asexual`[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"

###Here you go. A table with a summary stats on SNAP output
print(snap)
###save wherever you wish
write.csv(snap, "snappnpsINV_values&stats.csv")

##################
#now, the table requires a bit of reshaping - for more convenient plotting
##################
data <- read.csv("reshaped_snappnpsINV_values&stats.csv", sep = ",")
str(data)


data1 <- within(data,{
  p_value_factor <- NA
  p_value_factor[p.value <= 0.05] <- "Significant"
  p_value_factor[p.value >= 0.05] <- "Not significant" 
  p_value_factor[p.value == "No"] <- "Control"
  p_value_factor[p.value == "Control"] <- "Control"})
data1

SNAPpnps <- ggplot(data=data1, aes(x=gene, y=as.numeric(pn.ps_mean), group=X, color = X, shape = p_value_factor)) +
  geom_point(size=5)+
  #+ yscale("log2", .format = TRUE)+ 
  #geom_line(position = position_dodge(width=0.5))+
  #geom_pointrange(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)),size=0.8,position = position_dodge(width=0.5))+
  #geom_errorbar(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)), width=0.2, size=0.8, position = position_dodge(width=0.5))+
  xlab("Gene") + ylab("pn/ps ratio (SNAP algorithm)") +
  labs( color = "Group", shape = "Difference sexual \nvs asexual")+ scale_color_brewer(palette = "Paired")+theme_bw(base_size = 15)+
  scale_shape_manual(values=c(7:10))

ggsave("SNAPpnps_results.png", SNAPpnps,width = 15, height = 5)

#################################################################################################################################################################################

##################dn/ds function in "ape" - calculation and analysis########################

###the very basic idea
#sequence1 <- fasta2DNAbin("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/protein-coding_regions/split_clades_for_SNAP/c1_asexual_mer3_muscle.fas")
#dnds <- dnds(sequence1)
#mean(na.omit(dnds))

#################initiating a dataframe again, and setting a temprorary working directory where alignments are kept
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/Official_alignments/")
genes <- c("hop1", "hop2", "mnd1", "dmc1", "mer3", "msh4")
ape <-data.frame("genes" <- genes)
  

for(i in genes){
  print(i)
  sexual <- fasta2DNAbin(paste0("muscle_all_sexual_",i,".fas"))
  asexual <- fasta2DNAbin(paste0("muscle_all_asexual_",i,".fas"))
  asexual_vector <- c(as.numeric(na.omit(dnds(asexual))))
  sexual_vector <- c(as.numeric(na.omit(dnds(sexual))))
  ape[which(ape$X.genes.....genes == i),2] <- wilcox.test(sexual_vector, asexual_vector)[3]
  ape[which(ape$X.genes.....genes == i),3] <- mean(sexual_vector)
  ape[which(ape$X.genes.....genes == i),4] <- mean(asexual_vector)
  
  asexual <- fasta2DNAbin(paste0("c1_asexual_", i, "_muscle.fas"))
  asexual_vector1 <- as.numeric(na.omit(dnds(asexual)))
  ape[which(ape$X.genes.....genes == i),5] <- wilcox.test(sexual_vector, asexual_vector1)[3]
  ape[which(ape$X.genes.....genes == i),6] <- mean(asexual_vector1)
  
  asexual <- fasta2DNAbin(paste0("c2_asexual_", i, "_muscle.fas"))
  asexual_vector2 <- as.numeric(na.omit(dnds(asexual)))
  ape[which(ape$X.genes.....genes == i),7] <- wilcox.test(sexual_vector, asexual_vector2)[3]
  ape[which(ape$X.genes.....genes == i),8] <- mean(asexual_vector2)
  
  sumasexual <- as.numeric(c(asexual_vector1, asexual_vector2))
  ape[which(ape$X.genes.....genes == i),9] <- wilcox.test(sexual_vector, sumasexual)[3]
  ape[which(ape$X.genes.....genes == i),10] <- mean(sumasexual)
}
ape

colnames(ape) <- c("gene", "p-value_vs_all_asexual","dn.ds_mean_sexual","dn.ds_mean_all_asexual","p-value_vs_Core_asexual","dn.ds_mean_Core_asexual","p-value_vs_Prasino_asexual", "dn.ds_mean_Prasino_asexual", "p-value_vs_Core+Prasino_asexual", 
                    "dn.ds_mean_Core+Prasino_asexual")
ape$dn.ds_mean_Prasino_asexual[c(which(ape$gene == "hop1"),which(ape$gene == "msh4"))] <- "No"
ape$`p-value_vs_Prasino_asexual`[c(which(ape$gene == "hop1"),which(ape$gene == "msh4"))] <- "No"

###Here you go. A table with a summary stats on SNAP output
print(ape)
###save wherever you wish
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/")
write.csv(ape, "ape_values&stats.csv")
##################
#now, the table requires a bit of reshaping - for more convenient plotting
##################
data <- read.csv("reshaped_ape_values&stats.csv", sep = ",")
str(data)

data1 <- within(data,{
  p_value_factor <- NA
  p_value_factor[p.value <= 0.05] <- "Significant"
  p_value_factor[p.value >= 0.05] <- "Not significant" 
  p_value_factor[p.value == "No"] <- "Control"
  p_value_factor[p.value == "Control"] <- "Control"})
#data1

apednds <- ggplot(data=data1, aes(x=gene, y=as.numeric(dn.ds_mean), group=X, color = X, shape = p_value_factor)) +
  geom_point(size=5)+
  #+ yscale("log2", .format = TRUE)+ 
  #geom_line(position = position_dodge(width=0.5))+
  #geom_pointrange(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)),size=0.8,position = position_dodge(width=0.5))+
  #geom_errorbar(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)), width=0.2, size=0.8, position = position_dodge(width=0.5))+
  xlab("Gene") + ylab("dn/ds ratio (ape algorithm)") +
  labs( color = "Group", shape = "Difference sexual \nvs asexual")+ scale_color_brewer(palette = "Paired")+theme_bw(base_size = 15)+
  scale_shape_manual(values=c(7:10))

ggsave("apednds_results.png", apednds ,width = 15, height = 5)



######plots everything together(like in the article)
shared_legend <- extract_legend(apednds)
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

alltogather <- grid.arrange(arrangeGrob(apednds+ theme(legend.position = "none"), SNAPdnds + theme(legend.position = "none"), SNAPpnps + theme(legend.position = "none"), shared_legend, ncol = 2))


ggsave("Natural_selection_results.png", alltogather ,width = 10, height = 10)







#################################################################################################################################################################################
#################################################################################################################################################################################
#####The following are for SNAP non-inverted ratio (not used in the article but can be used just in case)
##################SNAP output analysis - ds/dn########################

###put the directory with SNAP raw output files below. Output in .tab format
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/SNAP/")


###make a list of genes and initiate empty dataframe
genes <- c("hop1", "hop2", "mnd1", "dmc1", "mer3", "msh4")

snap <-data.frame("genes" <- genes)

for(i in genes){
  x_sexual <- read.csv(paste0(i,"_sexual.tab"), sep = "\t")
  x_asexual <- read.csv(paste0(i,"_asexual.tab"), sep = "\t")
  asexual <- x_asexual$ds.dn
  sexual <- x_sexual$ds.dn
  print(wilcox.test(sexual, asexual)[3])
  snap[which(snap$X.genes.....genes == i),2] <- wilcox.test(sexual, asexual)[3]
  snap[which(snap$X.genes.....genes == i),3] <- mean(na.omit(sexual))
  snap[which(snap$X.genes.....genes == i),4] <- mean(na.omit(asexual))
  
  c1_asexual <- read.csv(paste0("c1_asexual_",i,".tab"), sep = "\t")
  c1asexual <- c1_asexual$ds.dn
  print(wilcox.test(sexual, c1asexual)[3])
  snap[which(snap$X.genes.....genes == i),5] <- wilcox.test(sexual, c1asexual)[3]
  snap[which(snap$X.genes.....genes == i),6] <- mean(na.omit(c1asexual))
  
  c2_asexual <- read.csv(paste0("c2_asexual_",i,".tab"), sep = "\t")
  c2asexual <- c2_asexual$ds.dn
  print(wilcox.test(sexual, c2asexual)[3])
  snap[which(snap$X.genes.....genes == i),7] <- wilcox.test(sexual, c2asexual)[3]
  snap[which(snap$X.genes.....genes == i),8] <- mean(na.omit(c2asexual))
  
  sumasexual <- mean(c(na.omit(c1_asexual$ds.dn), na.omit(c2_asexual$ds.dn)))
  print(wilcox.test(sexual, sumasexual)[3])
  snap[which(snap$X.genes.....genes == i),9] <- wilcox.test(sexual, sumasexual)[3]
  snap[which(snap$X.genes.....genes == i),10] <- mean(na.omit(sumasexual))
}
colnames(snap) <- c("gene", "p-value_vs_all_asexual","ds.dn_mean_sexual","ds.dn_mean_all_asexual","p-value_vs_Core_asexual","ds.dn_mean_Core_asexual","p-value_vs_Prasino_asexual", "ds.dn_mean_Prasino_asexual", "p-value_vs_Core+Prasino_asexual", 
                      "ds.dn_mean_Core+Prasino_asexual")
snap$ds.dn_mean_Prasino_asexual[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"
snap$`p-value_vs_Prasino_asexual`[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"

###Here you go. A table with a summary stats on SNAP output
print(snap)
###save wherever you wish
write.csv(snap, "snapdsdn_values&stats.csv")
##################
#now, the table requires a bit of reshaping - for more convenient plotting
##################
data <- read.csv("reshaped_snapdsdn_values&stats.csv", sep = ",")
str(data)


data1 <- within(data,{
  p_value_factor <- NA
  p_value_factor[p.value <= 0.05] <- "Significant"
  p_value_factor[p.value >= 0.05] <- "Not significant" 
  p_value_factor[p.value == "No"] <- "Absent"
  p_value_factor[p.value == "Control"] <- "Control"})
#data1

ggplot(data=data1, aes(x=gene, y=ds.dn_mean, group=X, color = X, shape = p_value_factor)) +
  geom_point(size=5)+
  #+ yscale("log2", .format = TRUE)+ 
  #geom_line(position = position_dodge(width=0.5))+
  #geom_pointrange(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)),size=0.8,position = position_dodge(width=0.5))+
  #geom_errorbar(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)), width=0.2, size=0.8, position = position_dodge(width=0.5))+
  xlab("Gene") + ylab("ds/dn ratio (SNAP algorithm)") +
  labs( color = "Group", shape = "Difference sexual \nvs asexual")+ scale_color_brewer(palette = "Paired")+theme_bw(base_size = 15)+
  scale_shape_manual(values=c(7:10))


##################SNAP output analysis - ps/pn########################

###put the directory with SNAP raw output files below. Output in .tab format
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/SNAP/")


###make a list of genes and initiate empty dataframe
genes <- c("hop1", "hop2", "mnd1", "dmc1", "mer3", "msh4")

snap <-data.frame("genes" <- genes)


for(i in genes){
  x_sexual <- read.csv(paste0(i,"_sexual.tab"), sep = "\t")
  x_asexual <- read.csv(paste0(i,"_asexual.tab"), sep = "\t")
  asexual <- x_asexual$ps.pn
  sexual <- x_sexual$ps.pn
  print(wilcox.test(sexual, asexual)[3])
  snap[which(snap$X.genes.....genes == i),2] <- wilcox.test(sexual, asexual)[3]
  snap[which(snap$X.genes.....genes == i),3] <- mean(na.omit(sexual))
  snap[which(snap$X.genes.....genes == i),4] <- mean(na.omit(asexual))
  
  c1_asexual <- read.csv(paste0("c1_asexual_",i,".tab"), sep = "\t")
  c1asexual <- c1_asexual$ps.pn
  print(wilcox.test(sexual, c1asexual)[3])
  snap[which(snap$X.genes.....genes == i),5] <- wilcox.test(sexual, c1asexual)[3]
  snap[which(snap$X.genes.....genes == i),6] <- mean(na.omit(c1asexual))
  
  c2_asexual <- read.csv(paste0("c2_asexual_",i,".tab"), sep = "\t")
  c2asexual <- c2_asexual$ps.pn
  print(wilcox.test(sexual, c2asexual)[3])
  snap[which(snap$X.genes.....genes == i),7] <- wilcox.test(sexual, c2asexual)[3]
  snap[which(snap$X.genes.....genes == i),8] <- mean(na.omit(c2asexual))
  
  sumasexual <- mean(c(na.omit(c1_asexual$ps.pn), na.omit(c2_asexual$ps.pn)))
  print(wilcox.test(sexual, sumasexual)[3])
  snap[which(snap$X.genes.....genes == i),9] <- wilcox.test(sexual, sumasexual)[3]
  snap[which(snap$X.genes.....genes == i),10] <- mean(na.omit(sumasexual))
}
colnames(snap) <- c("gene", "p-value_vs_all_asexual","ps.pn_mean_sexual","ps.pn_mean_all_asexual","p-value_vs_Core_asexual","ps.pn_mean_Core_asexual","p-value_vs_Prasino_asexual", "ps.pn_mean_Prasino_asexual", "p-value_vs_Core+Prasino_asexual", 
                    "ps.pn_mean_Core+Prasino_asexual")
snap$ps.pn_mean_Prasino_asexual[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"
snap$`p-value_vs_Prasino_asexual`[c(which(snap$gene == "hop1"),which(snap$gene == "msh4"))] <- "No"

###Here you go. A table with a summary stats on SNAP output
print(snap)
###save wherever you wish
write.csv(snap, "snappspn_values&stats.csv")
##################
#now, the table requires a bit of reshaping - for more convenient plotting
##################
data <- read.csv("reshaped_snappspn_values&stats.csv", sep = ",")
str(data)


data1 <- within(data,{
  p_value_factor <- NA
  p_value_factor[p.value <= 0.05] <- "Significant"
  p_value_factor[p.value >= 0.05] <- "Not significant" 
  p_value_factor[p.value == "No"] <- "Absent"
  p_value_factor[p.value == "Control"] <- "Control"})
#data1

ggplot(data=data1, aes(x=gene, y=as.numeric(ps.pn_mean), group=X, color = X, shape = p_value_factor)) +
  geom_point(size=5)+
  #+ yscale("log2", .format = TRUE)+ 
  #geom_line(position = position_dodge(width=0.5))+
  #geom_pointrange(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)),size=0.8,position = position_dodge(width=0.5))+
  #geom_errorbar(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)), width=0.2, size=0.8, position = position_dodge(width=0.5))+
  xlab("Gene") + ylab("ps/pn ratio (SNAP algorithm)") +
  labs( color = "Group", shape = "Difference control \nvs test condition")+ scale_color_brewer(palette = "Paired")+theme_bw(base_size = 15)+
  scale_shape_manual(values=c(7:10))



##################dn/ds function in "ape" - calculation and analysis########################

###the very basic idea
#sequence1 <- fasta2DNAbin("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/protein-coding_regions/split_clades_for_SNAP/c1_asexual_mer3_muscle.fas")
#dnds <- dnds(sequence1)
#mean(na.omit(dnds))

#################initiating a dataframe again, and setting a temprorary working directory where alignments are kept
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/Official_alignments/")
genes <- c("hop1", "hop2", "mnd1", "dmc1", "mer3", "msh4")
ape <-data.frame("genes" <- genes)

1/2*
  
  sexual <- fasta2DNAbin(paste0("muscle_all_sexual_",i,".fas"))
asexual <- fasta2DNAbin(paste0("muscle_all_asexual_","mer3",".fas"))
asexual_vector <- c(as.numeric(na.omit(dnds(asexual))))

for(i in genes){
  print(i)
  sexual <- fasta2DNAbin(paste0("muscle_all_sexual_",i,".fas"))
  asexual <- fasta2DNAbin(paste0("muscle_all_asexual_",i,".fas"))
  asexual_vector <- c(as.numeric(na.omit(dnds(asexual))))
  sexual_vector <- c(as.numeric(na.omit(dnds(sexual))))
  ape[which(ape$X.genes.....genes == i),2] <- wilcox.test(sexual_vector, asexual_vector)[3]
  ape[which(ape$X.genes.....genes == i),3] <- mean(sexual_vector)
  ape[which(ape$X.genes.....genes == i),4] <- mean(asexual_vector)
  
  asexual <- fasta2DNAbin(paste0("c1_asexual_", i, "_muscle.fas"))
  asexual_vector1 <- as.numeric(na.omit(dnds(asexual)))
  ape[which(ape$X.genes.....genes == i),5] <- wilcox.test(sexual_vector, asexual_vector1)[3]
  ape[which(ape$X.genes.....genes == i),6] <- mean(asexual_vector1)
  
  asexual <- fasta2DNAbin(paste0("c2_asexual_", i, "_muscle.fas"))
  asexual_vector2 <- as.numeric(na.omit(dnds(asexual)))
  ape[which(ape$X.genes.....genes == i),7] <- wilcox.test(sexual_vector, asexual_vector2)[3]
  ape[which(ape$X.genes.....genes == i),8] <- mean(asexual_vector2)
  
  sumasexual <- as.numeric(c(asexual_vector1, asexual_vector2))
  ape[which(ape$X.genes.....genes == i),9] <- wilcox.test(sexual_vector, sumasexual)[3]
  ape[which(ape$X.genes.....genes == i),10] <- mean(sumasexual)
}
ape

colnames(ape) <- c("gene", "p-value_vs_all_asexual","ds.dn_mean_sexual","ds.dn_mean_all_asexual","p-value_vs_Core_asexual","ds.dn_mean_Core_asexual","p-value_vs_Prasino_asexual", "ds.dn_mean_Prasino_asexual", "p-value_vs_Core+Prasino_asexual", 
                   "ds.dn_mean_Core+Prasino_asexual")
ape$dn.ds_mean_Prasino_asexual[c(which(ape$gene == "hop1"),which(ape$gene == "msh4"))] <- "No"
ape$`p-value_vs_Prasino_asexual`[c(which(ape$gene == "hop1"),which(ape$gene == "msh4"))] <- "No"

###Here you go. A table with a summary stats on SNAP output
print(ape)
###save wherever you wish
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/Phycocosm_annotation/")
write.csv(ape, "ape_values&stats.csv")
##################
#now, the table requires a bit of reshaping - for more convenient plotting
##################
data <- read.csv("reshaped_ape_values&stats.csv", sep = ",")
str(data)











