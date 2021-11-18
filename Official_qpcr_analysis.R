library(pcr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
setwd('/home/zavadska/Documents/biology/all_chlorophyta_its2/Konstanz-wetlab part/concentrated')

######pcr package estimates; formulas/explanation in Supplement#####

###put the directory with CT techical replicate means below.
ct1 <- read.csv("Official_CONC_forCTvalue_analysis_means.csv", sep = ",")

###making a list of genes, conditions and subsetting a dataframe to get rid of C.reinhardtii CT values
group_var <- rep(c("Nitrate deficiency","Standard culture","Darkness-light","Darkness"), each = 3)
ct1 <- ct1[1:12,]
genes <- c("mer3","hop1","hop2","mnd1")

# calculate all values and errors in one step
## mode == 'separate_tube' default
resi <- data.frame()
for(i in genes){
  resa <- pcr_analyze(ct1[,c(which(colnames(ct1) == i), which(colnames(ct1) == paste0('L5.',i)))],
                     group_var = group_var,
                     reference_gene = paste0('L5.',i),
                     reference_group = 'Standard culture')
  print(resa)
  resi <- rbind(resi, resa)
}
resi
#########for significance test by lm(multiple conditions)#########
resstati <- data.frame()
for(i in genes){
  resstata <- pcr_test(ct1[,c(which(colnames(ct1) == i), which(colnames(ct1) == paste0('L5.',i)))],
                     group_var = group_var,
                     reference_gene = paste0('L5.',i),
                     reference_group = 'Standard culture',
                  test = 'lm')
  print(resstata)
  resstati <- rbind(resstati, resstata)
}

resstati 

#now, the "resstati" and "resi" tables require a bit of reshaping - for more convenient plotting. Those tables, reshaped and combined together, are found in "ddCT_calc_res_publish.csv" - see "plotting"

###plotting
res <- read.csv("Official_ddCT_calc_res_publish.csv", sep = ",")
str(res)

res1 <- within(res,{p_value_factor <- NA
p_value_factor[p_value <= 0.05] <- "Significant"
p_value_factor[p_value >= 0.05] <- "Not significant" 
p_value_factor[p_value == "No"] <- "Control"
})

res1$relative_expression[12] <- NA
res1

qPCR <- ggplot(data=res1, aes(x=factor(group, levels = c("Standard culture","Nitrate deficiency","Darkness","Darkness-light")), y=relative_expression, group=gene, color = gene, shape = p_value_factor)) +
  yscale("log2", .format = TRUE)+ 
  geom_line(position = position_dodge(width=0.5))+
  geom_pointrange(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)),size=0.8, position = position_dodge(width=0.5))+
  #geom_errorbar(aes(ymin=as.numeric(relative_expression-error), ymax=as.numeric(relative_expression+error)), size=0.8, position = position_dodge(width=0.5))+
  xlab("Culture conditions") + ylab("Relative mRNA expression") +
  labs( color = "Gene", shape = "Difference control \nvs test condition")+ scale_color_brewer(palette = "Paired")+theme_bw(base_size = 15)+
  scale_shape_manual(values=c(7:10))

ggsave("qPCR_results.png", qPCR,width = 10, height = 7.5)

