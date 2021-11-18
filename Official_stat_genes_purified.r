library("dplyr")
library("pspearman")
library("ggplot2")
library("RColorBrewer")
setwd("/home/zavadska/Documents/biology/all_chlorophyta_its2/")
data <- read.csv("Official_gene_annotated.csv", sep = ",")
str(data)

#Normality test
with(data, shapiro.test(as.numeric(number.of.genes)[sexual.reproduction. == "N"]))
with(data, shapiro.test(as.numeric(number.of.genes)[sexual.reproduction. == "Y"]))

#Wilcox test
wilcox.test(as.numeric(data$number.of.genes[data$sexual.reproduction. == "Y"]), as.numeric(data$number.of.genes[data$sexual.reproduction. == "N"]))

#Test for the homogeneity of variations
library("lawstat")
levene.test(data$number.of.genes, group = data$sexual.reproduction.)

#Violin plot  
data1 <- within(data,{
  sexual.reproduction.factor <- NA
  sexual.reproduction.factor[sexual.reproduction. == "N"] <- "Unknown"
  sexual.reproduction.factor[sexual.reproduction. == "Y"] <- "Present" 
  })


violin_plot <- ggplot(data1, aes(x=factor(sexual.reproduction.factor,levels = c("Unknown","Present")), y=as.numeric(number.of.genes), fill=sexual.reproduction., color=sexual.reproduction.)) + 
  geom_violin(alpha = 0.5)+
  xlab("Sexual reproduction") + ylab("Number of MGHs per genome") +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text( size=15))+
  labs(fill = "Sexual reproduction")+
  geom_boxplot(width=0.1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)+theme_bw(base_size = 15)+ 
  theme(legend.position = "none")+scale_color_manual(values = c("#1F78B4","#1F78B4"))+scale_fill_manual(values = c("#A6CEE3","#FB9A99"))

ggsave("MGHs_violin.png",violin_plot,width = 3.5, height = 3.5)


#Boxplots (not included in the article)
ggplot(data, aes(x=sexual.reproduction., y=as.numeric(number.of.genes), fill=sexual.reproduction.)) + 
  geom_boxplot()+
  xlab("Sexual reproduction") + ylab("Number of meiotic gene homologs") +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text( size=15), legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(fill = "Sexual reproduction")+
  scale_fill_discrete(name = "Sexual reproduction", labels = c("Absent", "Present"))