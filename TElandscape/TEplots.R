### R script to create bar plots to represent TE superfamily compositions and the age distribution

library(ggplot2)
library(dplyr)  #has the merge function
library(tidyverse) #has the reduce function
library(reshape2) #This library has the melt function
library(readxl)
library(colorspace)
library(viridis)
library(hrbrthemes)
library(ggpubr) #for ggarrange

TECA <- read.csv('~/Documents/TEfasta/CATElandscape.txt', header =  TRUE, sep = "\t")
TErio <- read.csv('~/Documents/TEfasta/rioTElandscape.txt', header =  TRUE, sep = "\t")
TEleoti <- read.csv('~/Documents/TEfasta/leotiTElandscape.txt', header =  TRUE, sep = "\t")
TEpi229 <- read.csv('~/Documents/TEfasta/pi229841TElandscape.txt', header =  TRUE, sep = "\t")
TEpi297 <- read.csv('~/Documents/TEfasta/pi297155TElandscape.txt', header =  TRUE, sep = "\t")
TEpi300 <- read.csv('~/Documents/TEfasta/pi300119TElandscape.txt', header =  TRUE, sep = "\t")
TEpi329 <- read.csv('~/Documents/TEfasta/pi329311TElandscape.txt', header =  TRUE, sep = "\t")
TEpi506 <- read.csv('~/Documents/TEfasta/pi506069TElandscape.txt', header =  TRUE, sep = "\t")
TEpi510 <- read.csv('~/Documents/TEfasta/pi510757TElandscape.txt', header =  TRUE, sep = "\t")
TEpi655 <- read.csv('~/Documents/TEfasta/pi655972TElandscape.txt', header =  TRUE, sep = "\t")
TEgrassl <- read.csv('~/Documents/TEfasta/grasslTElandscape.txt', header = TRUE, sep = "\t")

geno <- list(TEpi655,TEpi297,TEpi229,TEleoti,TErio,TECA,TEpi510,TEpi506,TEpi329,TEgrassl)

mergedinfo <- geno %>% reduce(inner_join, by="Tefam")

datamelt <- melt(mergedinfo, id="Tefam")

### creating the bar plot
bplot <- ggplot(datamelt, aes(fill=Tefam, x=variable, y=value)) + geom_bar(position="fill",stat="identity") + theme_classic() + xlab("") + ylab("Percentage of TEs") + coord_flip() +  scale_fill_manual(values=c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")) +
 theme(text=element_text(size=12,  family="sans")) +
  labs(fill="TE families")
  
### creating the age distribution histogram for all 10 genotypes
TE_CA <- read.delim('CA_TEpos.txt')
TE_leoti <- read.delim('leoti_TEpos.txt')
TE_pi229841 <- read.delim('pi229841_TEpos.txt')
TE_pi297155 <- read.delim('pi297155_TEpos.txt')
TE_pi329311 <- read.delim('pi329311_TEpos.txt')
TE_pi506069 <- read.delim('pi506069_TEpos.txt')
TE_pi510757 <- read.delim('pi510757_TEpos.txt')
TE_pi655972 <- read.delim('pi655972_TEpos.txt')
TE_rio <- read.delim('rio_TEpos.txt')

TE_CA$genotype <- 'Chinese Amber'
TE_leoti$genotype <- 'Leoti'
TE_pi229841$genotype <- 'pi229841'
TE_pi297155$genotype <- 'pi297155'
TE_pi329311$genotype <- 'pi329311'
TE_pi506069$genotype <- 'pi506069'
TE_pi510757$genotype <- 'pi510757'
TE_pi655972$genotype <- 'pi655972'
TE_rio$genotype <- 'Rio'

master_data <- rbind(TE_CA,TE_leoti,TE_pi229841,TE_pi297155,TE_pi329311,TE_pi506069,TE_pi510757,TE_pi655972,TE_rio)

hplot <- ggplot(data=master_data, aes(x=Age_MY,fill = genotype)) + geom_histogram(bins = 50) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") + 
scale_fill_manual(values=c(rep(c("#F0E442","#117733"),5))) + scale_x_continuous(breaks=seq(0,8,by=.5),limits=c(0,8)) +
xlab("Age in Million years")

### Plotting both in the same figure
finalplot <- ggarrange(hplot,bplot, ncol=2,labels="AUTO")
ggsave('~/Desktop/TEplots.png',finalplot, width = 10, height = 4, units = "in", dpi = 300)
