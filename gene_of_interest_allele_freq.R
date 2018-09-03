################################################
####  The the frequency at which variants   ####
####  were detected by exome sequencing     ####
####  and some population genetics data     ####
################################################

################################################
####   *Caveat*                             ####
####    This file is specific to            ####
####    my gene of interest.                ####
####    For now, I replace the gene name    ####
####    with the holder:                    ####
####    Gene of Interest ("GoI")            ####
################################################

################################################
####   *Caveat*                             ####
####    In one or two places, I also         ####
####    replace the variand ID (rs number   ####
####    or protein position etc.) with      ####
####    the holders "variant1, variant2,    ####
####    variant3, etc.                      ####
################################################

################################################
####   *Caveat*                             ####
####    lastly, I use the column headers    ####
####    of "SC" and "VNP" which are         ####
####    specific to this project.           ####
################################################

#getwd()
#setwd("...")

################################################
#### Load the tools                         ####
################################################
library(ggplot2)
library(tidyr)
library(dplyr)
require(gridExtra)

################################################
####     Section 1                          ####
################################################
####      Exome quality -                   ####
####      I don't have fastq yet            ####
####      to do real QC. However,           ####
####      suspect index swaps               ####
################################################
####      First look at the basic data      ####
################################################
####      Import the datasets               ####
################################################

c <- read.csv(
  "GoI_allele_freq.csv", 
  stringsAsFactors = TRUE,
  header = TRUE)

################################################
####      gather data to long form          ####
################################################
colnames(c)
#### Gather to long form by specifying the range
#### of columns to condence
cf <- gather(c,
             "variant1":"variantN",
             key = "variant", value = "callfreq")

################################################
####  Select rows that contain a string ... ####
################################################
####  SC samples                            ####
####  5 samples needed name modification    ####
sc <- dplyr::filter(cf, grepl("SC",sample))

####  VNP samples                           ####
vnp <- dplyr::filter(cf, grepl("VNP",sample))

####     Everyone else                      ####
other <- dplyr::filter(cf, !grepl("SC|VNP",sample))

################################################
####  Here come the plots                   ####
################################################

p1 <- ggplot(sc, aes(y = callfreq, x = variant, group = sample, colour = sample))+ 
  geom_point() + scale_x_discrete()+
  geom_line(alpha=0.2) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 40, hjust = 1))+
  labs(x="SC samples",
       y="Frequency per individual")

p2 <- ggplot(vnp, aes(y = callfreq, x = variant, group = sample, colour = sample))+ 
  geom_point() + scale_x_discrete()+
  geom_line(alpha=0.2) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 40, hjust = 1))+
  labs(x="VNP samples",
       y="Frequency per individual")

p3 <- ggplot(other, aes(y = callfreq, x = variant, group = sample, colour = sample))+ 
  geom_point() + scale_x_discrete()+
  geom_line(alpha=0.2) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 40, hjust = 1))+
  labs(x="BNvs392 set",
       y="Frequency per individual")

grid.arrange(p1, p2, p3, ncol=1, top = "GoI variants per sample")

################################################
#### Import datasets                        ####
################################################

frequencies <- read.csv(
  "GoI_allele_freq_alt.csv", 
  stringsAsFactors = TRUE,
  header = TRUE)

################################################
#### Print col names                        ####
####  and gather data to long form          ####
################################################

colnames(frequencies)
df <- gather(frequencies, 
"variant1", "variant2", ... , "variantN",
key = "variant", value = "callfreq")

################################################
####  Select rows that contain a string ... ####
################################################
  ####  SC                                  ####
  sc <- dplyr::filter(df, grepl("SC",sample))

  ####  VNP   ####
  vnp <- dplyr::filter(df, !grepl("SC",sample))

################################################
####  The the frequency at which variants   ####
####  were detected by exome sequencing     ####
####  different plots arranged together     ####
################################################
colorder <- c("variant1", 
              "variant2", 
              ... , 
              "variantN")

p1 <- ggplot(sc, aes(y = callfreq, x = variant, group = sample, colour = sample))+ 
  geom_point()+
  geom_line(alpha=0.2) + scale_x_discrete()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 0.8))+
scale_x_discrete(limits=colorder,labels=c(
  "variant1", 
  "variant2", 
  ... , 
  "variantN"))+
  labs(x="Protein position",
       y="Frequency per individual")

p2 <- ggplot(vnp, aes(y = callfreq, x = variant, group = sample, colour = sample))+ 
  geom_point()+
  geom_line(alpha=0.2) + scale_x_discrete()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 0.8))+
  scale_x_discrete(limits=colorder,labels=c(
  "variant1", 
  "variant2", 
  ... , 
  "variantN"))+
  labs(x="Protein position",
       y="Frequency per individual")

grid.arrange(p1, p2, ncol=1, top = "GoI variants per sample")


################################################
####    Section 2                           ####
################################################
####    Population genetics                 ####
################################################
####    Import datasets                     ####
################################################
#### Updated the c data with gnomad ref (d) ####
################################################

d <- read.csv(
  "GoI_combined_freq_coordinates.csv", 
  stringsAsFactors = TRUE,
  header = TRUE)

g <- read.csv(
  "GoI_gnomad_including_novel.csv", 
  stringsAsFactors = TRUE,
  header = TRUE)

################################################
####  Gnomad genomic plot                   ####
################################################

ggplot(g, aes(
  y = (Allele_Frequency), 
  x = Variant, colour = Number_of_Homozygotes))+ 
  scale_colour_gradient(low="blue", high="red")+
  geom_point(size=2)+
  scale_x_discrete()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 0.7))+
  labs(x="Protein position",
       y="Frequency in population")+
  geom_text(data=subset(g, Allele_Frequency > 0.02), 
            colour="black", 
            alpha=0.8,
            aes(Variant, Allele_Frequency, label=Consequence), 
            position = position_dodge2(width = 1), 
            vjust = 1.5, size=3)+
  theme(panel.background = element_rect("#F7F7F7"))

################################################
####  Gnomad coding regions plot            ####
################################################
####  remove rows with no functional consequence
gc <- g[!(is.na(g$Consequence) | g$Consequence==""), ]
gc0.1 <- subset(gc, Allele_Frequency<0.1)

ggplot(gc0.1, aes(
  y = (Allele_Frequency), 
  x = Variant, colour = Number_of_Homozygotes))+ 
  scale_colour_gradient(low="blue", high="red")+
  geom_point(size=2)+
  scale_x_discrete()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 0.7))+
  labs(x="Protein position",
       y="Frequency in population")+
  geom_text(data=subset(gc0.1, Allele_Frequency > 0.005), 
            colour="black", 
            alpha=0.8,
            aes(Variant, Allele_Frequency, label=Consequence), 
            position = position_dodge2(width = 1), 
            vjust = -1.5, size=4)+
  theme(panel.background = element_rect("#F7F7F7"))


################################################
####    Section 3                           ####
################################################
####  Check the frequency of our candidates ####
################################################
####  Merge the d and g data                ####
dg <- merge(d, g, by = "Variant", all = TRUE)

################################################
####  Select rows that contain a string     ####
####  those flgged as "candidate"           ####
################################################

dgc <- dplyr::filter(dg, grepl("candidate",Flag))

################################################
####    Plot candidates vs gnomad freq      ####
################################################

g4 <- 
  ggplot(dgc, aes(
  y = (Allele_Frequency), 
  x = Variant, 
  colour = Allele_Frequency))+ 
  scale_colour_gradient(low="blue", high="red")+
  geom_point(size=2)+
    scale_x_discrete(
                     labels=(dgc$Consequence))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
  labs(x="Protein position",
       y="Frequency in population")+
  geom_text(data=subset(dgc, Allele_Frequency > 0.01),
            colour="black", 
            alpha=0.8,
            aes(Variant, Allele_Frequency, label=Consequence), 
            position = position_dodge2(width = 1), 
            vjust = -1, size=3)+
  theme(panel.background = element_rect("#F7F7F7"))+
        ggtitle("Candidate variant frequency in GnomAD")


################################################
####  The the frequency at which variants   ####
####  were detected by exome sequencing     ####
####  and the frequencies in GnomAD         ####
################################################
####    Gather to long form                 ####
################################################
gdgc <- gather(dgc,
             "Sample1":"Sample16",
             key = "ID", value = "callfreq")

####    Super controllers                   ####
dgc_sc <- dplyr::filter(gdgc, grepl("SC",ID))

####    VNP               ####
dgc_vnp <- dplyr::filter(gdgc, grepl("VNP",ID))

####    Everyone else                       ####
dgc_other <- dplyr::filter(gdgc, !grepl("SC|VNP",ID))

################################################
####    Plot all 3 and combine              ####
################################################

g1 <- 
  ggplot(dgc_sc, aes(group=ID,
  y = (callfreq), 
  x = Variant, 
  colour = ID))+ 
  geom_point(size=2)+
  geom_line(alpha=0.2)+
  scale_x_discrete(
    labels=(dgc$Consequence))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title.x=element_blank())+ #this removes x-axis title
  labs(#x="Protein position",
       y="Frequency per sample")+
  theme(panel.background = element_rect("#F7F7F7"))+
  ggtitle("SC")

g2 <- 
  ggplot(dgc_vnp, aes(group=ID,
                   y = (callfreq), 
                   x = Variant,
                   colour = ID))+ 
  geom_point(size=2)+
  geom_line(alpha=0.2)+
  scale_x_discrete(
    labels=(dgc$Consequence))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title.x=element_blank())+ #this removes x-axis title)+
  labs(#x="Protein position",
       y="Frequency per sample")+
  theme(panel.background = element_rect("#F7F7F7"))+
  ggtitle("VNP")

g3 <- 
  ggplot(dgc_other, aes(group=ID,
                   y = (callfreq), 
                   x = Variant, 
                   colour = ID))+ 
  geom_point(size=2)+
  geom_line(alpha=0.2)+
  scale_x_discrete(
    labels=(dgc$Consequence))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title.x=element_blank())+ #this removes x-axis title)+
  labs(#x="Protein position",
       y="Frequency per sample")+
  theme(panel.background = element_rect("#F7F7F7"))+
  ggtitle("controls")

grid.arrange(g1, g2, g3, g4, ncol=1, top = "GoI variants per sample \n and allele freq in GnomAD")

####  End ####

