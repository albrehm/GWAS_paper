#######Manhattan Plot########
install.packages("qqman")
library(qqman)
library(ggplot2)


#####Filter by p-value to avoid clutter in final plot, here in 95th percentile
i_17 <- filter(df_17_1_i, mean_pval> 11.291584)
i2_17 <- filter(df_17_2_i, mean_pval> 7.864891)
i_23 <- filter(df_23_1_i, mean_pval> 5.246754)
i2_23 <- filter(df_23_2_i, mean_pval> 8.261640)
i_a2 <- filter(df_a2_i, mean_pval> 2.912523)

#Remove columns
qq_geo_17_1<-(i_17[-c(4:14)])
qq_geo_17_2<-(i2_17[-c(4:14)])
qq_geo_23_1<-(i_23[-c(4:14)])
qq_geo_23_2<-(i2_23[-c(4:14)])
qq_geo_a2<-(i_a2[-c(4:14)])

#rename variables so they match qqman requirements
qq_geo_17_1 <- qq_geo_17_1 %>%
  rename(
    CHR = chrom, BP = window, P = mean_pval, SNP = n_snps)

qq_geo_17_2 <- qq_geo_17_2 %>%
  rename(
    CHR = chrom, BP = window, P = mean_pval, SNP = n_snps)

qq_geo_23_1 <- qq_geo_23_1 %>%
  rename(
    CHR = chrom, BP = window, P = mean_pval, SNP = n_snps)

qq_geo_23_2 <- qq_geo_23_2 %>%
  rename(
    CHR = chrom, BP = window, P = mean_pval, SNP = n_snps)

qq_geo_a2 <- qq_geo_a2 %>%
  rename(
    CHR = chrom, BP = window, P = mean_pval, SNP = n_snps)

#rename chromosomes in dataframe. 2L=1, 2R=2, 3L=3, 3R=4, X=5
qq_geo_17_1<- qq_geo_17_1 %>% mutate(CHR = recode(CHR, '2L' ='chr2L', '2R' = 'chr2R', '3L' = 'chr3L', '3R' = 'chr3R', 'X' = 'chrX'))
qq_geo_17_2<- qq_geo_17_2 %>% mutate(CHR = recode(CHR, '2L' ='chr2L', '2R' = 'chr2R', '3L' = 'chr3L', '3R' = 'chr3R', 'X' = 'chrX'))
qq_geo_23_1<- qq_geo_23_1 %>% mutate(CHR = recode(CHR, '2L' ='chr2L', '2R' = 'chr2R', '3L' = 'chr3L', '3R' = 'chr3R', 'X' = 'chrX'))
qq_geo_23_2<- qq_geo_23_2 %>% mutate(CHR = recode(CHR, '2L' ='chr2L', '2R' = 'chr2R', '3L' = 'chr3L', '3R' = 'chr3R', 'X' = 'chrX'))
qq_geo_a2<- qq_geo_a2 %>% mutate(CHR = recode(CHR, '2L' ='chr2L', '2R' = 'chr2R', '3L' = 'chr3L', '3R' = 'chr3R', 'X' = 'chrX'))

#convert columns to compatible types
qq_geo_17_1<- transform(qq_geo_17_1, CHR= as.integer(CHR))
qq_geo_17_1<-transform(qq_geo_17_1, SNP= as.character(SNP))
qq_geo_17_1<-transform(qq_geo_17_1, P = as.numeric(P))
qq_geo_17_1<- transform(qq_geo_17_1, BP=as.integer(BP))


#remove p values of 0
qq_geo_17_1<-na.omit(qq_geo_17_1)
qq_geo_17_1<-filter(qq_geo_17_1, P > 0)

qq_geo_17_2<-na.omit(qq_geo_17_2)
qq_geo_17_2<-filter(qq_geo_17_2, P > 0)

qq_geo_23_1<-na.omit(qq_geo_23_1)
qq_geo_23_1<-filter(qq_geo_23_1, P > 0)

qq_geo_23_2<-na.omit(qq_geo_23_2)
qq_geo_23_2<-filter(qq_geo_23_2, P > 0)

qq_geo_a2<-na.omit(qq_geo_a2)
qq_geo_a2<-filter(qq_geo_a2, P > 0)

###qqplot
#qq(qq_geo$P)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("karyoploteR")
library(karyoploteR)
install.packages("GenomicRanges")
library(GenomicRanges)
library(dplyr)


qq_geo_17_1 <- qq_geo_17_1 %>% mutate(BP2 = BP)
qq_geo_17_2 <- qq_geo_17_2 %>% mutate(BP2 = BP)
qq_geo_23_1 <- qq_geo_23_1 %>% mutate(BP2 = BP)
qq_geo_23_2 <- qq_geo_23_2 %>% mutate(BP2 = BP)
qq_geo_a2 <- qq_geo_a2 %>% mutate(BP2 = BP)

qq_geo_17_1<- makeGRangesFromDataFrame(qq_geo_17_1, keep.extra.columns = TRUE, ignore.strand=TRUE, seqnames.field = c("CHR"), start.field = "BP",
                         end.field = "BP2")
qq_geo_17_2<- makeGRangesFromDataFrame(qq_geo_17_2, keep.extra.columns = TRUE, ignore.strand=TRUE, seqnames.field = c("CHR"), start.field = "BP",
                                       end.field = "BP2")
qq_geo_23_1<- makeGRangesFromDataFrame(qq_geo_23_1, keep.extra.columns = TRUE, ignore.strand=TRUE, seqnames.field = c("CHR"), start.field = "BP",
                                       end.field = "BP2")
qq_geo_23_2<- makeGRangesFromDataFrame(qq_geo_23_2, keep.extra.columns = TRUE, ignore.strand=TRUE, seqnames.field = c("CHR"), start.field = "BP",
                                       end.field = "BP2")
qq_geo_a2<- makeGRangesFromDataFrame(qq_geo_a2, keep.extra.columns = TRUE, ignore.strand=TRUE, seqnames.field = c("CHR"), start.field = "BP",
                                       end.field = "BP2")
seqlevelsStyle(qq_geo_17_1) <- "UCSC"
seqlevelsStyle(qq_geo_17_2) <- "UCSC"
seqlevelsStyle(qq_geo_23_1) <- "UCSC"
seqlevelsStyle(qq_geo_23_2) <- "UCSC"
seqlevelsStyle(qq_geo_a2) <- "UCSC"

kp <- plotKaryotype(genome= "dm3", plot.type=4)
kpAddChromosomeSeparators(kp)
kpAddMainTitle(kp, main = "Mahattan plot of -log10 p-values by group per chromosome")


####Plot them all together!#####
kpAddLabels(kp, labels = "P1R1", srt=90, pos=3, r0=autotrack(1,5))
kpPlotManhattan(kp, data = qq_geo_17_1, pval = qq_geo_17_1$P, logp=FALSE, r0=autotrack(1,5), points.col = c("cyan","cyan2"))
kpAddLabels(kp, labels = "P1R2", srt=90, pos=3, r0=autotrack(2,5))
kp <- kpPlotManhattan(kp, data = qq_geo_17_2, pval = qq_geo_17_2$P, logp=FALSE,  r0=autotrack(2,5), points.col=c("cyan3","cyan4"))
kpAddLabels(kp, labels = "P2R1", srt=90, pos=3, r0=autotrack(3,5))
kp <- kpPlotManhattan(kp, data = qq_geo_23_1, pval = qq_geo_23_1$P, logp=FALSE,  r0=autotrack(3,5), points.col= c("coral", "coral2"))
kpAddLabels(kp, labels = "P2R2", srt=90, pos=3, r0=autotrack(4,5))
kp <- kpPlotManhattan(kp, data = qq_geo_23_2, pval = qq_geo_23_2$P, logp=FALSE,  r0=autotrack(4,5), points.col=c("coral3","coral4"))
kpAddLabels(kp, labels = "P3R1", srt=90, pos=3, r0=autotrack(5,5))
kp <- kpPlotManhattan(kp, data = qq_geo_a2, pval = qq_geo_a2$P, logp=FALSE,  r0=autotrack(5,5), points.col=c("deeppink", "deeppink3"))






