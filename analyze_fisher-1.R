library(tidyverse)
library(patchwork)
library(gggenes)
library(ggfittext)
# how big are windows?  This must match parameter used when running popoolation2
window_size <- 25000 
window_size_mb <- window_size / 1e6 

fet_wild <- read.delim ("./wild2.fet", sep="\t", header=F)
fet_foco17_1<- read.delim("./Foco17_1.fet", sep= "\t", header=F)
fet_foco17_2<- read.delim("./Foco17_2.fet", sep= "\t", header=F)
fet_foco23_1<- read.delim("./Foco23_1.fet", sep= "\t", header=F)
fet_foco23_2<- read.delim("./Foco23_2.fet", sep= "\t", header=F)



handle_input_fet <- function(df_to_handle, sample_set="wild") {

  # create column names (popoolation2 doesn't output column names)
  extra_col_names = paste0("X", seq(1,ncol(df_to_handle) - 5))
  colnames(df_to_handle) <- c("chrom", "window", "n_snps", "frac_sufficient_cov", "avg_min_coverage", extra_col_names) 
  
  df <- df_to_handle %>% pivot_longer(cols=starts_with("X"), names_to="unhelpful_name", values_to="p_p")  %>% select(-unhelpful_name)
  
  # chr2L	5000	1	0.002	59.3	1:2=2.29269731	1:3=2.91987693	1:4=9.43200886	1:5=33.37451954	1:6=0.94702188	1:7=1.51466859	1:8=0.17907030	2:3=0.50657442	2:4=7.84675129	2:5=70.85642166	2:6=1.11831056	2:7=0.53391643	2:8=4.34439480	3:4=5.34500299	3:5=43.77295891	3:6=1.88328402	3:7=1.24666426	3:8=5.12927387	4:5=5.51537397	4:6=9.84815294	4:7=9.47870196	4:8=14.57806963	5:6=59.53783747	5:7=75.64739283	5:8=83.47553973	6:7=0.38735115	6:8=1.09731604	7:8=2.61145212
  
  # parse metadata out of p-value value
  subparts <- str_match(df$p_p, "(\\d+):(\\d+)=([0-9\\.]+)")
  
  df$pop_1 <- as.character(subparts[,2])
  df$pop_2 <- as.character(subparts[,3])
  df$log_pval <- as.numeric(subparts[,4])
  
  # populations
if (sample_set == "wild") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "A2",
                           "2" ~ "Neg")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "A2",
                           "2" ~ "Neg")
    
  } else if (sample_set == "foco") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Negative2",
                           "3" ~ "Foco17Positive",
                           "4" ~ "Foco23Negative1",
                           "5" ~ "Foco23Negative2",
                           "6" ~ "Foco23Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Negative2",
                           "3" ~ "Foco17Positive",
                           "4" ~ "Foco23Negative1",
                           "5" ~ "Foco23Negative2",
                           "6" ~ "Foco23Positive")
  } else if (sample_set == "foco17_1") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Positive")
    
  } else if (sample_set == "foco17_2") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco17Negative2",
                           "2" ~ "Foco17Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco17Negative2",
                           "2" ~ "Foco17Positive")
  } else if (sample_set == "foco23_1") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco23Negative1",
                           "2" ~ "Foco23Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco23Negative1",
                           "2" ~ "Foco23Positive")
    
  } else if (sample_set == "foco23_2") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco23Negative2",
                           "2" ~ "Foco23Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco23Negative2",
                           "2" ~ "Foco23Positive")
  }
  
  df <- df %>% mutate(pop_1_2 = paste0(pop_1, "-", pop_2))
  
  # keep track of position in megabases
  df <- df %>% mutate(window_mb = window / 1e6)
  
  df <- df %>% mutate(window_start = window - window_size / 2,
                      window_end   = window + window_size / 2)
  
  # what does distribution of fraction sufficient coverage look like?
  # ggplot(df) + 
    # geom_histogram(aes(x=frac_sufficient_cov), bins=100)  +
    # theme_bw()
  
  # filter out windows with insufficient coverage
  df <- df %>% filter(frac_sufficient_cov > 0.05)
  
  # ggplot(df) + 
    # geom_histogram(aes(x=frac_sufficient_cov), bins=100)  +
    # theme_bw()
  
  # only keep big chromosomes
  df <- filter(df, str_detect(chrom, "[23X]"))
  
  # calculate geometric mean of p-values
  # which is the nth root of the product of p-values 
  # (the product is what popoolation2 reports by default)
  # n = the # of SNPs 
  df <- df %>% mutate(geo_mean_pval = log_pval ^ (1/n_snps),
                     mean_pval = log_pval / n_snps)
  
  # return the df
  df
}


df_wild <- handle_input_fet(fet_wild, sample_set="wild")

df_17_1 <- handle_input_fet(fet_foco17_1, sample_set="foco17_1")

df_17_2 <- handle_input_fet(fet_foco17_2, sample_set="foco17_2")

df_23_1 <- handle_input_fet(fet_foco23_1, sample_set="foco23_1")

df_23_2 <- handle_input_fet(fet_foco23_2, sample_set="foco23_2")

######Add new column to say which pool it is
df_17_1$pool<- "Foco17 Pool 1"
df_17_2$pool<-"Foco17 Pool 2"
df_23_1$pool<- "Foco23 Pool 1"
df_23_2$pool<-"Foco23 Pool 2"
df_wild$pool<-"Wild"

#####combine all files

df_combined <- rbind(df_17_1, df_17_2, df_23_1, df_23_2, df_wild)

combined_filter<- filter(df_combined, 
                         pop_1_2 == "A2-Neg" | 
                           pop_1_2 == "Foco17Negative1-Foco17Positive" | 
                          # pop_1_2 == "Foco17Negative1-Foco17Negative2" | 
                           pop_1_2 == "Foco17Negative2-Foco17Positive" |
                           pop_1_2 == "Foco23Negative1-Foco23Positive" |
                           pop_1_2 == "Foco23Negative2-Foco23Positive" 
                          # pop_1_2 == "Foco23Negative1-Foco23Negative2"
                      )
combined_filter_sort<- combined_filter[order(-combined_filter$log_pval),]

foco23<-filter(df_combined,
                 pop_1_2 == "Foco23Negative1-Foco23Positive" |
                 pop_1_2 == "Foco23Negative2-Foco23Positive")
foco17<-filter(df_combined,
               pop_1_2 == "Foco17Negative1-Foco17Positive" | 
                 pop_1_2 == "Foco17Negative2-Foco17Positive")

#####cutoff by pval#######
p_val_cutoffs_window<- df_combined %>%
  group_by(pop_1_2)%>%
  reframe(log_pval=quantile(log_pval,c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)), q=c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) %>%
  filter(q==0.95)%>% 
  rename(log_pval_cutoff=log_pval)

df_combined<- left_join(df_combined,p_val_cutoffs_window)

above_cutoff_window<- filter(df_combined, log_pval > log_pval_cutoff)

above_cutoff_window<-above_cutoff_window%>%
  group_by(chrom, window)%>%
  mutate(n_groups= n())



###save to file
write.table(above_cutoff, file= "above_cutoff.txt", sep=",") 

#Filter dataframes so that they are sorted highest p-value to lowest
df_foco_sort<- df_foco[order(-df_foco$log_pval),]

ggplot(filter(df_combined, 
              pop_1_2 == "A2-Neg" | 
              pop_1_2 == "Foco17Negative1-Foco17Positive" | 
              pop_1_2 == "Foco17Negative1-Foco17Negative2" | 
                pop_1_2 == "Foco17Negative2-Foco17Positive" |
         #       pop_1_2 == "Foco23Negative1-Foco23Positive" |
                pop_1_2 == "Foco23Negative2-Foco23Positive"
          )) +
  geom_histogram(aes(x=log_pval), bins=100) +
  facet_wrap(~pop_1_2, ncol = 1) +
  scale_y_log10()



# a helper function to plot p-values across the genome
plot_pvals <- function(df_to_plot, 
                       chr = NA, 
                       start = NA, 
                       end = NA, 
                       p1 = NA, 
                       p2 = NA, 
                       title=NA,
                       facet_1 = "pop_1_2",
                       facet_2 = NA,
                       y_var = "log_pval",
                       ylabel = "product(-log10(pvalue))") { 
  
  # do successive filtering based on parameters 
  if (!is.na(chr)) {
    df_to_plot <- filter(df_to_plot, chrom==chr)
  }
  
  # must specify start & end
  if (!is.na(end) & is.na(start)) {
    message ("Error: must specify start and end position")
    quit(status = 1)
  }
  
  # plot part of one or more chrom
  if (!is.na(start)) {
    if (is.na(end)) { message ("Error: must specify start and end position")
      quit(status = 1)
    }
    df_to_plot <- filter(df_to_plot, window > start & window < end)
  }
  
  # plot one or more populations
  if (!is.na(p1)) {
    df_to_plot <- filter(df_to_plot, pop_1 == p1)
  }
    
  if (!is.na(p2)) {
    df_to_plot <- filter(df_to_plot, pop_2 == p2)
  }
  
  # region of interest labels
  start_mb <- sprintf("%0.2f", start/1e6)
  end_mb   <- sprintf("%0.2f", end/1e6)
  
  # can define title or use default
  if (!is.na(title)) {
    title_text = title
  } else {
    title_text =  paste0("Fisher p-vals: chr. ", chr, " : ", start_mb, "-", end_mb, "Mbp")
  }
  
  y_var_sym <- sym(y_var)
  
  # the plot
  p <- ggplot(df_to_plot) +
    # geom_point(aes(x=window_mb, y=log_pval), fill="slateblue", 
    geom_point(aes(x=window, y=!!y_var_sym, fill = pool), 
               shape=21, size=2, color="black", stroke=0.25) + 
    #scale_fill_discrete()
    theme_bw() +
    xlab("Position in genome (Mb)") +
    ylab(ylabel) + facet_wrap(chrom~pool, nrow=1)
    ggtitle(title_text, subtitle = today())
  
  # # do 1D or 2D faceting
  # # see this section on using this .data [[]] syntax
  # # https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html
 #  if (is.na(facet_2)) {
  #   p <- p + facet_wrap( vars(.data[[facet_1]] ) , scales = "free_x")
  # } else {
   #  p <- p + facet_grid( vars(.data[[facet_1]]), vars(.data[[facet_2]]), scales = "free_x" )
  # }
   
  
  # save a pdf of plot with dynamic filename
 # filename = paste0("fisher_pval_chr_", chr, "_", start_mb, "-", end_mb, ".pdf")
  #ggsave(filename, units="in", width=10, height=7.5)
  
  # return the plot
  p
}


# -------------------------------------------------
# plot some particular ROI
# from wild A2/Neg GWAS
#
#  Note: plotting 100k bases on either side of ROI 
# -------------------------------------------------

# General test code to try out different regions

plot_pvals(above_cutoff, chr = "3R", start=565000, end=2020000+1e5,
          ylabel="SNP p-values")#+theme(axis.text.x=element_blank())




# random region
plot_pvals(combined_filter_sort, chr = "2R", start = 12020000, end=12050000)


# A function to narrow in on interesting (low p-val) peaks
filter_peaks_one_pair <- function(pval_df, p1, p2, pval_cutoff=500) {
  
  p1_2_name = paste0(p1, "-", p2)
  
  ggplot(filter(pval_df, pop_1 == p1 & pop_2 == p2)) +
    geom_point(aes(x=window_mb, y=log_pval, fill=chrom), shape=21, size=1.5, color="black", stroke=0.25)+
    theme_bw() +
    xlab("Chromosome position, Mb") +
    ylab("-log10(p-value)") +
    facet_wrap(~chrom, scales="free_x") +
    theme(legend.position = "none")
  
  ggsave(paste0(p1_2_name, "_pval.pdf"), width=10, height=7.5, units="in")
  
  # filter interesting peaks
  peaks <- filter(pval_df, pop_1 == p1 & pop_2 == p2 & log_pval > pval_cutoff)
  print(peaks, n=50)
  peaks
}

# all peaks
all_peaks <- filter(df, log_pval > 900) %>%
  mutate(chr_win = paste0(chrom, window))


# lots of noise from GWAS (F populations) vs. negative comparisons
# filter them out but keep A2-Neg
# all_peaks <- all_peaks %>% filter(pop_2 != "Neg" | pop_1_2 == "A2-Neg")

ggplot(all_peaks) +
  geom_tile(aes(x=window_mb, y=pop_1_2, fill=log_pval)) +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~chrom, nrow = 1)
  
ggsave("all_peaks_tiles.pdf", width=10, height=7.5, units="in")

library(pheatmap)

# convert to a numeric matrix for pheatmap
all_peaks_mat <- all_peaks %>% select(chr_win, log_pval, pop_1_2)
all_peaks_mat <- all_peaks_mat %>% pivot_wider(id_cols = chr_win, names_from = "pop_1_2", values_from = log_pval)
row.names(all_peaks_mat) <- all_peaks_mat$chr_win
all_peaks_mat <- all_peaks_mat %>% select(-chr_win)
all_peaks_mat <- as.matrix(all_peaks_mat)
all_peaks_mat[is.na(all_peaks_mat)] <- 0

pheatmap(all_peaks_mat)

?pheatmap


d <- filter_peaks_one_pair(df, "A2", "Neg")
d <- filter_peaks_one_pair(df, "F2", "F5")

A2_neg_peaks_400 <-  filter(df, pop_1_2 == "A2-Neg" & log_pval > 400)
write.table(A2_neg_peaks_400, file="A2_neg_peaks_gt_400.txt", sep="\t", quote=F, row.names=F)
A2_neg_peaks_400


F2_F3_peaks_400 <-  filter(df, pop_1_2 == "F2-F3" & log_pval > 400)
print(F2_F3_peaks_400, n=50)


write.table(F2_F3_peaks_400, file="F2-F3_peaks_gt_400.txt", sep="\t", quote=F, row.names=F)

ggplot(filter(df, pop_1_2 == "F3-F6")) +
  geom_point(aes(x=window_mb, y=log_pval, fill=chrom), shape=21, size=1.5, color="black", stroke=0.25)+
  theme_bw() +
  xlab("Chromosome position, Mb") +
  ylab("-log10(p-value)") +
  facet_wrap(~chrom, scales="free_x") +
  theme(legend.position = "none",)

ggsave("F3-F6_pval.pdf", width=10, height=7.5, units="in")

F3_F6_peaks_400 <-  filter(df, pop_1_2 == "F3-F6" & log_pval > 400)
print(F3_F6_peaks_400, n=50)

write.table(F3_F6_peaks_400, file="F3-F6_peaks_gt_400.txt", sep="\t", quote=F, row.names=F)


# from popoolation2 manual
# https://code.google.com/archive/p/popoolation2/wikis/Manual.wiki
#
# column 1: reference contig
# column 2: position of the window (middle)
# column 3: SNPs identified in the window
# column 4: fraction of window having sufficient coverage
# column 5: average minimum coverage
# column >5: Fst values for all pairwise comparisons; For example "2:3=0.02" states that the Fst for comparing population 2 with population 3 (using the order of the synchronized file) is 0.2



###plot geo mean in one specific ROI
# a helper function to plot p-values across the genome
plot_geos <- function(df_to_plot, 
                       chr = NA, 
                       start = NA, 
                       end = NA, 
                       p1 = NA, 
                       p2 = NA, 
                       title=NA,
                       facet_1 = "pop_1_2",
                       facet_2 = NA,
                       y_var = "geo_mean_pval",
                       ylabel = "product(-log10(geo_mean_pval))") { 
  
  # do successive filtering based on parameters 
  if (!is.na(chr)) {
    df_to_plot <- filter(df_to_plot, chrom==chr)
  }
  
  # must specify start & end
  if (!is.na(end) & is.na(start)) {
    message ("Error: must specify start and end position")
    quit(status = 1)
  }
  
  # plot part of one or more chrom
  if (!is.na(start)) {
    if (is.na(end)) { message ("Error: must specify start and end position")
      quit(status = 1)
    }
    df_to_plot <- filter(df_to_plot, window > start & window < end)
  }
  
  # plot one or more populations
  if (!is.na(p1)) {
    df_to_plot <- filter(df_to_plot, pop_1 == p1)
  }
  
  if (!is.na(p2)) {
    df_to_plot <- filter(df_to_plot, pop_2 == p2)
  }
  
  # region of interest labels
  start_mb <- sprintf("%0.2f", start/1e6)
  end_mb   <- sprintf("%0.2f", end/1e6)
  
  # can define title or use default
  if (!is.na(title)) {
    title_text = title
  } else {
    title_text =  paste0("Geometric mean P-vals for specified region with corresponding gene map")
      #"Geo mean p-vals: chr. ", chr, " : ", start_mb, "-", end_mb, "Mbp")
  }
  
  y_var_sym <- sym(y_var)
  
  # the plot
  p <- ggplot(df_to_plot) +
    # geom_point(aes(x=window_mb, y=log_pval), fill="slateblue", 
    geom_point(aes(x=window_mb, y=!!y_var_sym, fill = pool), 
               shape=21, size=2, color="black", stroke=0.25, show.legend=FALSE) + 
    #scale_fill_discrete()
    theme_bw() +
    xlab('') +
    ylab(ylabel) +
    ggtitle(title_text)
  
  # # do 1D or 2D faceting
  # # see this section on using this .data [[]] syntax
  # # https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html
  # if (is.na(facet_2)) {
  #   p <- p + facet_wrap( vars(.data[[facet_1]] ) , scales = "free_x")
  # } else {
  #   p <- p + facet_grid( vars(.data[[facet_1]]), vars(.data[[facet_2]]), scales = "free_x" )
  # }
  # 
  
  # save a pdf of plot with dynamic filename
  filename = paste0("geo_pval_chr_", chr, "_", start_mb, "-", end_mb, ".pdf")
  ggsave(filename, units="in", width=10, height=7.5)
  
  # return the plot
  p
}


#######Individual Fisher's exact instead of windowed
fet_17_1_i<- read.delim("./17_1_indiv.fet", sep= "\t", header=F)
fet_17_2_i<- read.delim("./17_2_indiv.fet", sep= "\t", header=F)
fet_23_1_i<- read.delim("./23_1_indiv.fet", sep= "\t", header=F)
fet_23_2_i<- read.delim("./23_2_indiv.fet", sep= "\t", header=F)
fet_a2_i<- read.delim("./A2_indiv.fet", sep="\t", header =F)

df_a2_i <- handle_input_fet(fet_a2_i, sample_set="wild")
df_17_1_i <- handle_input_fet(fet_17_1_i, sample_set="foco17_1")
df_17_2_i <- handle_input_fet(fet_17_2_i, sample_set="foco17_2")
df_23_1_i <- handle_input_fet(fet_23_1_i, sample_set="foco23_1")
df_23_2_i <- handle_input_fet(fet_23_2_i, sample_set="foco23_2")
df_combined_i <- rbind(df_17_1_i, df_17_2_i, df_23_1_i, df_23_2_i, df_a2_i)



####cutoff by p-val
p_val_cutoffs_i<- df_combined_i %>%
  group_by(pop_1_2)%>%
  reframe(log_pval=quantile(log_pval,c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)), q=c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) %>%
  filter(q==0.95)%>% 
  rename(log_pval_cutoff=log_pval)

df_combined_i<- left_join(df_combined_i,p_val_cutoffs_i)

above_cutoff_i<- filter(df_combined_i, log_pval > log_pval_cutoff)

above_cutoff_i<-above_cutoff_i%>%
  group_by(chrom, window)%>%
  mutate(n_groups= n())


plot_geos(above_cutoff_i, chr = "3R", start=32057546-1e5, end=32040891+1e5,
               ylabel="Geometric mean of p-values for SNPs in window")#+theme(axis.text.x=element_blank())

