##R code to plot a region of interest with a gggenes plot underneath

# a helper function to plot p-values across the genome (not geo mean pval!)
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
    geom_point(aes(x=window_mb, y=!!y_var_sym, fill = pop_1_2), 
               shape=21, size=2, color="black", stroke=0.25) + 
    #scale_fill_discrete()
    theme_bw() +
    xlab("Position in genome (Mb)") +
    ylab(ylabel) +
    ggtitle(title_text, subtitle = today())
  
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
  filename = paste0("fisher_pval_chr_", chr, "_", start_mb, "-", end_mb, ".pdf")
  ggsave(filename, units="in", width=10, height=7.5)
  
  # return the plot
  p
}

#######Plot of region of interest with gggenes figure at the bottom showing a gene map


###ROI p-vals
exp<-plot_pvals(above_cutoff2, chr = "2L", start=21500000, end=21545000+1e3,
               ylabel="Geometric mean of p-values for SNPs in window")

###gggenes gene map
genes <- read.delim ("./myo.txt", sep="\t", header=T)
###gene file reads weird so had to specify that it's one row with multiple columns of data with headers
genes[1,6]=1


exp2<-
  ggplot(genes, aes(xmin=start,xmax=end, y = molecule, fill= gene, forward=orientation, label=gene, show.legend=FALSE))+
  geom_gene_arrow(show.legend = FALSE,
                  arrowhead_height = grid::unit(12, "mm"),
                  arrowhead_width = grid::unit(6, "mm"),
                  arrow_body_height = grid::unit(6, "mm")
  )+
  geom_gene_label(height = grid::unit(6, "mm"), grow = TRUE)+theme_genes()+labs(x='Position on genome', y='', title='')+xlim(500000, 2540000)

exp / exp2 + plot_layout(heights=c(1.15, .2))