# ------create cutoffs------
cutoffs_combined<- df_combined %>%
  group_by(pool)%>%
  reframe(log_p=quantile(log_p,c(0.75, 0.95, 0.99, 0.995, 0.999)), q=c(0.75, 0.95, 0.99, 0.995, 0.999)) %>%
  filter(q==0.995)%>% 
  rename(log_pval_cutoff=log_p)

df_combined_group<- left_join(df_combined,cutoffs_combined)

df_combined_group<- filter(df_combined_group, log_p > log_pval_cutoff)

df_combined_group<-df_combined_group%>%
  group_by(chrom, pos)%>%
  mutate(n_groups= n())

cutoffs<-filter(df_combined_group, n_groups>1)
cutoffs$pos2=cutoffs$pos

#####
###read in gene file
genes <- read.delim("dmel-all-r6.58.gggenes.txt", sep="\t", header=T)

# Extract parameters from cutoffs
# Assuming cutoffs has columns 'param1', 'param2', and 'param3' to be used as parameters
parameters <- cutoffs %>%
  select(chrom, pos, pos2)

# Define a function to look up information
lookup_info <- function(chrom, pos, pos2) {
  filtered_genes <- genes %>%
    filter(molecule == chrom,
           start <= pos2,
           end >= pos)
  
  
  result <- filtered_genes %>%
    mutate(
      is_nc = ifelse(grepl("lnc", gene, ignore.case = TRUE), gene, NA),
      is_tRNA = ifelse(grepl("tRNA", gene, ignore.case = TRUE), gene, NA),
      is_snoRNA = ifelse(grepl("snoRNA", gene, ignore.case = TRUE), gene, NA),
      is_asRNA = ifelse(grepl("asRNA", gene, ignore.case=TRUE), gene, NA),
      is_other = ifelse(!grepl("nc|tRNA|snoRNA|asRNA", gene, ignore.case = TRUE), gene, NA)
    ) %>%
    select(is_nc,is_other, is_tRNA, is_snoRNA, is_asRNA)  # Select the processed columns
  
  return(result)
}

# Apply the lookup function and combine with parameters
result_list <- lapply(1:nrow(parameters), function(i) {
  row <- parameters[i, ]
  lookup_df <- lookup_info(row[['chrom']], row[['pos']], row[['pos2']])
  
  if (nrow(lookup_df) == 0) {
    # Return a data frame with NA values if no results found
    lookup_df <- data.frame(is_nc = NA, is_asRNA=NA, is_other=NA, is_snoRNA=NA)
  }
  
  # Add parameters to the results
  cbind(row, lookup_df)
})

# Combine all results into a single data frame
result_df <- bind_rows(result_list)


#-----separate by type------
df <- result_df %>% separate(is_other, into = c("Fbgn_coding","coding" ), sep = "_")
df <- df %>% separate(is_nc, into = c("Fbgn_lncRNA", "lncRNA"), sep = "_")
df <- df %>% separate(is_asRNA, into = c( "Fbgn_asRNA","asRNA"), sep = ":")

df <- df %>% separate(lncRNA, into = c("type", "lncRNA"), sep = ":")
df <- df %>% separate(coding, into = c("coding", "type"), sep = ":")
df<-df%>% distinct()


# export
write.csv(df, "gene_ont.csv", row.names = TRUE)




#
merged <- df %>%
  inner_join(df_combined_group, by = c("pos", "chrom"), relationship = "many-to-many")

#Clean up background data (nonsignificant)
df_combined_clean <- df_combined %>%
  filter(!is.na(pos), !is.na(log_p))

df_positions <- df %>%
  select(chrom, pos) %>%
  distinct()

# Keep only background positions
nonsignificant_df <- df_combined_clean %>%
  anti_join(df_positions, by = c("chrom", "pos"))

# Sample 10% for plotting background
set.seed(123)
nonsignificant_sample <- nonsignificant_df %>%
  sample_frac(0.1)



# Combine RNA types into plot_df
lncRNA_df <- merged %>%
  filter(!is.na(lncRNA)) %>%
  transmute(pool, chrom, pos, log_p, Gene = lncRNA, Type = "lncRNA")

gene_df <- merged %>%
  filter(!is.na(coding)) %>%
  transmute(pool, chrom, pos, log_p, Gene = coding, Type = "Coding")

asRNA_df <- merged %>%
  filter(!is.na(asRNA)) %>%
  transmute(pool, chrom, pos, log_p, Gene = asRNA, Type = "asRNA")

snoRNA_df <- merged %>%
  filter(!is.na(snoRNA)) %>%
  transmute(pool, chrom, pos, log_p, Gene = snoRNA, Type = "snoRNA")

plot_df <- bind_rows(lncRNA_df, gene_df, asRNA_df, snoRNA_df)

# Add significance info
plot_df <- merged %>%
  distinct()%>%
  group_by(pool, chrom) %>%
  mutate(
    logp_threshold = quantile(log_p, 0.995, na.rm = TRUE),
    is_top_10pct = log_p >= logp_threshold
  ) %>%
  ungroup()





#### Groups separate and by chrom
ggplot() +
  # Background grey points (non-significant, sampled)
  geom_point(data = nonsignificant_sample, aes(x = pos / 1e6, y = log_p), 
             color = "grey60", alpha = 0.3, size = 1.2, show.legend = FALSE) +
  
  # Main data: colored points
  geom_point(data = merged, aes(x = pos / 1e6, y = log_p, color = "red"), size = 2, show.legend = FALSE) +

  facet_grid(pool ~ chrom, scales = "free", space = "free") +
  labs(
    x = "Chromosome Position (Mb)",
    y = "-log10(P-value)"
  ) +
  theme_bw()+
  geom_hline(data= cutoffs_combined, aes(yintercept=log_pval_cutoff), linetype="dashed")


