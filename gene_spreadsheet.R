cutoffs<- df_combined_i %>%
  group_by(pop_1_2)%>%
  reframe(p_val=quantile(p_val,c(0.75, 0.95, 0.99, 0.995, 0.999)), q=c(0.75, 0.95, 0.99, 0.995, 0.999)) %>%
  filter(q==0.99)%>% 
  rename(log_pval_cutoff=p_val)

df_combined_i<- left_join(df_combined_i,cutoffs)

cutoffs<- filter(df_combined_i, p_val > log_pval_cutoff)

cutoffs<-cutoffs%>%
  group_by(chrom, window)%>%
  mutate(n_groups= n())


cutoffs<-filter(cutoffs, n_groups>1)

#####

###read in gene file
genes <- read.delim("dmel-all-r6.58.gggenes.txt", sep="\t", header=T)


# Load necessary library
library(dplyr)

# Step 1: Read the data frames
cutoffs <- read.csv("data_frame1.csv")
df2 <- read.csv("data_frame2.csv")

# Step 2: Extract parameters from cutoffs
# Assuming cutoffs has columns 'param1', 'param2', and 'param3' to be used as parameters
parameters <- cutoffs %>%
  select(chrom, window_start, window_end)

# Step 3: Define a function to look up information in df2
lookup_info <- function(chrom, window_start, window_end) {
  # Filter df2 based on parameters
  filtered_genes <- genes %>%
    filter(molecule == chrom, # replace with actual conditions
           start > window_start,
           end < window_end)
  
  result <- filtered_genes %>%
    mutate(
      # Replace 'specific_string' with the actual string to check
      result_column1 = ifelse(grepl("nc", gene), gene, NA),
      result_column2 = ifelse(!grepl("nc", gene), gene, NA)
    ) %>%
    select(result_column1, result_column2)  # Select the processed columns
  
  return(result)
}

# Step 4: Apply the lookup function and combine with parameters
result_list <- lapply(1:nrow(parameters), function(i) {
  row <- parameters[i, ]
  lookup_df <- lookup_info(row[['chrom']], row[['window_start']], row[['window_end']])
  
  if (nrow(lookup_df) == 0) {
    # Return a data frame with NA values if no results found
    lookup_df <- data.frame(result_column1 = NA, result_column2 = NA)
  }
  
  # Add parameters to the results
  cbind(row, lookup_df)
})

# Combine all results into a single data frame
result_df <- bind_rows(result_list)

df <- result_df %>% separate(result_column2, into = c("Gene_Name", "Annotation_Name"), sep = "_")
df <- df %>% separate(result_column1, into = c("Non-coding", "Annotation"), sep = "_")

# Step 5: Export the resulting data frame
write.csv(df, "gene_ont.csv", row.names = FALSE)
