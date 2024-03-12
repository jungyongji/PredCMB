# Load libraries
library(readxl)
library(dplyr)
library(argparse)


# Set up argument parsing
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help="Input data file", dest="input_data")
parser$add_argument("-s", "--size", type="integer", default=5, help="Minimum number of metabolites")
parser$add_argument("-o", "--output", help="Output file")

args <- parser$parse_args()


# Read input data
met_data <- read.csv(args$input_data, sep='\t')
met_data <- met_data[c('metabolite_KEGG_ID', 'metabolite_z_value')]
colnames(met_data) <- c('cpd', 'stat')

# Standardize statistics
met_data$stat <- (met_data$stat - mean(met_data$stat)) / sd(met_data$stat)


# Read and preprocess HMDB class data
hmdb_class <- read.csv('interactions/raw/hmdb_class_curation.tsv', sep='\t', stringsAsFactors=FALSE)
hmdb_class <- unique(hmdb_class[c('cpd', 'class')])

# Merge and filter by class occurrence
merged_data <- merge(met_data, hmdb_class, by='cpd')
class_count <- merged_data %>% count(class) %>% filter(n >= args$size)

# Apply criteria for class inclusion
final_data <- merged_data %>% filter(class %in% class_count$class)

# Wilcoxon signed-rank test
results <- list()
unique_classes <- unique(final_data$class)

for (class in unique_classes) {
  in_class <- final_data$stat[final_data$class == class]
  
  if (length(in_class) < args$size) next
  
  out_class <- final_data$stat[final_data$class != class]
  
  if (length(out_class) < args$size) next
  
  test_res <- wilcox.test(out_class, in_class, alternative="two.sided", exact=FALSE)
  
  if (!is.na(test_res$p.value)) {
    results <- rbind(results, data.frame(
      class_name=class,
      class_z_median=median(in_class),
      class_z_mean=mean(in_class),
      class_p_value=test_res$p.value
    ))
  }
}

# Adjust for multiple hypothesis testing
results$class_adjusted_p_value <- p.adjust(results$class_p_value, method='BH')

# Output
write.table(results, args$output, row.names=FALSE, sep='\t')


# Rscript run_summary.R -i output_ex.tsv -s 5 -o summary_results.tsv