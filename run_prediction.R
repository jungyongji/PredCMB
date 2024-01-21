## Load library ##
library_requirement <- c("argparse","DESeq2","data.table","BiocParallel","ggplot2","dplyr","scales","piano")
for (libr in library_requirement){
  if (libr %in% installed.packages()[,"Package"] == FALSE){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install(libr)
  }
  suppressPackageStartupMessages(require(libr, character.only = TRUE))
}

#if (length(args) <= 2) {
#  stop("At least three arguments must be supplied", call.=FALSE)
#}


## Set argument ##
parser <- ArgumentParser()

parser$add_argument("-i", "--input", action='store',
                    help="input_data",dest = 'gene_families')
parser$add_argument("-m", "--metadata", action='store',
                    help="input_metadata")
parser$add_argument("-o", "--output", action='store',
                    help="output")
parser$add_argument("-r", "--reference", type='character', action = 'store',
                    help="reference")
parser$add_argument("-c", "--core", type='integer', default=1, 
                    help="Number of threads")
parser$add_argument("-gb", "--mgx_mbx", action='store',
                    help="output", dest = 'gene_metabolite')
parser$add_argument("-p", "--pvalue", type='double', action='store',
                    help="set the p-value", dest = 'pvalue')

args <- parser$parse_args()

input_data<- args$gene_families
input_metadata <- args$metadata
output <- args$output
ref <- args$reference
cores <- args$core
gml <- args$gene_metabolite
pvalue <- args$pvalue


# DESeq2 #
## Loading DESeq2's input file (genefamily & metadata) ##
abundance_table <- data.frame(fread(input_data), row.names=1)
meta_table <- read.csv(input_metadata, sep = '\t',row.names = 1)



## Metadata processing ##
meta_table <- meta_table[1]
meta_table <- arrange(meta_table, eval(parse(text=colnames(meta_table)[1])))
abundance_table <- abundance_table[,rownames(meta_table) ]


## Excluding gene families with zero of 75% or more for all samples ##
abundance_table <- abundance_table[rowSums(abundance_table) != 0, ]
abundance_table <- abundance_table[rowMeans(abundance_table==0) < 0.75, ]


## Pick reference ##
meta_table[[1]] <- factor(meta_table[[1]])
meta_table[[1]] <- relevel(meta_table[[1]], ref)
colnames(meta_table)[1] <- "feature"



## Create DESeq2 data set ##
dds <- DESeqDataSetFromMatrix(countData = abundance_table,
                                 colData = meta_table,
                                 design = ~ 1 + feature)


## Run DESeq2 ##
# Use SnowParam on Windows
if (.Platform$OS.type == "windows") {
	  BPPARAM <- SnowParam(cores)
} else {
	  BPPARAM <- MulticoreParam(cores)
}

dds_res <- DESeq(dds, sfType = "poscounts", parallel = T, BPPARAM = BPPARAM)
res <- results(dds_res, alpha = pvalue, name = resultsNames(dds_res)[2], parallel = T, BPPARAM = BPPARAM )
res <- as.data.frame(res)
#write.table(res, output, sep = '\t')


# Reporter analysis #
## Loading Reporter analysis's input file (p-value & fold_change) ##
signres <- subset(res, pvalue > 0 & pvalue < 0.05 & abs(log2FoldChange) > 1)
signres_pv <- signres['pvalue']
signres_fc <- signres['log2FoldChange']


## Loading the gene-metabolite information file ##
gene_metabolite <- fread(gml)
gene_metabolite <- loadGSC(gene_metabolite)

#saveRDS(met2, file ='./interactions/uniref_com_e_oc.rds')
#met3 <- readRDS('./interactions/uniref_com_kerk_oc.rds')



## Run Reporter analysis ##
gsa_rep <- runGSA(signres_pv,
                  signres_fc,
                  gsc=gene_metabolite,
                  geneSetStat = "reporter",
                  signifMethod = "nullDist", nPerm=1000,
                  gsSizeLim = c(1, Inf),
                  adjMethod = "BH")



## Save predicted metabolite's statistics ##
gsa_table <- GSAsummaryTable(gsa_rep)
metabolite_name <- read.csv('./interactions/raw/kegg_cpd_name_v2021-11-24.tsv', sep = '\t')
merged_data <- merge(gsa_table, metabolite_name, by.x = "Name", by.y = "cpd", all.x = TRUE)
selected_data <- merged_data[c(1, ncol(merged_data), 2:5) ]
colnames(selected_data) <- c("metabolite", "name", "genes", "statistics", "p-value", "adj-pvalue")
GSAsummaryTable(selected_data, save = TRUE, output)


# Example #
#Rscript run_prediction.R -i abun_int_prism_cpm.tsv -m prism_cohort_CD.tsv -o output.tsv -r "Control" -c 15 -gb ../reporter/uniref_com_kerk_oc.tsv -p 0.05
