library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  counts_raw <- read_tsv(counts_csv) %>% 
    column_to_rownames("gene") %>%
    as.matrix()
  
  sample_metadata <- read_csv(metafile_csv) %>% 
    dplyr::select(-1) %>% 
    dplyr::filter(timepoint %in% selected_times) %>% 
    dplyr::filter(samplename %in% colnames(counts_raw))
  
  counts_raw <- counts_raw[, sample_metadata$samplename]
  
  coldata <- sample_metadata %>%
    mutate(timepoint = factor(timepoint)) %>%
    as("DataFrame")
  
  rownames(coldata) <- coldata$samplename
  
  counts_raw <- counts_raw[, rownames(coldata)]
  
  coldata$timepoint <- relevel(coldata$timepoint, ref = "vP0")
  
  se <- SummarizedExperiment(
    assays = list(counts = counts_raw),
    colData = coldata
  )
  
    return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  res <- results(dds)
  res_df <- as.data.frame(res)
  
  return(list(
    dds = dds,
    results = res_df
  ))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  deseq_tib <- deseq2_res %>% 
    as.data.frame() %>% 
    rownames_to_column("genes") %>%
    as_tibble() %>% 
    mutate(volc_plot_status = case_when(
      padj < padj_threshold & log2FoldChange > 0  ~ "UP",
      padj < padj_threshold & log2FoldChange < 0~ "DOWN",
      TRUE ~ "NS",
    ))
    
    return(deseq_tib)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  p1 <- ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(binwidth = .1, fill = "coral", color = "black") +
    labs(title = "Raw P-Value Histogram of DESeq2 Results",
         x = "P-Value",
         y = "Count") +
    theme_minimal()
    return(p1)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig_genes <- labeled_results %>% 
    filter(padj < padj_threshold)
  
  p2 <- ggplot(sig_genes, aes(x = log2FoldChange)) +
    geom_histogram(binwidth = .1, fill = "coral", color = "black") +
    labs(title = "Raw P-Value Histogram of DESeq2 Results",
         x = "Log2FoldChange",
         y = "Count") +
    theme_minimal()
    return(p2)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  
  top_genes <- labeled_results %>% 
    filter(!is.na(padj)) %>%
    arrange(padj) %>% 
    head(num_genes) %>%
    pull(genes)
  
  norm_counts <- counts(dds_obj, normalized = TRUE)
  
  norm_counts_top <- norm_counts[top_genes, ]
  
  plot_df <- as.data.frame(norm_counts_top) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      -gene,
      names_to = "samplename",
      values_to = "normalized_count"
    )
  
  # --- 6. Add metadata ---
  plot_df <- plot_df %>%
    dplyr::left_join(
      as.data.frame(SummarizedExperiment::colData(dds_obj)) %>%
        dplyr::select(samplename, timepoint),
      by = "samplename"
    )
  
  # --- 7. Plot ---
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = gene, y = normalized_count)
  ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = timepoint),
      width = 0.2,
      size = 2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = samplename),
      vjust = -0.5,
      size = 3
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Gene",
      y = "Normalized Counts",
      title = paste("Top", num_genes, "Differentially Expressed Genes")
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
    return(p)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  p <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point(size = .5)
    return(p)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene <- read_tsv(id2gene_path, col_names = c("genes", "symbol")) 
  
  ranked_results <- labeled_results %>% 
    left_join(id2gene, by = "genes") %>% 
    filter(!is.na(symbol) & symbol != "" & !is.na(log2FoldChange)) %>% 
    group_by(symbol) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>% 
    ungroup() %>%
    arrange(desc(log2FoldChange))
  
  log2fc_vec <- ranked_results$log2FoldChange
  names(log2fc_vec) <- ranked_results$symbol
  
    return(log2fc_vec)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- gmtPathways(gmt_file_path)
  
  fgsea_results <- fgsea(
    pathways  = pathways,
    stats     = rnk_list,
    minSize   = min_size,
    maxSize   = max_size
  )
  
  fgsea_results <- fgsea_results %>%
    as_tibble() %>%
    arrange(desc(padj))
  
    return(fgsea_results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  top <- bind_rows(
    fgsea_results %>% slice_max(order_by = NES, n = num_paths),
    fgsea_results %>% slice_min(order_by = NES, n = num_paths)
  ) %>%
    mutate(
      direction = ifelse(NES > 0, "Positive", "Negative"),
      pathway = str_replace_all(pathway, "_", " ")
    ) %>%
    arrange(NES)
  
  p <- ggplot(top, aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c("Positive" = "Blue", "Negative" = "Red")) +
    theme(axis.text.y = element_text(size = 5),) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) 
    theme_minimal()
  
    return(p)
}

