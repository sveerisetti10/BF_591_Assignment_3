library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples

read_data <- function(intensity_data, delimiter) {
  reading_intensity <- read.csv(intensity_data, sep = delimiter)
  return(reading_intensity)
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples

calculate_variance_explained <- function(pca_results) {
  pca_variance <- pca_results$sdev^2
  pca_variation_percentage <- pca_variance/cumsum(pca_variance)*100

  return(pca_variation_percentage)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PCA names, variance explained by each PC, and the
#' cumulative sum of variance explained
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#'
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained and the cumulative variance explained
#' @export
#'
#' @examples 

make_variance_tibble <- function(pca_ve, pca_results) {
  pca_tibble <- tibble(pca_variance_cumulative = pca_ve) %>%
    dplyr::mutate(names = colnames(pca_results$x),
                           individual_variance = colnames(pca_results$x),
                           cumulative_variance_explained = cumsum(pca_variance_cumulative))
                  
  return(pca_tibble)
}

#' Define a function that creates a bar plot of the variance explained by each
#' PC along with a scatter plot showing the cumulative sum of variance explained
#' using ggplot2
#' @param variance_tibble (tibble): the tibble generated in the previous function
#' that contains each PC label, the variance explained by each PC, and the 
#' cumulative sum of variance explained
#'
#' @return A ggplot with a barchart representing individual variance
#'   explained and a scatterplot (connected with a line) that represents the
#'   cumulative sum of PCs
#' @export
#'
#' @examples

plot_pca_variance <- function(variance_tibble) {
 
  variance_tibble %>%
    ggplot()+
    geom_bar(aes(x=names, y=pca_variance_cumulative, fill='blue'), stat = 'identity', color = 'black') +
    geom_point(aes(x=names, y=cumulative_variance_explained, group=1))+
    geom_line(aes(x=names, y=cumulative_variance_explained, group=1))+
    labs(x='PC', y='percent variance')+

}

#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples

make_biplot <- function(metadata, pca_results) {
  
  meta_data <- read_csv(metadata) %>%
    select(SixSubtypesClassification, geo_accession)
  
  subtype_classification <- pca_results$x %>%
    as.tibble(rownames = 'geo_accession')%>% left_join(meta_data, by='geo_accession') 
  
  PC_plot <- subtype_classification %>%
    ggplot() + geom_point(aes(x = PC1, y = PC2, color = SixSubtypesClassification))
  
  return(PC_plot)
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_csv (str): The path to the differential expression results
#'   file we have provided
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the CSV.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples

list_significant_probes <- function(diff_exp_csv, fdr_threshold) {
  differential_e <- read_csv(diff_exp_csv) %>%
    as_tibble(rownames = 'name_of_probeids') %>%
    filter(padj < 0.01)
    
    differential_e$name_of_probeids
  
  return(differential_e)
  
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples

return_de_intensity <- function(intensity, sig_ids_list) {
  #I tried many times to make my data work here but I keep getting an error. 
  #However, I still intend to write code to make a heatmap in the next question
  #even though there likely will not be a resulting heatmap. 
  differentially_expressed <- intensity %>%
    as_tibble(rownames = 'name_of_probeids') %>%
    filter(intensity$name_of_probeids %in% sig_ids_list) %>%
    as.matrix()
  
  return(differentially_expressed)
}


#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples

library(RColorBrewer)
plot_heatmap <- function(de_intensity, num_colors, palette) {
   
  display.brewer.all()
  color_scheme <- colorRampPalette(brewer.pal(num_colors, palette))
  
  return(heatmap(de_intensity, col = color_scheme ))
  
  #Here even though I do not have a de_intensity parameter (since I 
  #was not able to get a result in the last problem) I tried my very best 
  #to look at our class website and try to make a generic heat map 
  #as if I actually had the de_intensity data.
}
