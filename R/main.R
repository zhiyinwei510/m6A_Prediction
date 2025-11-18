#' Encode DNA sequences into position-wise nucleotide factors
#'
#' This function converts DNA strings (e.g. "ATCGA") into a data frame,
#' where each column corresponds to a nucleotide position, and each value
#' is a factor with levels A, T, C, G. This encoding is useful for machine
#' learning model inputs.
#'
#' @param dna_strings A character vector containing one or more DNA sequences.
#' @return A data frame where each column represents a nucleotide position.
#' @examples
#' dna_encoding(c("ATCGA", "GGTAC"))
#' @import randomForest
#' @export
dna_encoding <- function(dna_strings){ #在每个函数上方添加help。
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A methylation probability and status for multiple samples
#'
#' This function takes a trained random forest model and a feature data frame
#' as input, and outputs predicted m6A probabilities and classification results
#' ("Positive" or "Negative") for each sample.
#'
#' @param ml_fit A trained random forest model object.
#' @param feature_df A data frame containing features such as gc_content,
#'   RNA_type, RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation, and DNA_5mer.
#' @param positive_threshold A numeric threshold (0–1) that determines whether
#'   a site is classified as "Positive" (above threshold) or "Negative". Default value is 0.5.
#' @return A data frame containing the original features and two new columns:
#'   predicted_m6A_prob and predicted_m6A_status.
#' @examples
#' # Load data objects included in the package
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_data <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#'
#' # Run prediction example
#' prediction_multiple(rf_fit, example_data, 0.6)
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold=0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  #Check errors if incorrect column names of input data.frame

  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  #feature_df$m6A_status <- factor(feature_df$m6A_status, levels = c("Negative", "Positive"))

  seq_df<-dna_encoding(feature_df$DNA_5mer)
  newdata<-cbind(feature_df,seq_df)

  prob_test<-predict(ml_fit,newdata=newdata, type="prob")

  feature_df$predicted_m6A_prob<-prob_test[,"Positive"]
  feature_df$predicted_m6A_status<-ifelse(prob_test[,"Positive"]>positive_threshold, "Positive","Negative")

  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}


#' Predict m6A methylation probability and status for a single RNA site
#'
#' This function predicts m6A probability and classification for one
#' individual observation using user-specified input feature values.
#'
#' @param ml_fit A trained random forest model object.
#' @param gc_content Numeric GC content value between 0 and 1.
#' @param RNA_type RNA type as a string ("mRNA", "lincRNA", "lncRNA", "pseudogene").
#' @param RNA_region RNA region as a string ("CDS", "intron", "3'UTR", "5'UTR").
#' @param exon_length Numeric exon length value.
#' @param distance_to_junction Numeric distance from exon junction.
#' @param evolutionary_conservation Numeric conservation score between 0 and 1.
#' @param DNA_5mer A 5-character DNA sequence string (A, T, C, G).
#' @param positive_threshold Numeric classification threshold (0–1).Default value is 0.5.
#' @return A named vector with two values:
#'   predicted_m6A_prob and predicted_m6A_status.
#' @examples
#' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(rf_fit, 0.6, "mRNA", "CDS", 12, 5, 0.8, "ATCGAT", 0.5)
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold){
  feature_df<-data.frame(gc_content=gc_content,RNA_type=RNA_type,RNA_region=RNA_region,exon_length=exon_length,distance_to_junction=distance_to_junction,evolutionary_conservation=evolutionary_conservation,DNA_5mer=DNA_5mer)
  predict_df<-prediction_multiple(ml_fit, feature_df, positive_threshold=0.5)
  returned_vector<-unlist(predict_df[,c("predicted_m6A_prob","predicted_m6A_status")])
  names(returned_vector)<-c("predicted_m6A_prob", "predicted_m6A_status")
  #注意上面两行改成vectors的格式。
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}
