

library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id

load("intermediate_data/go_annotations_propagated_filtered_matrix.Rdata")
load("intermediate_data/go_annotations_propagated_filtered_matrix_by_ontology.Rdata")

load("product/coexp_intra_study_coexp_network_spearman_20210812.Rdata")
load("intermediate_data/blastp_network.Rdata")


## Summary
gba <- CalCEN::network_gba(
    genes = genes,
    annotation_sets = c(
        list(all = go_annotations_propagated_filtered_matrix),
        go_annotations_propagated_filtered_matrix_by_ontology),
    networks = list(
        blastp = blastp_network,
        coexp = coexp_intra_study_coexp_network_spearman),
    nfold = 3,
    degrees = 1:2,
    verbose = TRUE)

gba %>%
    readr::write_tsv(
        file = paste0("product/go_annotation_prediction_summary_", CalCEN::date_code(), ".tsv"))
        

coexp_go_term_prediction <- EGAD::run_GBA(
    network = coexp_intra_study_coexp_network_spearman,
    labels = go_annotations_propagated_filtered_matrix)
save(coexp_go_term_prediction, file = "product/coexp_go_term_prediction.Rdata")

blastp_go_term_prediction <- EGAD::run_GBA(
    network = blastp_network,
    labels = go_annotations_propagated_filtered_matrix)
save(blastp_go_term_prediction, file = "product/blastp_go_term_prediction.Rdata")

blastp_coexp_go_term_prediction <- EGAD::run_GBA(
    network = blastp_network + coexp_intra_study_coexp_network_spearman,
    labels = go_annotations_propagated_filtered_matrix)
save(blastp_coexp_go_term_prediction, file = "product/blastp_coexp_go_term_prediction.Rdata")



##################3
coexp_full_direct_bootstrap$path %>%
    purrr::map_dfr(function(network){
    network <- as.matrix(network)
    gba <- CalCEN::network_gba(
        genes = genes,
        annotation_sets = c(
            list(all = go_annotations_propagated_filtered_matrix)),
        networks = network,
        nfold = 3,
        degrees = 1,
        verbose = TRUE)
    })

gba %>%
    readr::write_tsv(
        file = paste0("product/go_annotation_prediction_summary_", CalCEN::date_code(), ".tsv"))



gba <- CalCEN::network_gba(
    genes = genes,
    annotation_sets = c(
        list(all = go_annotations_propagated_filtered_matrix),
        go_annotations_propagated_filtered_matrix_by_ontology),
    networks = list(
        blastp = blastp_network,
        coexp = coexp_intra_study_coexp_network_spearman,
        coexp_direct = coexp_full_direct_bootstrap$merge[[30]] %>% as.matrix()),
    nfold = 3,
    degrees = 1:3,
    verbose = TRUE)

gba %>%
    readr::write_tsv(
        file = paste0("product/go_annotation_prediction_summary_", CalCEN::date_code(), ".tsv"))
