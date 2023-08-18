

library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id

load("intermediate_data/go_annotations_propagated_filtered_matrix.Rdata")
load(paste0(
    "intermediate_data/",
    "go_annotations_propagated_filtered_matrix_by_ontology.Rdata"))

load("product/coexp_intra_study_coexp_network_spearman_20210812.Rdata")
load("intermediate_data/blastp_network.Rdata")
load("intermediate_data/sac_physical_network.Rdata")
load("intermediate_data/sac_genetic_network.Rdata")
load("intermediate_data/coevolution_network.Rdata")

load("intermediate_data/CryptoNetV1.Rdata")
CryptoNetV1_wide <- CryptoNetV1 |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::transmute(
            cnag_id_1 = cnag_id,
            target1 = cnag_id |> stringr::str_extract("^[A-Z_0-9]+")),
        by = "target1",
        relationship = "many-to-many") |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::arrange(cnag_id) |>
        dplyr::transmute(
            cnag_id_2 = cnag_id,
            target2 = cnag_id |> stringr::str_extract("^[A-Z_0-9]+")),
        by = "target2",
        relationship = "many-to-many") |>
    dplyr::select(cnag_id_1, cnag_id_2, score) |>
    data.frame() |>
    EGAD::build_weighted_network(genes)
CryptoNetV1_wide[is.na(CryptoNetV1_wide)] <- 0

## Summary
gba <- CalCEN::network_gba(
    genes = genes,
    annotation_sets = c(
        list(all = go_annotations_propagated_filtered_matrix),
        go_annotations_propagated_filtered_matrix_by_ontology),
    networks = list(
        blastp = blastp_network,
        CryptoNetV1 = CryptoNetV1_wide,
        coexp = coexp_intra_study_coexp_network_spearman,
        SacPhysical = sac_physical_network,
        SacGenetic = sac_genetic_network),
    nfold = 10,
    verbose = TRUE)

gba |>
    readr::write_tsv(
        file = paste0(
            "product/go_annotation_prediction_summary_",
            CalCEN::date_code(), ".tsv"))


## Summary
gba_v2 <- CalCEN::network_gba(
    genes = genes,
    annotation_sets = c(
        list(all = go_annotations_propagated_filtered_matrix),
        go_annotations_propagated_filtered_matrix_by_ontology),
    networks = list(
        blastp = blastp_network,
        CryptoNetV1 = CryptoNetV1_wide,
        coexp = coexp_intra_study_coexp_network_spearman,
        coevo = coevolution_network),
    nfold = 10,
    verbose = TRUE)

# coevo
# Computing GBA for 'all' annotations, using the coevo network: mean:  0.2670345 std:  0.07612189
# Computing GBA for 'BP' annotations, using the coevo network: mean:  0.2484704 std:  0.06052518
# Computing GBA for 'CC' annotations, using the coevo network: mean:  0.296897 std:  0.08360653
# Computing GBA for 'MF' annotations, using the coevo network: mean:  0.2579861 std:  0.08418503

gba_v2 |>
    readr::write_tsv(
        file = paste0(
            "product/go_annotation_prediction_summary_v2_",
            CalCEN::date_code(), ".tsv"))


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
save(
    blastp_coexp_go_term_prediction,
    file = "product/blastp_coexp_go_term_prediction.Rdata")



##################3
coexp_full_direct_bootstrap$path |>
    purrr::map_dfr(function(network) {
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

gba |>
    readr::write_tsv(
        file = paste0(
            "product/go_annotation_prediction_summary_",
            CalCEN::date_code(), ".tsv"))



gba <- CalCEN::network_gba(
    genes = genes,
    annotation_sets = c(
        list(all = go_annotations_propagated_filtered_matrix),
        go_annotations_propagated_filtered_matrix_by_ontology),
    networks = list(
        blastp = blastp_network,
        coexp = coexp_intra_study_coexp_network_spearman,
        coexp_direct = coexp_full_direct_bootstrap$merge[[30]] |> as.matrix()),
    nfold = 3,
    degrees = 1:3,
    verbose = TRUE)

gba |>
    readr::write_tsv(
        file = paste0(
            "product/go_annotation_prediction_summary_",
            CalCEN::date_code(), ".tsv"))
