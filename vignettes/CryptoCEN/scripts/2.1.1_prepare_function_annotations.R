library(plyr)
library(tidyverse)
library(GO.db)
library(CalCEN)

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")
load("intermediate_data/go_annotations.Rdata")

load("intermediate_data/go_annotations_sac.Rdata")

go_term_metadata <- GO.db::GO_dbconn() |>
    dplyr::tbl("metadata") |>
    dplyr::collect(n = Inf)
go_term_metadata |>
    readr::write_tsv(
        file = paste0(
            "raw_data/go_term_metadata_", CalCEN::date_code(), ".tsv"))

go_terms <- GO.db::GO_dbconn() |>
    dplyr::tbl("go_term") |>
    dplyr::collect(n = Inf)

go_terms |>
    readr::write_tsv(
        file = paste0("raw_data/go_term_", CalCEN::date_code(), ".tsv"))


# note that the same gene and annotation can come from
# difference envidence sources
go_annotations_propagated <- CalCEN::propagate_go_annotations(
    go_annotations = go_annotations,
    verbose = TRUE)

save(
    go_annotations_propagated,
    file = "intermediate_data/go_annotations_propagated.Rdata")

go_annotations_propagated_filtered <- go_annotations_propagated |>
    dplyr::filter(is.na(qualifier) | qualifier != "NOT") |>
    dplyr::distinct(db_object_id, go_id) |>
    dplyr::semi_join(
        go_annotations_propagated |>
        dplyr::filter(is.na(qualifier) | qualifier != "NOT") |>
        dplyr::distinct(db_object_id, go_id) |>
        dplyr::count(go_id) |>
        dplyr::filter(20 <= n, n <= 1000),
        by = c("go_id"))

go_annotations_propagated_filtered |> nrow()
# 14215 total annotations

save(
    go_annotations_propagated_filtered,
    file = "intermediate_data/go_annotations_propagated_filtered.Rdata")

go_annotations_propagated_filtered_matrix <- go_annotations_propagated_filtered |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::select(cnag_id, gene_id),
        by = c(db_object_id = "gene_id")) |>
    dplyr::transmute(cnag_id, go_id, value = 1) |>
    tidyr::spread(go_id, value, fill = 0) |>
    as.data.frame() |>
    tibble::column_to_rownames("cnag_id")

save(
    go_annotations_propagated_filtered_matrix,
    file = "intermediate_data/go_annotations_propagated_filtered_matrix.Rdata")



go_annotations_propagated_filtered_matrix_by_ontology <-
    go_terms |>
    dplyr::filter(ontology != "universal") |>
    split(f = ~ as.factor(ontology)) |>
    purrr::map(function(ontology_terms) {
        go_annotations_propagated_filtered_matrix[
          , colnames(go_annotations_propagated_filtered_matrix) %in%
            ontology_terms$go_id]
    })

save(
    go_annotations_propagated_filtered_matrix_by_ontology,
    file = "intermediate_data/go_annotations_propagated_filtered_matrix_by_ontology.Rdata")



###################
# SAC Annotations #
###################
go_annotations_sac_propagated <- CalCEN::propagate_go_annotations(
    go_annotations = go_annotations_sac |>
        dplyr::filter(`Feature type` == "ORF|Verified") |>
        dplyr::select(ORF, SGDID, go_id),
    verbose = TRUE)

save(
    go_annotations_sac_propagated,
    file = "intermediate_data/go_annotations_sac_propagated.Rdata")

go_annotations_sac_propagated_filtered <- go_annotations_sac_propagated |>
    dplyr::distinct(SGDID, go_id) |>
    dplyr::semi_join(
        go_annotations_sac_propagated |>
        dplyr::distinct(SGDID, go_id) |>
        dplyr::count(go_id) |>
        dplyr::filter(20 <= n, n <= 1000),
        by = c("go_id"))

# 42129 total annotations
save(
    go_annotations_sac_propagated_filtered,
    file = "intermediate_data/go_annotations_sac_propagated_filtered.Rdata")


go_annotations_sac_propagated_filtered_matrix <- go_annotations_sac_propagated_filtered |>
    dplyr::transmute(sgd_id = SGDID, go_id, value = 1) |>
    tidyr::spread(go_id, value, fill = 0) |>
    as.data.frame() |>
    tibble::column_to_rownames("sgd_id")

save(
    go_annotations_sac_propagated_filtered_matrix,
    file = "intermediate_data/go_annotations_sac_propagated_filtered_matrix.Rdata")


load("intermediate_data/go_terms_sac.Rdata")
go_annotations_sac_propagated_filtered_matrix_by_ontology <-
    go_terms_sac |>
    split(f = ~ as.factor(ontology)) |>
    purrr::map(function(ontology_terms) {
        go_annotations_sac_propagated_filtered_matrix[
          , colnames(go_annotations_sac_propagated_filtered_matrix) %in%
            ontology_terms$go_id]
    })

save(
    go_annotations_sac_propagated_filtered_matrix_by_ontology,
    file = "intermediate_data/go_annotations_sac_propagated_filtered_matrix_by_ontology.Rdata")


