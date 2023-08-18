

library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")
load("intermediate_data/coexp_intra_study_coexp_spearman.Rdata")
madhani_mutant_collection <- readxl::read_excel(
    path = "raw_data/madhani mutant collections.xlsx",
    skip = 1) |>
    dplyr::filter(
        !is.na(`Quality control colony PCR`),
        `Quality control colony PCR` == "PASS",
        CNAG |> stringr::str_detect("^CNAG")) |>
    dplyr::select(gene_id = CNAG) |>
    dplyr::left_join(
        h99_transcript_annotations |>
            dplyr::arrange(gene_id, cnag_id) |> # take the first variant only
            dplyr::distinct(gene_id),
        by = "gene_id")

save(
    madhani_mutant_collection,
    file = "intermediate_data/madhani_mutant_collection.Rdata")
   
load("intermediate_data/madhani_mutant_collection.Rdata")

dna_damage_hits <- readr::read_tsv(
    file = "intermediate_data/dna_damage_hits.tsv",
    show_col_types = FALSE) |>
    dplyr::left_join(node_weight, by = "gene_id") |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::select(gene_id, description, gene_symbol),
        by = "gene_id")


node_weight <- coexp_intra_study_coexp_spearman |>
    dplyr::group_by(cnag_id_1) |>
    dplyr::summarize(
        weight = sum(score),
        .groups = "drop") |>
    dplyr::rename(cnag_id = cnag_id_1) |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::select(
            cnag_id,
            gene_id),
        by = "cnag_id") |>
    dplyr::arrange(weight) |>
    dplyr::mutate(
        weight_rank = rank(weight) / dplyr::n())

madhani_mutant_collection <- madhani_mutant_collection |>
    dplyr::left_join(
        node_weight,
        by = "gene_id")

dna_damage_matched_controls <- dna_damage_hits |>
    dplyr::rowwise() |>
    dplyr::do({
        data <- .
        matched_control <- madhani_mutant_collection |>
            dplyr::filter(gene_id != data$gene_id[1]) |>
            dplyr::mutate(weight_diff = abs(weight - data$weight[1])) |>
            dplyr::arrange(weight_diff) |>
            head(10)
        data.frame(
            hit_gene_id = data$gene_id[1],
            hit_weight = data$weight[1],
            hit_weight_rank = data$weight_rank[1],
            matched_control_gene_id = matched_control$gene_id,
            matched_control_weight = matched_control$weight,
            matched_control_weight_rank = matched_control$weight_rank,
            matched_control_rank = seq(1, nrow(matched_control)),
            weight_diff = matched_control$weight_diff)
    }) |>
    dplyr::ungroup()

dna_damage_matched_controls |>
    readr::write_tsv(
        file = "product/dna_damage_matched_controls_20230404.tsv")
