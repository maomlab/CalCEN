
library(tidyverse)
library(CalCEN)

load("intermediate_data/h99_transcript_annotations.Rdata")
load("intermediate_data/coexp_intra_study_coexp_spearman.Rdata")
load("intermediate_data/blastp.Rdata")


gene_gene_report <- function(data) {
    report <- data |>
        dplyr::select(cnag_id_1, cnag_id_2)

    # gene annotations
    report <- report |>
        left_join(
            h99_transcript_annotations |>
            dplyr::select(
                cnag_id_1 = cnag_id,
                gene_id_1 = gene_id,
                gene_symbol_1 = gene_symbol,
                gene_product_1 = gene_product,
                description_1 = description),
            by = "cnag_id_1") |>
        left_join(
            h99_transcript_annotations |>
            dplyr::select(
                cnag_id_2 = cnag_id,
                gene_id_2 = gene_id,
                gene_symbol_2 = gene_symbol,
                gene_product_2 = gene_product,
                description_2 = description),
            by = "cnag_id_2")

    #Co-expression
    report <- report |>
        dplyr::left_join(
            coexp_intra_study_coexp_spearman |>
            dplyr::rename(coexp_score = score),
            by = c("cnag_id_1", "cnag_id_2"))

    #BlastP
    report <- report |>
        dplyr::left_join(
            blastp |>
            dplyr::transmute(
                cnag_id_1 = ref_target |> stringr::str_replace("-[^-]+$", ""),
                cnag_id_2 = query_target |> stringr::str_replace("-[^-]+$", ""),
                blastp_EValue = EValue),
            by = c("cnag_id_1", "cnag_id_2"))
}
