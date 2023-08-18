

library(plyr)
library(tidyverse)

source("scripts/gene_gene_report.R")

top_coexp_hits <- coexp_intra_study_coexp_spearman |>
    dplyr::filter(cnag_id_1 != cnag_id_2) |>
    dplyr::filter(score > .95) |>
    dplyr::group_by(cnag_id_1) |>
    dplyr::arrange(desc(score)) |>
    dplyr::filter(dplyr::row_number() <= 50) |>
    dplyr::ungroup() |>
    gene_gene_report()


top_coexp_hits |>
    readr::write_tsv(
        file = paste0(
            "product/top_coexp_hits_top0.05_",
            CalCEN::date_code(), ".tsv"))

top_coexp_hits <- coexp_intra_study_coexp_spearman |>
    dplyr::filter(cnag_id_1 != cnag_id_2) |>
    dplyr::group_by(cnag_id_1) |>
    dplyr::arrange(desc(score)) |>
    dplyr::filter(dplyr::row_number() <= 50) |>
    dplyr::ungroup() |>
    gene_gene_report()


top_coexp_hits |>
    readr::write_tsv(
        file = paste0("product/top_coexp_hits_", CalCEN::date_code(), ".tsv"))



################
# CNAG_05968-t26_1 cdc420
# CNAG_05348-t26_1 cdc42

cdc420_hits <- coexp_intra_study_coexp_spearman |>
    dplyr::filter(cnag_id_1 == "CNAG_05968-t26_1") |>
    dplyr::arrange(desc(score)) |>
    gene_gene_report()
cdc420_hits |> readr::write_tsv("product/cdc420_hits_20230403.tsv")

cdc42_hits <- coexp_intra_study_coexp_spearman |>
    dplyr::filter(cnag_id_1 == "CNAG_05348-t26_1") |>
    dplyr::arrange(desc(score)) |>
    gene_gene_report()
cdc42_hits |> readr::write_tsv("product/cdc42_hits_20230403.tsv")
