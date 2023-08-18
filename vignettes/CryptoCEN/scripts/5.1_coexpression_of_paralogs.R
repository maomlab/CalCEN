library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(CalCEN)

source("scripts/gene_gene_report.R")

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id

load("intermediate_data/Cneo_genes_per_og.Rdata")
load("intermediate_data/Cneo_coevo_res.Rdata")
load("intermediate_data/Cneo_coevo_gene_identifiers.Rdata")


paralogs <- dplyr::inner_join(
    Cneo_genes_per_og |>
        dplyr::rename(gene_id_1 = gene_id),
    Cneo_genes_per_og |>
        dplyr::rename(gene_id_2 = gene_id),
    by = "orthologous_group") |>
    dplyr::filter(gene_id_1 < gene_id_2) |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::select(
            cnag_id_1 = cnag_id,
            gene_id_1 = gene_id),
        by = "gene_id_1") |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::select(
            cnag_id_2 = cnag_id,
            gene_id_2 = gene_id),
        by = "gene_id_2") |>
    dplyr::select(
        -gene_id_1,
        -gene_id_2)

paralogs <- paralogs |>
    dplyr::left_join(
        paralogs |> gene_gene_report(),
        by = c("cnag_id_1", "cnag_id_2"))

paralogs <- paralogs |>
    dplyr::arrange(
        dplyr::desc(coexp_score)) |>
    dplyr::distinct(
        gene_id_1,
        gene_id_2,
        .keep_all = TRUE)


paralogs |>
    readr::write_tsv(
        "product/paralogs_report_20220607.tsv")


# https://stats.stackexchange.com/questions/264358/how-to-specify-a-distribution-for-left-skewed-data
model <- brms::brm(
    formula = brms::bf(coexp_score ~ 1),
    data = paralogs,
    family = brms::skew_normal(),
    chains = 4)
save(
    model,
    file = paste0(
        "product/figures/coevolution_analysis/",
        "model_paralogs_coexp_score_distribution_skew_normal_20220708.Rdata")

plot <- ggplot2::ggplot(
    data = paralogs) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram(
        data = paralogs,
        mapping = ggplot2::aes(
            x = coexp_score,
            y = ggplot2::after_stat(density)),
        bins = 30) +
    ggplot2::coord_cartesian(
        xlim = c(-.02, 1.02), expand = FALSE) +
    ggplot2::ggtitle(
        "Distribution of paralog co-expression scores") +
    ggplot2::scale_x_continuous("Co-expression score") +
    ggplot2::scale_y_continuous("Density")

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/",
        "paralogs_coexp_score_20230506.pdf"),
    plot = plot,
    width = 5,
    height = 4)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/",
        "paralogs_coexp_score_20230506.png"),
    plot = plot,
    width = 5,
    height = 4)
