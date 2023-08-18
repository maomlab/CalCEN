library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(CalCEN)

source("scripts/gene_gene_report.R")

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id

load("intermediate_data/coexp_intra_study_coexp_spearman.Rdata")
load("intermediate_data/blastp.Rdata")
load("intermediate_data/sac_physical_network_long.Rdata")
load("intermediate_data/sac_genetic_network_long.Rdata")
load("intermediate_data/coevolution_network_long.Rdata")


source("scripts/gene_gene_report.R")


############

top_coevolution_hits <- coevolution_network_long |>
    dplyr::filter(
        coevolution_coefficient == 1) |>
    gene_gene_report() |>
    dplyr::select(-gene_product_1, -gene_product_2) |>
    readr::write_tsv(
        "product/figures/coevolution_analysis/top_coevolution_hits.tsv")

sample_top_coevolution_hits <- coevolution_network_long |>
    dplyr::filter(
        coevolution_coefficient == 1) |>
    gene_gene_report() |>
    dplyr::filter(
       !(gene_product_1 |> stringr::str_detect("hypothetical protein")),
       !(gene_product_2 |> stringr::str_detect("hypothetical protein")))

sample_top_coevolution_hits |>
    dplyr::mutate(
        description_1 = description_1 |> stringr::str_sub(1, 40),
        description_2 = description_2 |> stringr::str_sub(1, 40),
        gene_product_1 = gene_product_1 |> stringr::str_sub(1, 40),
        gene_product_2 = gene_product_2 |> stringr::str_sub(1, 40)) |>
    data.frame()

################


coexp_coevo <- dplyr::inner_join(
    coexp_intra_study_coexp_spearman |>
        dplyr::rename(coexp_score = score),
    coevolution_network_long |>
    dplyr::select(
        cnag_id_1, cnag_id_2, coevo_score = coevolution_coefficient),
    by = c("cnag_id_1", "cnag_id_2"))

 plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Co-Expression vs. Co-Evolution Scores") +
    ggplot2::scale_x_continuous(
        "Co-Expression Spearman Rank Correlation (Bigger is Better)",
        limits = c(-0.02, 1.02),
        expand = FALSE) +
    ggplot2::scale_y_continuous(
        "Co-Evolution Coefficient (Bigger is Better)",
        limits = c(-1.02, 1.02),
        expand = FALSE) +
    ggrastr::geom_point_rast(
        data = coexp_coevo,
        mapping = ggplot2::aes(
            x = coexp_score,
            y = coevo_score),
        alpha = 0.1,
        size = .2,
        shape = 16) +
    ggplot2::geom_density_2d(
        data = coexp_coevo |>
            dplyr::sample_n(size = 1000000),
        mapping = ggplot2::aes(
            x = coexp_score,
            y = coevo_score),
        linewidth = 0.8,
        color = "gray50") +
    ggplot2::geom_smooth(
        data = coexp_coevo |>
            dplyr::sample_n(size = 100000),
        mapping = ggplot2::aes(
            x = coexp_score,
            y = coevo_score),
        method = "lm",
        formula = y ~ x)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/coexp-vs-coevo_",
        CalCEN::date_code(), ".pdf"),
    width = 5,
    height = 5,
    useDingbats = FALSE)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/coexp-vs-coevo_",
        CalCEN::date_code(), ".png"),
    width = 5,
    height = 5)



 plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Co-Expression of Significant Co-Evolved Genes") +
    ggplot2::scale_x_continuous(
        "Co-Expression Spearman Rank Correlation (Bigger is Better)") +
    ggplot2::scale_y_continuous("Density") +
    ggplot2::geom_histogram(
        data = coexp_coevo |>
            dplyr::filter(coevo_score == 1),
        mapping = ggplot2::aes(
            x = coexp_score,
            y = ggplot2::after_stat(density)),
        bins = 30) +
    ggplot2::coord_cartesian(
        xlim = c(-.02, 1.02), expand = FALSE)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/coexp_signif_coevo_",
        CalCEN::date_code(), ".pdf"),
    width = 5,
    height = 4)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/coexp_signif_coevo_",
        CalCEN::date_code(), ".png"),
    width = 5,
    height = 4)

#########################################
blastp_coevo <- dplyr::inner_join(
    blastp |>
    dplyr::transmute(
        cnag_id_1 = ref_target |> stringr::str_replace("-[^-]+$", ""),
        cnag_id_2 = query_target |> stringr::str_replace("-[^-]+$", ""),
        blastp_score = bit_score),
    coevolution_network_long |>
    dplyr::select(
        cnag_id_1, cnag_id_2, coevo_score = coevolution_coefficient),
    by = c("cnag_id_1", "cnag_id_2"))

 plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("BlastP vs. Co-Evolutionion Scores") +
    ggplot2::scale_x_log10(
        "BlastP Bit Score (Bigger is Better)") +
    ggplot2::scale_y_continuous(
        "Co-Evolution Coefficient (Bigger is Better)",
        limits = c(-1.02, 1.02),
        expand = c(0, 0)) +
    ggrastr::geom_point_rast(
        data = blastp_coevo,
        mapping = ggplot2::aes(
            x = blastp_score,
            y = coevo_score),
        alpha = 0.1,
        size = .4,
        shape = 16) +
    ggplot2::geom_smooth(
        data = blastp_coevo,
        mapping = ggplot2::aes(
            x = blastp_score,
            y = coevo_score),
        method = "lm",
        formula = y ~ x)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/blastp-vs-coevo_",
        CalCEN::date_code(), ".pdf"),
    width = 5,
    height = 5,
    useDingbats = FALSE)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/blastp-vs-coevo_",
        CalCEN::date_code(), ".png"),
    width = 5,
    height = 5)




#####################################################################33
sac_physical_coevo <- dplyr::inner_join(
    sac_physical_network_long |>
    dplyr::transmute(
        cnag_id_1 = target1,
        cnag_id_2 = target2,
        sac_phys_score = max_EValue),
    coevolution_network_long |>
    dplyr::select(
        cnag_id_1, cnag_id_2, coevo_score = coevolution_coefficient),
    by = c("cnag_id_1", "cnag_id_2"))

 plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Sac Physical vs. Co-Evolutionion Scores") +
    ggplot2::scale_x_log10(
        "Sac Physical EValue (Smaller is Better)") +
    ggplot2::scale_y_continuous(
        "Co-Evolution Coefficient (Bigger is Better)",
        limits = c(-1.02, 1.02),
        expand = c(0, 0)) +
    ggrastr::geom_point_rast(
        data = sac_physical_coevo,
        mapping = ggplot2::aes(
            x = sac_phys_score,
            y = coevo_score),
        alpha = 0.1,
        size = .4,
        shape = 16) +
    ggplot2::geom_smooth(
        data = sac_physical_coevo,
        mapping = ggplot2::aes(
            x = sac_phys_score,
            y = coevo_score),
        method = "lm",
        formula = y ~ x)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/sac_phys-vs-coevo_",
        CalCEN::date_code(), ".pdf"),
    width = 5,
    height = 5,
    useDingbats = FALSE)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coevolution_analysis/sac_phys-vs-coevo_",
        CalCEN::date_code(), ".png"),
    width = 5,
    height = 5)
