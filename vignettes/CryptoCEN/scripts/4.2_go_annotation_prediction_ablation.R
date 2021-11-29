# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

library(plyr)
library(reshape2)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(EGAD)

runs <- readr::read_tsv("product/runs_20210803.tsv")

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id

load("intermediate_data/estimated_expression.Rdata")

load("intermediate_data/go_annotations_propagated_filtered_matrix.Rdata")

exprs <- reshape2::acast(
    data = estimated_expression,
    formula = cnag_id ~ run_accession,
    value.var = "FPKM") %>%
    magrittr::extract(genes, )


coexp_intra_study_coexp_network_spearman <- estimated_expression %>%
    dplyr::semi_join(runs, by = c("run_accession")) %>%
    dplyr::group_by(study_accession) %>%
    dplyr::group_map(function(study_expression, study) {
        exprs <- reshape2::acast(
            data = study_expression,
            formula = cnag_id ~ run_accession,
            value.var = "FPKM") %>%
            magrittr::extract(genes, )
        cat(
            "building network for study: ",
            study$study_accession[1],
            " with ", ncol(exprs), " samples ...\n", sep = "")
        network <- EGAD::build_coexp_network(
            exprs = exprs,
            gene.list = genes)
        ifelse(is.na(network), .5, network)
    })


#############################
# spearman rank correlation #
#############################

studies <- unique(runs$study_accession)
n_studies <- length(studies)

########
# test #
########
orders <- 3:3
networks_per_order <- 3

orders <- 1:n_studies
networks_per_order <- n_studies

ablation_scores <- orders %>%
    purrr::map_dfr(function(order) {
        cat("Building network sets of order ", order, ":\n", sep = "")
        utils::combn(x = 1:n_studies, m = order, simplify = FALSE) %>%
            sample(size = networks_per_order, replace = FALSE) %>%
            purrr::map_dfr(function(study_set) {
                run_set <- runs %>%
                    dplyr::filter(study_accession %in% studies[study_set])
                cat("  set: [", paste0(studies[study_set], collapse = ", "), "]",
                    " nruns: ", nrow(run_set), sep = "")
                network <- coexp_intra_study_coexp_network_spearman[study_set] %>%
                    purrr::reduce(`+`)
                network <- network / length(study_set)
                gba <- EGAD::run_GBA(
                    network = network,
                    labels = go_annotations_propagated_filtered_matrix)
                cat(
                    "  GBA: mean: ", gba[[1]][, 1] %>% mean(na.rm = TRUE),
                    " std: ", gba[[1]][, 1] %>% sd(na.rm = TRUE), "\n", sep = "")
                tibble::tibble(
                    order = order,
                    study_set = studies[study_set] %>% paste0(collapse = "|"),
                    gba_mean = gba[[1]][, 1] %>% mean(na.rm = TRUE),
                    gba_std = gba[[1]][, 1] %>% sd(na.rm = TRUE))
            })
    })



ablation_scores %>%
    readr::write_tsv(
        file = "product/coexp_go_prediction_summary_study_ablation_", CalCEN::date_code(), ".tsv"))

plot <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_smooth(
        data = ablation_scores,
        mapping = ggplot2::aes(
            x = order,
            y = gba_mean)) +
    ggplot2::geom_jitter(
        data = ablation_scores,
        mapping = ggplot2::aes(
            x = order,
            y = gba_mean),
        height = 0,
        width = .1,
        color = "grey20") +
    ggplot2::ggtitle(
        label = "Co-Expression GBA by study subset size") +
    ggplot2::scale_x_continuous(
        "Number of Studies",
        breaks = c(1, 5, 10, 15, 18)) +
    ggplot2::scale_y_continuous(
        "Go Prediction AUROC",
        limits = c(.625, .8))

ggplot2::ggsave(
    paste0("product/figures/coexp_go_prediction_study_ablation_", CalCEN::date_code(), ".pdf"),
    useDingbats=FALSE,
    width = 5,
    height = 4)

ggplot2::ggsave(
    "product/figures/go_pred_coexp_study_ablation_20201121.png",
    width = 5,
    height = 4)
