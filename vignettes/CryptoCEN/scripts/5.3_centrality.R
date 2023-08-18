

library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")
load("intermediate_data/coexp_intra_study_coexp_spearman.Rdata")

load("product/coexp_intra_study_coexp_network_spearman_20210812.Rdata")

load("intermediate_data/go_annotations_propagated_filtered_matrix.Rdata")

go_terms <- readr::read_tsv(
    file = "raw_data/go_term_20230506.tsv",
    show_col_types = FALSE)


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



gba <- EGAD::run_GBA(
    network = coexp_intra_study_coexp_network_spearman,
    labels = go_annotations_propagated_filtered_matrix,
    min = 20,
    max = 1000,
    nfold = 10)

scores <- gba[[1]] |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "go_id") |>
    dplyr::left_join(
        go_terms,
        by = "go_id")



p <- ggplot2::ggplot(
    data = scores) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(.8, .3)) +
    ggplot2::geom_abline(
        slope = 1,
        color = "grey70") +
    ggplot2::geom_point(
        mapping = ggplot2::aes(
            x = degree_null_auc,
            y = auc,
            color = ontology),
        shape = 16,
        alpha = 0.9,
        size = 3) +
    ggrepel::geom_label_repel(
        data = scores |>
            dplyr::mutate(
                auc_rank = dplyr::n() - rank(auc) + 1,
                degree_null_auc_rank = dplyr::n() - rank(degree_null_auc) + 1) |>
            dplyr::filter(
                term == "DNA repair" |
                auc_rank <= 5 |
                degree_null_auc_rank <= 5),
        mapping = ggplot2::aes(
            x = degree_null_auc,
            y = auc,
            label = term),
        min.segment.length = 0,
        box.padding = 0.5,
        size = 2) +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(
        "Degree-Null AUROC",
        limits = c(0, 1)) +
    ggplot2::scale_y_continuous(
        "Co-Expression Network AUROC",
        limits = c(0, 1)) +
    ggplot2::ggtitle(
        label = "GO term prediction : Co-Expession")

ggplot2::ggsave(
    file = paste0(
        "product/figures/go_pred_coexp_scatter_",
        CalCEN::date_code(), ".pdf"),
  height = 7,
  width = 7,
  useDingbats = FALSE)

ggplot2::ggsave(
    file = paste0(
        "product/figures/go_pred_coexp_scatter_",
        CalCEN::date_code(), ".png"),
  height = 7,
  width = 7)

