
library(tidyverse)
load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id


load("intermediate_data/go_annotations_propagated_filtered_matrix.Rdata")
load(paste0(
    "intermediate_data/",
    "go_annotations_propagated_filtered_matrix_by_ontology.Rdata"))


load("product/coexp_network_full_spearman_20230612.Rdata")
load("product/coexp_intra_study_coexp_network_spearman_20210812.Rdata")
load("intermediate_data/coexp_full_direct_bootstrap_60.Rdata")

if (!dir.exists("product/figures/direct_network")) {
    cat("creating 'product/figures/direct_network' ...\n")
    dir.create("product/figures/direct_network")
}


gba <- CalCEN::network_gba(
    genes = genes,
    annotation_sets = c(
        list(all = go_annotations_propagated_filtered_matrix)),
    networks = coexp_intra_study_coexp_network_spearman,
    nfold = 3,
    degrees = 1,
    verbose = TRUE)

n_total_edges <- sum(coexp_intra_study_coexp_network_spearman > 0) / 2

spearman_sparsity <- data.frame(network_index = 1:60) |>
    dplyr::rowwise() |>
    dplyr::do({
        network_index <- .$network_index[1]
        network <- coexp_intra_study_coexp_network_spearman
        network <- ifelse(network >= network_index / 60, network, 0)
        n_edges <- sum(network > 0) / 2

        cat(
            "network_index = ", network_index,
            " n_edges: ", n_edges, " fraction = ",
            n_edges / n_total_edges, "\n", sep = "")
        gba <- CalCEN::network_gba(
            genes = genes,
            annotation_sets = c(
                list(all = go_annotations_propagated_filtered_matrix)),
            networks = network,
            nfold = 3,
            degrees = 1,
            verbose = TRUE) |>
            dplyr::mutate(n_edges = n_edges)

    }) |>
    dplyr::mutate(
        bootstrap_index = dplyr::row_number(),
        .before = 3)
spearman_sparsity |> readr::write_tsv(
    file = paste0(
        "product/figures/direct_network/sparman_sparsity_gba_",
        CalCEN::date_code(), ".tsv"))


# sparsity of the direct network
direct_sparsity <- coexp_full_direct_bootstrap_60$merge |>
    purrr::map_dfr(function(network) {
        network <- as.matrix(network)
        gba <- CalCEN::network_gba(
            genes = genes,
            annotation_sets = c(
                list(all = go_annotations_propagated_filtered_matrix)),
            networks = network,
            nfold = 3,
            degrees = 1,
            verbose = TRUE) |>
            dplyr::mutate(
                n_edges = sum(network  > 0) / 2)
    }) |>
    dplyr::mutate(
        bootstrap_index = dplyr::row_number(),
        .before = 3)
direct_sparsity |> readr::write_tsv(
    file = paste0(
        "product/figures/direct_network/direct_sparsity_gba_",
        CalCEN::date_code(), ".tsv"))




###########33

plot <- ggplot2::ggplot(
    data = dplyr::bind_rows(
        spearman_sparsity |> dplyr::mutate(method = "Spearman"),
        direct_sparsity |> dplyr::mutate(method = "Direct"))) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = log10(n_edges),
            y = degree_null_auroc_mean,
            color = method),
        linewidth = 0.8) +
    ggplot2::geom_ribbon(
        mapping = ggplot2::aes(
            x = log10(n_edges),
            y = auroc_mean,
            ymin = auroc_mean - auroc_std,
            ymax = auroc_mean + auroc_std,
            fill = method),
        alpha = 0.2) +
    ggplot2::geom_line(
        mapping = ggplot2::aes(
            x = log10(n_edges),
            y = auroc_mean,
            color = method),
        linewidth = 1.2) +
    ggplot2::ggtitle(
        "Direct Network GO term prediction AUROC") +
    ggplot2::scale_x_continuous(
        "log10(N edges)",
        breaks = log10(c(1000, 10000, 100000, 1000000, 10000000, 100000000)),
        labels = c("1k", "10k", "100k", "1M", "10M", "100M")) +
    ggplot2::scale_y_continuous("AUROC", limits = c(0, 1), expand = c(0, 0))

ggplot2::ggsave(
    filename = paste0(
        "product/figures/direct_network/AUROC_by_bootstrap_index_",
        CalCEN::date_code(), ".pdf"),
    width = 5,
    height = 5)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/direct_network/AUROC_by_bootstrap_index_",
        CalCEN::date_code(), ".png"),
    width = 5,
    height = 5)
