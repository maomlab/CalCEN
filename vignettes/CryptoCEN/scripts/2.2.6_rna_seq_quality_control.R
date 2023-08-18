library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(ggplot2)
library(scales)
library(CalCEN)

parameters <- CalCEN::load_parameters()

runs <- readr::read_tsv(
    "product/runs_20210803.tsv",
    show_col_types = FALSE)

load("intermediate_data/estimated_expression.Rdata")
load("intermediate_data/estimated_expression_meta.Rdata")

#####################################
# check that all runs have metadata #
#####################################
run_progress <- runs |>
    dplyr::left_join(
        estimated_expression |>
        dplyr::distinct(run_accession) |>
        dplyr::mutate(got_expression = TRUE),
        by = c("run_accession")) |>
    dplyr::left_join(
        estimated_expression_meta |>
        dplyr::distinct(run_accession) |>
        dplyr::mutate(got_meta = TRUE),
        by = c("run_accession"))

run_progress |>
    dplyr::filter(got_expression, is.na(got_meta))


###########################################################
# Filter runs by number of genes with non-zero expression #
###########################################################

n_zero_genes <- estimated_expression |>
    dplyr::group_by(study_accession, run_accession) |>
    dplyr::summarize(
        n_nonzero_expression = sum(FPKM != 0),
        .groups = "drop")


runs_final <- runs |>
    dplyr::semi_join(
        runs |>
        dplyr::count(study_accession) |>
        dplyr::filter(n >= 20),
        by = "study_accession") |>
    dplyr::semi_join(
        n_zero_genes |>
        dplyr::filter(n_nonzero_expression >= 3113), # half the genes
        by = c("study_accession", "run_accession"))

save(runs_final, file = "intermediate_data/runs_final.Rdata")

runs_final |>
    readr::write_tsv("product/runs_20201024.tsv")



###########################################################################
# get the minimum fraction of reads aligned exactly 1 time for each study #
###########################################################################
n_reads_aligned <- runs |>
    dplyr::mutate(layout = ifelse(is_paired, "Paired", "Single")) |>
    dplyr::select(study_accession, run_accession, layout) |>
    dplyr::left_join(
        estimated_expression_meta,
        by = c("study_accession", "run_accession")) |>
    #dplyr::filter(n_reads |> is.na() |> magrittr::not()) |>
    dplyr::mutate(frac_exact_align = n_reads_aligned_1_time / n_reads) |>
    dplyr::arrange(frac_exact_align) |>
    dplyr::select(
        study_accession,
        run_accession,
        layout,
        n_reads_aligned_1_time,
        frac_exact_align)

p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 180)) +
    ggplot2::geom_dotplot(
        data = n_reads_aligned,
        mapping = ggplot2::aes(
            x = study_accession,
            y = frac_exact_align),
        binaxis = "y", stackdir = "center",
        stackratio = .4,
        dotsize = .3) +
    ggplot2::facet_wrap(
        ~layout,
        nrow = 2,
        scales = "free_x") +
    ggplot2::scale_x_discrete(name = "SRA Study Accession") +
    ggplot2::scale_y_continuous(
        name = "Percent reads mapped exactly once to the reference genome",
        label = scales::percent)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/percent_mapped_once_by_study_", CalCEN::date_code(), ".pdf"),
    height = 5, width = 10,
    useDingbats = FALSE)
    
ggplot2::ggsave(
    filename = paste0(
        "product/figures/percent_mapped_once_by_study_", CalCEN::date_code(), ".png"),
    height = 5, width = 10)

n_reads_aligned |>
    readr::write_tsv(
        file = paste0(
            "product/figures/percent_mapped_once_by_study_source_data_",
            CalCEN::date_code(), ".tsv"))

###################################################
# N zero expression genes vs percent reads mapped #
###################################################


n_zero_genes <- estimated_expression |>
    dplyr::group_by(study_accession, run_accession) |>
    dplyr::summarize(
        n_nonzero_expression = sum(FPKM != 0),
        .groups = "drop")

data <- n_zero_genes |>
    dplyr::left_join(
        n_reads_aligned,
        by = c("study_accession", "run_accession"))


ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::geom_point(
        data = data |> dplyr::filter(n_nonzero_expression < 5000),
        mapping = ggplot2::aes(
            x = n_reads_aligned_1_time + 1,
            y = n_nonzero_expression + 1),
        color = "grey50") +
    ggplot2::geom_point(
        data = data |> dplyr::filter(n_nonzero_expression >= 5000),
        mapping = ggplot2::aes(
            x = n_reads_aligned_1_time + 1,
            y = n_nonzero_expression + 1)) +
    ggplot2::ggtitle(
        label = "Expression depth vs coverage by RNA-seq run") +
    ggplot2::scale_x_log10(
        "Number reads that map exactly once to reference genome") +
    ggplot2::scale_y_log10(
        "Number of genes that have non-zero expression")

ggplot2::ggsave(
    filename = paste0(
        "product/figures/expression_depth_vs_coverage_by_study_",
        CalCEN::date_code(), ".pdf"),
    height = 8,
    width = 8,
    useDingbats = FALSE)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/expression_depth_vs_coverage_by_study_",
        CalCEN::date_code(), ".png"),
    height = 8,
    width = 8)

data |>
    readr::write_tsv(
        file = paste0(
            "product/figures/expression-depth_vs_coverage_by_study_",
            CalCEN::date_code(), ".tsv"))
