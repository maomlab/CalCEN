
library(tidyverse)
library(CalCEN)


runs <- readr::read_tsv("product/runs_20210803.tsv")

load("intermediate_data/h99_transcript_annotations.Rdata")


######## Collect Expression data ##############
estimated_expression <- runs |>
    dplyr::select(study_accession, run_accession, is_paired) |>
    dplyr::left_join(
        CalCEN::gather_estimated_expression(
            results_dir = "intermediate_data/estimated_expression_20210804"),
        by = "run_accession") |>
    dplyr::rename(cnag_id = gene_id)
 
# gene_id and transcript_ids are 1-1
estimated_expression |> BioChemPantry::summarize_map(
    x_cols = cnag_id,
    y_cols = transcript_ids)

estiamted_expression <- estimated_expression |>
    dplyr::select(-transcript_ids)

# are all the expressed genes in the h99_transcript_annotations?
estimated_expression |>
    dplyr::distinct(cnag_id) |>
    dplyr::mutate(in_expression = TRUE) |>
    dplyr::full_join(
        h99_transcript_annotations |>
        dplyr::transmute(
            cnag_id,
            in_transcripts = TRUE),
        by = "cnag_id") |>
    dplyr::count(in_expression, in_transcripts)


save(
    estimated_expression,
    file = "intermediate_data/estimated_expression.Rdata")

estimated_expression |>
    readr::write_tsv("product/estimated_expression_20210807.tsv")


######## Collect Meta data ###################
estimated_expression_meta <- runs |>
    dplyr::select(study_accession, run_accession, is_paired) |>
    dplyr::left_join(
        CalCEN::gather_estimated_expression_meta(
            results_dir = "intermediate_data/estimated_expression_20210804"),
        by = "run_accession")

save(
    estimated_expression_meta,
    file = "intermediate_data/estimated_expression_meta.Rdata")

estimated_expression_meta |>
    readr::write_tsv("product/estimated_expression_meta_20210807.tsv")
