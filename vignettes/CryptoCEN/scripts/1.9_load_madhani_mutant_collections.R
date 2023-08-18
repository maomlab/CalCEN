library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()

load("intermediate_data/h99_transcript_annotations.Rdata")

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
   
