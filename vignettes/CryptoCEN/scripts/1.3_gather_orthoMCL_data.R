library(tidyverse)
library(CalCEN)

parameters <- CalCEN::load_parameters()

cmd <- paste0(
    "mkdir -p raw_data/OrthoMCL && ",
    "cd raw_data/OrthoMCL && ",
    "wget ", parameters$source_data$OrthoMCL$ortho_groups_url, " && ",
    "gunzip ", parameters$source_data$OrthoMCL$ortho_groups_url |>
        basename())
cat(cmd, "\n", sep = "")
system(cmd)

ortho_groups <- readr::read_delim(
    parameters$source_data$OrthoMCL$ortho_groups_path,
    delim = ":",
    col_names = c(
        "group_id",
        "accession"),
    col_types = readr::cols(
        group_id = readr::col_character(),
        accession = readr::col_character()),
    trim_ws = TRUE) |>
    dplyr::filter(
        accession |> stringr::str_detect("cneq[|]")) |>
    tidyr::separate_rows(accession, sep = " ") |>
    dplyr::transmute(
        group_id,
        species_mnemonic = accession |> stringr::str_extract("^[^|]+"),
        protein_accession = accession |> stringr::str_extract("[^|]+$"))

# cneq: Cryptococcus neoformans var. grubii H99
# scer: Saccharomyces cerevisiae S288C
# calb: Candida albicans SC5314

H99_orthologs <- ortho_groups |>
    dplyr::filter(species_mnemonic %in% c("cneq", "scer", "calb")) |>
    tidyr::pivot_wider(
        names_from = "species_mnemonic",
        values_from = "protein_accession",
        values_fn = list)

save(H99_orthologs, file = "intermediate_data/H99_orthologs.Rdata")


H99_to_sac_orthoMCL <- H99_orthologs |>
    dplyr::rowwise() |>
    dplyr::mutate(n_cneq = length(cneq), n_scer = length(scer)) |>
    dplyr::ungroup() |>
    dplyr::filter(n_cneq == 1, n_scer == 1)

save(
    H99_to_sac_orthoMCL,
    file = "intermediate_data/H99_to_sac_orthoMCL.Rdata")
