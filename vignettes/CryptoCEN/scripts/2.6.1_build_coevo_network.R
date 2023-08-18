library(tidyverse)
library(EGAD)

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id


load("intermediate_data/Cneo_genes_per_og.Rdata")
load("intermediate_data/Cneo_coevo_res.Rdata")
load("intermediate_data/Cneo_coevo_gene_identifiers.Rdata")



coevolution_network_long <- Cneo_coevo_res |>
    dplyr::inner_join(
        Cneo_coevo_gene_identifiers |>
        dplyr::transmute(
            coevo_gene_id_1 = coevo_gene_id,
            gene_id_1 = gene_id),
        by = "coevo_gene_id_1") |>
    dplyr::inner_join(
        Cneo_coevo_gene_identifiers |>
        dplyr::transmute(
            coevo_gene_id_2 = coevo_gene_id,
            gene_id_2 = gene_id),
        by = "coevo_gene_id_2") |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::distinct(gene_id, .keep_all = TRUE) |>
        dplyr::select(
            gene_id_1 = gene_id,
            cnag_id_1 = cnag_id),
        by = "gene_id_1") |>
    dplyr::left_join(
        h99_transcript_annotations |>
        dplyr::distinct(gene_id, .keep_all = TRUE) |>
        dplyr::select(
            gene_id_2 = gene_id,
            cnag_id_2 = cnag_id),
        by = "gene_id_2")

save(
    coevolution_network_long,
    file = "intermediate_data/coevolution_network_long.Rdata")
coevolution_network_long |>
    readr::write_tsv(
        file = "product/coevolution_network_long_20230712.tsv")


# check that the gene-gene interactions are 1-1
coevolution_network_long |>
    dplyr::count(cnag_id_1, cnag_id_2, sort = TRUE) |>
    dplyr::filter(n > 1)

coevolution_network <- coevolution_network_long |>
  dplyr::transmute(
    cnag_id_1,
    cnag_id_2,
    weight = (-coevolution_coefficient + 1) / 2) |>
  as.matrix() |>
  EGAD::build_weighted_network(
    list = genes)

coevolution_network[is.na(coevolution_network)] <- 0
    

save(coevolution_network, file = "intermediate_data/coevolution_network.Rdata")

# over 5264 genes
# 27511028 non-zero entries
# 32.58144% sparse (27511028 / (9189 * 9189))


coevolution_network <- coevolution_network_long |>
  dplyr::transmute(
    cnag_id_1,
    cnag_id_2,
    weight = p_value < 1e-5) |>
  as.matrix() |>
  EGAD::build_binary_network(list = genes) |>
  EGAD::extend_network()

coevolution_network[is.na(coevolution_network)] <- 0


save(coevolution_network, file = "intermediate_data/coevolution_network.Rdata")

# over 5264 genes
# 27511028 non-zero entries
# 32.58144% sparse (27511028 / (9189 * 9189))
