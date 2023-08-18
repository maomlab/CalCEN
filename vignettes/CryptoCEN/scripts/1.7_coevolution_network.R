
library(tidyverse)

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id


# Cneo genes per orthologous group.  The first column is the
# orthologous group identifier and the second column is the Cneo gene.
Cneo_genes_per_og <- readr::read_tsv(
    "raw_data/coevolution_networks/Cneo_genes_per_og.txt",
    col_names = c(
        "orthologous_group",
        "gene_id"))
save(Cneo_genes_per_og, file = "intermediate_data/Cneo_genes_per_og.Rdata")

# I have also attached the coevolution results with the same format as
# before wherein columns 1 and 2 are gene identifiers and columns 3
# and 4 are the coevolution coefficients and the p-value
# (uncorrected).

Cneo_coevo_res <- readr::read_tsv(
    "raw_data/coevolution_networks/Cneo_coevo_res.txt",
    col_names = c(
        "coevo_gene_id_1",
        "coevo_gene_id_2",
        "coevolution_coefficient",
        "p_value"))
save(Cneo_coevo_res, file = "intermediate_data/Cneo_coevo_res.Rdata")


# To help navigate which gene identifiers correspond to
# which Cneo gene, I generated a Cneo_coevo_gene_identifiers.txt file
# with this information. In Cneo_coevo_gene_identifiers.txt, the first
# column is the orthogroup identifier, followed by the CNAG gene
# identifier, followed by the NCBI identifier.

Cneo_coevo_gene_identifiers <- readr::read_tsv(
    "raw_data/coevolution_networks/Cneo_coevo_gene_identifiers.txt",
    col_names = c(
        "coevo_gene_id",
        "gene_id",
        "ncbi_id"))
save(
    Cneo_coevo_gene_identifiers,
    file = "intermediate_data/Cneo_coevo_gene_identifiers.Rdata")
