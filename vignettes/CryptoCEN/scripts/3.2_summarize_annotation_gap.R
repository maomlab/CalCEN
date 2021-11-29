


library(tidyverse)
library(ggplot2)
library(CalCEN)

parameters <- CalCEN::load_parameters()

# H99 genes
# "hypothetical product" or "unspecified product" in description
# "sac ortholog"
#   reciprical best
#   OrthoMCL
# "C. albicans ortholog"
#   reciprical best
#   OrthoMCL
# "Biological Process" annotation
#   propagate
#   filter out "NOT" qualified annotations
#   filter for terms with < 1000 annotations
#   
# CryptoNet Edge
# Co-expression cluster
# BioGRID PPI
# Co-evolution



load("intermediate_data/h99_transcript_annotations.Rdata")
load("intermediate_data/go_annotations_propagated.Rdata")
load("intermediate_data/H99_to_sac_best_hits.Rdata")
load("intermediate_data/H99_to_sac_orthoMCL.Rdata")
load("intermediate_data/CryptoNetV1.Rdata")

UMAP_genes_cluster_labels <- readr::read_tsv(
    "product/figures/coexp_embedding/UMAP_genes_cluster_labels_20210815.tsv")



data <- h99_transcript_annotations %>%
    dplyr::filter(!is_pseudo) %>%
    # is orphan
    dplyr::transmute(
        cnag_id,
        gene_id,
        is_orphan =
            (description %>% stringr::str_detect("hypothetical protein")) |
            (description %>% stringr::str_detect("unspecified product"))) %>%
    # un annotated
    dplyr::left_join(
        go_annotations_propagated %>%
            dplyr::rename(gene_id = db_object_id) %>%
            dplyr::count(gene_id, ontology) %>%
            tidyr::pivot_wider(
                id_cols = "gene_id",
                names_prefix = "count_",
                names_from = "ontology",
                values_from = "n",
                values_fill = 0),
        by = "gene_id") %>%
    dplyr::mutate(across(
        .cols = tidyselect::starts_with("count_"),
        ~ifelse(is.na(.), 0, .))) %>%
    # no sac ortholog blast
    dplyr::left_join(
        H99_to_sac_best_hits %>%
            dplyr::transmute(
                gene_id,
                has_sac_ortholog_blast = !is.na(sac_standard_name)),
        by = "gene_id") %>%
    dplyr::mutate(has_sac_ortholog_blast = !is.na(has_sac_ortholog_blast)) %>%
    # no sac ortholog orthoMCL
    dplyr::left_join(
        H99_to_sac_orthoMCL %>%
        dplyr::transmute(
            gene_id = cneq %>% unlist(),
            has_sac_ortholog_orthoMCL = TRUE),
        by = "gene_id") %>%
    dplyr::mutate(has_sac_ortholog_orthoMCL = !is.na(has_sac_ortholog_orthoMCL)) %>%
    #in CryptoNet
    dplyr::left_join(
        CryptoNetV1 %>%
        dplyr::distinct(target1) %>%
        dplyr::transmute(
            gene_id = target1,
            in_CryptoNetV1 = TRUE),
        by = "gene_id") %>%
    dplyr::mutate(in_CryptoNetV1 = !is.na(in_CryptoNetV1)) %>%
    # expression cluster label
    dplyr::left_join(
        UMAP_genes_cluster_labels %>%
        dplyr::transmute(
            cnag_id,
            cluster_label),
        by = "cnag_id")

# how many gene products
data %>% nrow()
# 9185

# how many genes
data %>% dplyr::distinct(gene_id) %>% nrow
# 8334

# totally un-annotated
data %>%
    dplyr::distinct(
        gene_id,
        .keep_all = TRUE) %>%
    dplyr::filter(
        is_orphan,
        count_BP == 0,
        !has_sac_ortholog_blast,
        !has_sac_ortholog_orthoMCL,
        !in_CryptoNetV1) %>%
    nrow()
# 1359

data %>%
    dplyr::distinct(
        gene_id,
        .keep_all = TRUE) %>%
    dplyr::filter(
        !is_orphan) %>%
    dplyr::count(
        count_BP == 0,
        !has_sac_ortholog_blast,
        !has_sac_ortholog_orthoMCL,
        !in_CryptoNetV1)        

# Genes                                8334
#  Has BP                              2662
#  Sac ortholog                        1481
#  In CryptoNetV1                      1488
#  Has gene description                1344
#  Totally un-annotated                1359

#  BP Sac CryptoNet Name
#  Y
#  N  Y
#  N  N   Y
#  N  N   N         Y
#  N  N   N         N
