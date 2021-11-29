

library(tidyverse)
library(ggplot2)
library(dutchmasters) # color palette
library(CalCEN)


load("intermediate_data/h99_transcript_annotations.Rdata")
load("intermediate_data/sac_transcript_annotations.Rdata")
load("intermediate_data/sac_to_uniprot.Rdata")
load("intermediate_data/H99_to_sac_best_hits.Rdata")
load("intermediate_data/H99_to_sac_orthoMCL.Rdata")
load("intermediate_data/sac_complexes.Rdata")
load("intermediate_data/coexp_intra_study_coexp_spearman.Rdata")


if (!dir.exists("product/crypto_coexp_of_sac_complexes")) {
    cat("Creating output directory 'product/crypto_coexp_of_sac_complexes' ...\n")
    dir.create("product/crypto_coexp_of_sac_complexes")
}


# convert from 1 row per complex, to 1 row per member
sac_complex_members <- sac_complexes %>%
    dplyr::transmute(
        complex_id = `#Complex ac`,
        member = `Identifiers (and stoichiometry) of molecules in complex`) %>%
    tidyr::separate_rows(member, sep = "[|]")

# filter down to the proteins
sac_complex_members <- sac_complex_members %>%
    # trim off misc parts of the member name
    dplyr::mutate(
        member = member %>%
            stringr::str_replace("[(][0-9]+[)]$", "") %>% # stochiometry
            stringr::str_replace("-PRO_[0-9]+$", "")) %>%  # pro-peptides
    # handle the ribosome CPX-1599, CPX-1601
    dplyr::mutate(
        member = member %>%
            stringr::str_replace("^[\\[]", "") %>%
            stringr::str_replace("[\\]]$", "")) %>%
    tidyr::separate_rows(member, sep = ",") %>%
    # handle non-protein members
    dplyr::filter(
        !(member %>% stringr::str_detect("^URS")),    # RNA (rnacentral.org id)
        !(member %>% stringr::str_detect("^EBI-")),   # mRNA
        !(member %>% stringr::str_detect("^CHEBI")),  # small molecule
        !(member %>% stringr::str_detect("^CPX-")),   # other complex
        member != "P12294") %>%                       # mitochondrially encoded proteins
    dplyr::distinct(
        complex_id, member) %>%                       # P49687 is in CPX-824 twice
    dplyr::rename(uniprot_accession = member)

# traverse complexes to crypto orthologs
sac_complex_members <- sac_complex_members %>%
    dplyr::left_join(
        sgd_to_uniprot,
        by = "uniprot_accession") %>%
    dplyr::left_join(
        sac_transcript_annotations %>%
            dplyr::select(
                sgd_id,
                systematic_name,
                standard_name,
                description),
        by = "sgd_id") %>%
    dplyr::rename(
        sac_uniprot_accession = uniprot_accession,
        sac_uniprot_entry = uniprot_entry,
        sac_systematic_name = systematic_name,
        sac_standard_name = standard_name,
        sac_description = description) %>%
    dplyr::left_join(
        H99_to_sac_orthoMCL %>%
        dplyr::transmute(
            sac_systematic_name = scer %>% unlist(),
            h99_gene_id_orthoMCL = cneq %>% unlist()),
        by = "sac_systematic_name") %>%
    dplyr::left_join(
        H99_to_sac_best_hits %>%
        dplyr::transmute(
            sac_standard_name,
            h99_gene_id_blast = gene_id),
        by = "sac_standard_name")

# Where does blast and orthoMCL orthology differ?
# => 8 genes
# I don't have a great way to deside, so take the OrthoMCL
sac_complex_members %>%
    dplyr::left_join(
        h99_transcript_annotations %>%
        dplyr::select(
            h99_gene_id_orthoMCL = gene_id,
            h99_description_orthoMCL = description),
        by = "h99_gene_id_orthoMCL") %>%
    dplyr::left_join(
        h99_transcript_annotations %>%
        dplyr::select(
            h99_gene_id_blast = gene_id,
            h99_description_blast = description),
        by = "h99_gene_id_blast") %>%
    dplyr::filter(
        !is.na(h99_gene_id_orthoMCL),
        !is.na(h99_gene_id_blast),
        h99_gene_id_orthoMCL != h99_gene_id_blast) %>%
    dplyr::distinct(sac_systematic_name, .keep_all = TRUE) %>%
    dplyr::transmute(
        sac_systematic_name,
        sac_description = sac_description %>% stringr::str_sub(1, 60),
        h99_description_orthoMCL = h99_description_orthoMCL %>% stringr::str_sub(1, 60),
        h99_description_blast = h99_description_blast %>% stringr::str_sub(1, 60)) %>%
    data.frame()

sac_complex_members <- sac_complex_members %>%
    dplyr::mutate(
        h99_gene_id = ifelse(
            !is.na(h99_gene_id_orthoMCL),
            h99_gene_id_orthoMCL,
            h99_gene_id_blast)) %>%
    dplyr::select(
        -h99_gene_id_orthoMCL,
        -h99_gene_id_blast) %>%
    dplyr::left_join(
        h99_transcript_annotations %>%
        dplyr::select(
            h99_gene_id = gene_id,
            h99_description = description),
        by = "h99_gene_id")
    
# convert from 1-row per member to 1-row per pair members in a complex
sac_complex_interactions <- sac_complex_members %>%
    dplyr::rename(
        sac_uniprot_accession_1 = sac_uniprot_accession,
        sac_uniprot_entry_1 = sac_uniprot_entry,
        sgd_id_1 = sgd_id,
        sac_systematic_name_1 = sac_systematic_name,
        sac_standard_name_1 = sac_standard_name,
        sac_description_1 = sac_description,
        h99_gene_id_1 = h99_gene_id,
        h99_description_1 = h99_description) %>%
    dplyr::left_join(
        sac_complex_members %>%
            dplyr::rename(
                sac_uniprot_accession_2 = sac_uniprot_accession,
                sac_uniprot_entry_2 = sac_uniprot_entry,
                sgd_id_2 = sgd_id,
                sac_systematic_name_2 = sac_systematic_name,
                sac_standard_name_2 = sac_standard_name,
                sac_description_2 = sac_description,
                h99_gene_id_2 = h99_gene_id,
                h99_description_2 = h99_description),
        by = c("complex_id")) %>%
    dplyr::filter(
        sac_uniprot_accession_1 != sac_uniprot_accession_2)

# prepend info about each complex
sac_complex_interactions <- sac_complexes %>%
    dplyr::select(
        complex_id = `#Complex ac`,
        complex_name = `Recommended name`,
        complex_aliases = `Aliases for complex`,
        complex_description  = Description) %>%
    dplyr::left_join(
        sac_complex_interactions,
        by = "complex_id")
    
# For each pair of genes,
# filter down to the best scoring pair of transcripts
coexp_genes <- coexp_intra_study_coexp_spearman %>%
    dplyr::left_join(
        h99_transcript_annotations %>%
        dplyr::select(
            cnag_id_1 = cnag_id,
            gene_id_1 = gene_id),
        by = "cnag_id_1") %>%
    dplyr::left_join(
        h99_transcript_annotations %>%
        dplyr::select(
            cnag_id_2 = cnag_id,
            gene_id_2 = gene_id),
        by = "cnag_id_2") %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    dplyr::distinct(
        gene_id_1,
        gene_id_2,
        .keep_all = TRUE) %>%
    dplyr::select(
        h99_gene_id_1 = gene_id_1,
        h99_gene_id_2 = gene_id_2,
        coexp_score = score)

sac_complex_interactions <- sac_complex_interactions %>%
    dplyr::left_join(
        coexp_genes,
        by = c("h99_gene_id_1", "h99_gene_id_2"))
sac_complex_interactions %>%
    readr::write_tsv(
        file = paste0(
            "product/crypto_coexp_of_sac_complexes/sac_complex_interactions_",
            CalCEN::date_code(), ".tsv"))
save(
    sac_complex_interactions,
    file = "intermediate_data/sac_complex_interactions.Rdata")


sac_complex_summary <- sac_complex_interactions %>%
    dplyr::group_by(
        complex_id,
        complex_name,
        complex_description) %>%
    dplyr::summarize(
        n_sac_members = sac_uniprot_accession_1 %>%
            unique() %>%
            length(),
        n_H99_members = sum(!is.na(unique(h99_gene_id_1))),
        mean_coexp_score = mean(coexp_score, na.rm = TRUE)) %>%
    dplyr::ungroup()
sac_complex_summary %>%
    readr::write_tsv(
        file = paste0(
            "product/crypto_coexp_of_sac_complexes/sac_complex_summary_",
            CalCEN::date_code(), ".tsv"))


################################################

# compare distribution of co-expression across co-complex pairs
# relative to all pairs
plot <- ggplot2::ggplot(
    data = sac_complex_interactions %>%
        dplyr::filter(
            !is.na(coexp_score),
            gene_id_1 < gene_id_2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_density(
        mapping = ggplot2::aes(
            x = coexp_score),
        fill = dutchmasters$pearl_earring[1]) +
    ggplot2::ggtitle(
        "Distribution of co-expression of co-complex genes") +
    ggplot2::scale_x_continuous(
        "Co-expression rank",
        limits = c(0, 1))

ggplot2::ggsave(
    filename = paste0(
        "product/crypto_coexp_of_sac_complexes/coexp_distribution_of_co-complex_genes_",
        CalCEN::date_code(), ".pdf"),
    plot = plot,
    width = 6,
    height = 4)

ggplot2::ggsave(
    filename = paste0(
        "product/crypto_coexp_of_sac_complexes/coexp_distribution_of_co-complex_genes_",
        CalCEN::date_code(), ".png"),
    plot = plot,
    width = 6,
    height = 4)

########################
#
#

# n complexes
# 616
sac_complexes %>% nrow()

# number of h99 associations among complex members
# 13950
sac_complex_interactions %>%
    dplyr::filter(!is.na(h99_gene_id_1), !is.na(h99_gene_id_2)) %>%
    dplyr::distinct(h99_gene_id_1, h99_gene_id_2) %>%
    nrow()


# number of associations among ribosomal proteins
# 7142
sac_complex_interactions %>%
    dplyr::filter(!is.na(h99_gene_id_1), !is.na(h99_gene_id_2)) %>%
    dplyr::filter(
        complex_name %>% stringr::str_detect("ribosomal") |
        complex_name %>% stringr::str_detect("ribonucleoprotein")) %>%
    dplyr::distinct(h99_gene_id_1, h99_gene_id_2) %>%
    nrow()

# number of ribosomal complexes
# 17
sac_complex_interactions %>%
    dplyr::filter(!is.na(h99_gene_id_1), !is.na(h99_gene_id_2)) %>%
    dplyr::filter(
        complex_name %>% stringr::str_detect("ribosomal") |
        complex_name %>% stringr::str_detect("ribonucleoprotein")) %>%
    dplyr::distinct(complex_id) %>%
    nrow()
