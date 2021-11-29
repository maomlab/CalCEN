
library(tidyverse)
library(seqinr)
library(CalCEN)

parameters <- CalCEN::load_parameters()

#######
# General Feature Format (GFF) file
# a simple tab-delimited text file for describing genomic features

cmd <- paste0(
    "cd raw_data && ",
    "wget ", parameters$source_data$genome$gff_annotations_url)
cat(cmd, "\n", sep = "")
system(cmd)

gff_annotations <- readr::read_tsv(
    file = parameters$source_data$genome$gff_annotations_path,
    col_names = c(
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute"),
    skip = 17) %>%
    dplyr::filter(feature == "gene") %>%
    dplyr::mutate(
        gene_id = attribute %>%
            stringr::str_extract("ID=[^;]+") %>%
            stringr::str_replace("ID=", ""),
        description = attribute %>%
            stringr::str_extract("description=[^;]+") %>%
            stringr::str_replace("description=", ""),
        .before = 1) %>%
    dplyr::select(-attribute, -frame)

save(gff_annotations, file = "intermediate_data/gff_annotations.Rdata")

# H99 transcripts
cmd <- paste0(
    "cd raw_data && ",
    "wget ", parameters$source_data$genome$transcript_fasta_url)
cat(cmd, "\n", sep = "")
system(cmd)


h99_transcripts <- seqinr::read.fasta(
    file = parameters$source_data$genome$transcript_fasta_path,
    seqtype = "AA")
    
h99_transcript_annotations <- h99_transcripts %>%
    purrr::map_chr(~seqinr::getAnnot(.)) %>%
    data.frame(annotation = .) %>%
    tidyr::separate(
        col = annotation,
        into = c("cnag_id", "gene", "organism", "gene_product", "transcript_product", "location", "length", "sequence_SO", "SO", "is_pseudo"),
        sep = " [|] ") %>%
    dplyr::mutate(
        dplyr::across(
            .cols = everything(),
            ~stringr::str_replace(., "^[a-zA-Z_]+=", ""))) %>%
    dplyr::mutate(
        cnag_id = cnag_id %>% stringr::str_replace("^>", ""),
        variant = cnag_id %>% stringr::str_extract("t[0-9]+_[0-9]+$"),
        gene_id = cnag_id %>% stringr::str_replace("-t[0-9]+_[0-9]+$", ""),
        is_pseudo = ifelse(is_pseudo == "true", TRUE, FALSE)) %>%
    dplyr::select(-gene)

h99_transcript_annotations <- h99_transcript_annotations %>%
    dplyr::left_join(
        gff_annotations,
        by = "gene_id")


# get gene-names from fungidb.org
gene_names <- readr::read_tsv("raw_data/gene_names_20210830.tsv") %>%
    dplyr::select(
        cnag_id = source_id,
        gene_symbol = `Gene Name or Symbol`) %>%
    dplyr::mutate(
        gene_symbol = ifelse(gene_symbol == "N/A", NA, gene_symbol))

h99_transcript_annotations <- h99_transcript_annotations %>%
    dplyr::left_join(gene_names, by = "cnag_id")

save(
    h99_transcript_annotations,
    file = "intermediate_data/h99_transcript_annotations.Rdata")

h99_transcript_annotations %>%
    readr::write_tsv(file = "product/h99_transcript_annotations_20210724.tsv")


####
# Annotated Protein sequences
cmd <- paste0(
    "cd raw_data && ",
    "wget ", parameters$source_data$genome$protein_fasta_url)
cat(cmd, "\n", sep = "")
system(cmd)



#######
# General Feature Format (GFF) file
# a simple tab-delimited text file for describing genomic features

cmd <- paste0(
    "cd raw_data && ",
    "wget ", parameters$source_data$genome$gaf_annotations_url)
cat(cmd, "\n", sep = "")
system(cmd)

go_annotations <- readr::read_tsv(
    file = parameters$source_data$genome$gaf_annotations_path,
    skip = 1,
    col_names = c(
        "db",
        "db_object_id",
        "db_object_symbol",
        "qualifier",
        "go_id",
        "db_reference_id",
        "evidence_code",
        "with_or_from",
        "aspect",
        "db_object_name",
        "db_object_synonym",
        "db_object_type",
        "taxon",
        "date",
        "assigned_by",
        "annotation_extension",
        "gene_product_form_id"))

save(go_annotations, file = "intermediate_data/go_annotations.Rdata")

