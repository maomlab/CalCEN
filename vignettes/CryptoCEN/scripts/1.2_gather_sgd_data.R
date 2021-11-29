
library(tidyverse)
library(seqinr)
library(CalCEN)

parameters <- CalCEN::load_parameters()


cmd <- paste0(
    "cd raw_data && ",
    "wget ", parameters$source_data$genome_sac$reference_genome_url, " && ",
    "tar xzvf ", parameters$source_data$genome_sac$reference_genome_url %>% basename(), " && ",
    "rm -rf ", parameters$source_data$genome_sac$reference_genome_url %>% basename()), " && ",
    "cd ", parameters$source_data$genome_sac$reference_genome_url %>% basename() %>% stringr::str_replace(".tgz", ""), " && ",
    "gunzip orf_trans_all_*.fasta.gz")
cat(cmd, "\n", sep = "")
system(cmd)

sac_transcripts <- parameters$source_data$genome_sac$transcript_fasta_path %>%
    seqinr::read.fasta(seqtype = "AA")

sac_transcript_annotations <- sac_transcripts %>%
    purrr::map_chr(~seqinr::getAnnot(.)) %>%
    stringr::str_replace("ADE5,7", "ADE57") %>%
    stringr::str_replace("DUR1,2", "DUR12") %>%
    stringr::str_replace("ARG5,6", "ARG56") %>%
    stringr::str_replace("(reverse complement), ([^,]+)", "\\1; \\2") %>%
    stringr::str_replace_all("(Chr .*[0-9]+), Genome Release", "\"\\1\", Genome Release") %>%
    readr::read_csv(
        col_names = c(
            "identifiers",
            "genomic_location",
            "genome_release",
            "feature_type",
            "description")) %>%
    tidyr::separate(
        col = identifiers,
        into = c("systematic_name", "standard_name", "sgd_id"),
        sep = " ") %>%
    dplyr::mutate(
        is_reverse_complement = feature_type %>% stringr::str_detect("reverse complement"),
        feature_type = feature_type %>% stringr::str_replace("reverse complement; ", ""),
        systematic_name = systematic_name %>% stringr::str_replace("^>", ""),
        sgd_id = sgd_id %>% stringr::str_replace("^SGDID:", ""),
        genome_release = genome_release %>% stringr::str_replace("Genome Release ", ""))
    
save(
    sac_transcript_annotations,
    file = "intermediate_data/sac_transcript_annotations.Rdata")

sac_transcript_annotations %>%
    readr::write_tsv(file = "product/sac_transcript_annotations_20211106.tsv")
