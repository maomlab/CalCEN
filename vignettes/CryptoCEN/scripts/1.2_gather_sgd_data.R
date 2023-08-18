library(tidyverse)
library(seqinr)
library(CalCEN)

parameters <- CalCEN::load_parameters()

#####################
# gather SAC genome #
#####################
cmd <- paste0(
    "cd raw_data && ",
    "wget ", parameters$source_data$genome_sac$reference_genome_url, " && ",
    "tar xzvf ", parameters$source_data$genome_sac$reference_genome_url |>
        basename(), " && ",
    "rm -rf ", parameters$source_data$genome_sac$reference_genome_url |>
        basename(), " && ",
    "cd ", parameters$source_data$genome_sac$reference_genome_url |>
        basename() |> stringr::str_replace(".tgz", ""), " && ",
    "gunzip orf_trans_all_*.fasta.gz")
cat(cmd, "\n", sep = "")
system(cmd)


##########################
# gather SAC transcripts #
##########################
sac_transcripts <- parameters$source_data$genome_sac$transcript_fasta_path |>
    seqinr::read.fasta(seqtype = "AA")


#####################################
# gather SAC transcript annotations #
#####################################
sac_transcript_annotations <- sac_transcripts |>
    purrr::map_chr(~seqinr::getAnnot(.)) |>
    stringr::str_replace("ADE5,7", "ADE57") |>
    stringr::str_replace("DUR1,2", "DUR12") |>
    stringr::str_replace("ARG5,6", "ARG56") |>
    stringr::str_replace("(reverse complement), ([^,]+)", "\\1; \\2") |>
    stringr::str_replace_all(
        "(Chr .*[0-9]+), Genome Release", "\"\\1\", Genome Release") |>
    readr::read_csv(
        col_names = c(
            "identifiers",
            "genomic_location",
            "genome_release",
            "feature_type",
            "description")) |>
    tidyr::separate(
        col = identifiers,
        into = c("systematic_name", "standard_name", "sgd_id"),
        sep = " ") |>
    dplyr::mutate(
        is_reverse_complement = feature_type |>
            stringr::str_detect("reverse complement"),
        feature_type = feature_type |> stringr::str_replace("reverse complement; ", ""),
        systematic_name = systematic_name |> stringr::str_replace("^>", ""),
        sgd_id = sgd_id |> stringr::str_replace("^SGDID:", ""),
        genome_release = genome_release |> stringr::str_replace("Genome Release ", ""))

save(
    sac_transcript_annotations,
    file = "intermediate_data/sac_transcript_annotations.Rdata")

sac_transcript_annotations |>
    readr::write_tsv(file = "product/sac_transcript_annotations_20211106.tsv")


##################################
# gather SAC GO Term annotations #
##################################

cmd <- paste0(
    "wget -cO ", parameters$source_data$genome_sac$go_slim_mapping_path, " ",
    parameters$source_data$genome_sac$go_slim_mapping_url)
cat(cmd, "\n", sep = "")
system(cmd)


go_annotations_sac <- readr::read_tsv(
    file = parameters$source_data$genome_sac$go_slim_mapping_path,
    col_names = c(
        "ORF", "Gene", "SGDID", "GO_Aspect", "GO Slim term",
        "ID", "Feature type"),
    show_col_types = FALSE) |>
    dplyr::rename(
        term = `GO Slim term`,
        go_id = ID,
        ontology = GO_Aspect)
save(go_annotations_sac, file = "intermediate_data/go_annotations_sac.Rdata")

cmd <- paste0(
    "wget -cO ", parameters$source_data$genome_sac$go_terms_path, " ",
    parameters$source_data$genome_sac$go_terms_url)
cat(cmd, "\n", sep = "")
system(cmd)

go_terms_sac <- readr::read_tsv(
    file = parameters$source_data$genome_sac$go_terms_path,
    col_names = c("GOID", "GO_Term", "GO_Aspect", "GO_Term_Definition"),
    show_col_types = FALSE) |>
    dplyr::transmute(
        go_id = paste0(
            "GO:",
            GOID |> stringr::str_pad(width = 7, side = "left", pad = "0")),
        term = GO_Term,
        ontology = GO_Aspect,
        definition = GO_Term_Definition)

save(go_terms_sac, file = "intermediate_data/go_terms_sac.Rdata")
