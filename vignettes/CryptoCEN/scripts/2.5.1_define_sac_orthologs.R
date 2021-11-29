library(plyr)
library(tidyverse)
library(future)
library(future.batchtools)
library(CalCEN)
library(Bethany)
library(EGAD)

future::plan(
    future.batchtools::batchtools_slurm(
        template = system.file(
            "hpc",
            "batchtools.greatlakes.tmpl",
            package = "CalCEN")),
    resources = list(
        ncpus = 20,
        memory = 5000))


parameters <- CalCEN::load_parameters()


load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id


# use the future package to distribute computing the H99-by-sac blastp
tasks <- tibble::tibble(
    sequence = seqinr::read.fasta(
        file = parameters$source_data$genome$protein_fasta_path,
        seqtype = "AA")) %>%
    dplyr::mutate(
        batch = rep(1:100, length.out = dplyr::n())) %>%
    dplyr::group_by(batch) %>%
    dplyr::do({
        batch <- .$batch[1]
        query <- .$sequence
        cat("Submitting blastp task for batch '", batch, "' querying '", length(query), "' H99 genes against the Sac genome\n", sep = "")
        blastp_network <- future::future({
            Bethany::blastp(
                ref = parameters$source_data$genome_sac$protein_fasta_path,
                query = query,
                run_id = paste0("H99-vs-sac_", batch),
                blastp_num_threads = 20,
                verbose = TRUE)
        })
        tibble::tibble(
            job = list(blastp_network))
    })

# wait till all the jobs have completed in the queue
H99_to_sac_blastp <- tasks$job %>%
    purrr::map_dfr(function(task) {
        network <- task %>% future::value()
    })
save(H99_to_sac_blastp, file = "intermediate_data/H99_to_sac_blastp.Rdata")

######
# use the future package to distribute computing the sac-by-H99 blastp
tasks <- tibble::tibble(
    sequence = seqinr::read.fasta(
        file = parameters$source_data$genome_sac$protein_fasta_path,
        seqtype = "AA")) %>%
    dplyr::mutate(
        batch = rep(1:100, length.out = dplyr::n())) %>%
    dplyr::group_by(batch) %>%
    dplyr::do({
        batch <- .$batch[1]
        query <- .$sequence
        cat("Submitting blastp task for batch '", batch, "' querying '", length(query), "' Sac genes against the H99 genome\n", sep = "")
        blastp_network <- future::future({
            Bethany::blastp(
                ref = parameters$source_data$genome$protein_fasta_path,
                query = query,
                run_id = paste0("H99-vs-sac_", batch),
                blastp_num_threads = 20,
                verbose = TRUE)
        })
        tibble::tibble(
            job = list(blastp_network))
    })

# wait till all the jobs have completed in the queue
sac_to_H99_blastp <- tasks$job %>%
    purrr::map_dfr(function(task) {
        network <- task %>% future::value()
    })
save(sac_to_H99_blastp, file = "intermediate_data/sac_to_H99_blastp.Rdata")


####

H99_to_sac_best_hits <- H99_to_sac_blastp %>%
    dplyr::arrange(EValue) %>%
    dplyr::group_by(query_target) %>%
    dplyr::slice(1) %>%
    dplyr::rename(
        H99_target = ref_target,
        sac_target = query_target,
        H99_to_sac_bit_score = bit_score,
        H99_to_sac_EValue = EValue) %>%
    dplyr::inner_join(
        sac_to_H99_blastp %>%
        dplyr::arrange(EValue) %>%
        dplyr::group_by(query_target) %>%
        dplyr::slice(1) %>%
        dplyr::rename(
            sac_target = ref_target,
            H99_target = query_target,
            sac_to_H99_bit_score = bit_score,
            sac_to_H99_EValue = EValue),
        by = c("sac_target", "H99_target")) %>%
    dplyr::filter(
        H99_to_sac_EValue < 1e-5,
        sac_to_H99_EValue < 1e-5) %>%
    dplyr::mutate(
        cnag_id = H99_target %>% stringr::str_replace("-p[0-9]+$", ""))

load("intermediate_data/h99_transcript_annotations.Rdata")

H99_to_sac_best_hits <- H99_to_sac_best_hits %>%
    dplyr::inner_join(
        h99_transcript_annotations,
        by = c("cnag_id"))

load("intermediate_data/sac_transcript_annotations.Rdata")
H99_to_sac_best_hits <- H99_to_sac_best_hits %>%
    dplyr::inner_join(
        sac_transcript_annotations %>%
        dplyr::select(
            systematic_name,
            sac_standard_name = standard_name,
            sac_feature_type = feature_type,
            sac_description = description),
        by = c("sac_target" = "systematic_name")) %>%
    dplyr::ungroup()


save(H99_to_sac_best_hits, file = "intermediate_data/H99_to_sac_best_hits.Rdata")

H99_to_sac_best_hits %>%
    readr::write_tsv("product/H99_to_sac_best_hits_20211106.tsv")
