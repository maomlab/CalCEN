# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:

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


# use the future package to distribute computing the all-by-all blastp
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
				cat("Submitting blastp task for batch '", batch, "' querying '", length(query), "' genes against the genome\n", sep = "")
				blastp_network <- future::future({
						Bethany::blastp(
								ref = parameters$source_data$genome$protein_fasta_path,
								query = query,
								run_id = paste0("gene-vs-gene_", batch),
								blastp_num_threads = 20,
								verbose = TRUE)
				})
				tibble::tibble(
						job = list(blastp_network))
		})

# wait till all the jobs have completed in the queue
blastp <- tasks$job %>%
		purrr::map_dfr(function(task) {
				network <- task %>% future::value()
		})
save(blastp, file = "intermediate_data/blastp.Rdata")

############################################
# compute the blastp scores into a network #
############################################
blastp_network_long <- blastp %>%
		dplyr::transmute(
				# CNAG_00001-t26_1-p1 --> CNAG_00001-t26_1
				ref_target = ref_target %>% stringr::str_replace("-[^-]+$", ""),
				query_target = query_target %>% stringr::str_replace("-[^-]+$", ""),
				rank_score = rank(bit_score, na.last = "keep", ties.method = "average") / dplyr::n()) %>%
		dplyr::select(ref_target, query_target, rank_score) %>%
		as.data.frame()

blastp_network_long %>%
		as.data.frame() %>%
		readr::write_tsv(
				file = paste0("product/blastp_network_long_", CalCEN::date_code(), ".tsv"))

blastp_network <- blastp_network_long %>%
		EGAD::build_weighted_network(genes)

save(blastp_network, file = "intermediate_data/blastp_network.Rdata")

