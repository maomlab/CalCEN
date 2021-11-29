


#' Gather estimated expression data for analysis
#'
#' After estimating expression locally or in parallel
#' gather_estimated_expression parses files of the form
#'
#'      <results_dir>/<run_accession>.genes.results
#' 
#' @param results_dir path used in e.g. estimate_expression(...)
#' @param verbose (default: true)
#' @return a data.frame with columns (see RSEM for interpretation)
#'     run_accession
#'     gene_id
#'     transcript_ids
#'     length
#'     effective_length
#'     expected_count
#'     TPM
#'     FPKM
#' and one row per (run_accession, gene_id)
#' @export
gather_estimated_expression <- function(
    results_dir,
    verbose = TRUE) {
    if (!dir.exists(results_dir)) {
        cat("ERROR: results_dir '", results_dir, "', does not exist\n")
    }

    results_files <- list.files(
        path = results_dir,
        pattern = "*genes.results")
    
    if (verbose) {
        cat(
            "Gathering estimated expression for ", length(results_files), " runs in ",
            " '", results_dir, "' ...\n", sep = "")
    }
    
    estimated_expression <- results_files %>%
        plyr::ldply(function(results_fname) {
            results <- readr::read_tsv(
                file = paste0(results_dir, "/", results_fname),
                col_types = readr::cols(
                    gene_id = readr::col_character(),
                    `transcript_id(s)` = readr::col_character(),
                    length = readr::col_double(),
                    effective_length = readr::col_double(),
                    expected_count = readr::col_double(),
                    TPM = readr::col_double(),
                    FPKM = readr::col_double())) %>%
                dplyr::rename(transcript_ids = `transcript_id(s)`) %>%
                dplyr::mutate(
                    run_accession = results_fname %>% stringr::str_extract("^[^.]+"),
                    .before = 1)
        })
}
#'Gather estimated expression data metadata
#'
#' @param results_dir path used in e.g. estimate_expression(...)
#' @param verbose (default: true)
#' @return tibble::tibble() with rows representing runs and columns
#'     n_reads
#'     n_reads_aligned_0_times
#'     n_reads_aligned_1_time
#'     n_reads_aligned_mutiple_times
#'     overall_alignment_percent
#'     run_time
#'
#' 
#'@export
gather_estimated_expression_meta <- function(
    results_dir,
    verbose = TRUE) {

    if (!dir.exists(results_dir)) {
        cat("ERROR: results_dir '", results_dir, "', does not exist\n")
    }

    log_files <- list.files(
        path = paste0(results_dir, "/logs"),
        pattern = "*log")
    
    if (verbose) {
        cat(
            "Gathering estimated expression for ", length(log_files), " runs in ",
            " '", results_dir, "' ...\n", sep = "")
    }


    estimated_expression_meta <- log_files  %>%
	plyr::ldply(function(log_fname) {
            cat("processing log file: ", log_fname, " ...\n", sep = "")
            tryCatch({
                lines <- readr::read_lines(
                    file = paste0(results_dir, "/logs/", log_fname))
                n_reads <- lines %>%
                    stringr::str_detect(" reads; of these:") %>%
                    magrittr::extract(lines, .) %>%
                    stringr::str_extract("^[0-9]+") %>%
                    as.numeric()
                n_reads_aligned_0_times <- lines %>%
                    stringr::str_detect(" aligned[ a-z]* 0 times") %>%
                    magrittr::extract(lines, .) %>%
                    stringr::str_match("^ +([0-9]+)") %>%
                    unlist() %>%
                    magrittr::extract2(2) %>%
                    as.numeric()
                n_reads_aligned_1_time <- lines %>%
                    stringr::str_detect(" aligned[ a-z]* exactly 1 time") %>%
                    magrittr::extract(lines, .) %>%
                    stringr::str_match("^ +([0-9]+)") %>%
                    unlist() %>%
                    magrittr::extract2(2) %>%
                    as.numeric()
                n_reads_aligned_mutiple_times <- lines %>%
                    stringr::str_detect(" aligned[ a-z]* >1 times") %>%
                    magrittr::extract(lines, .) %>%
                    stringr::str_match("^ +([0-9]+)") %>%
                    unlist() %>%
                    magrittr::extract2(2) %>%
                    as.numeric()
                overall_alingment_percent <- lines %>%
                    stringr::str_detect("overall alignment rate") %>%
                    magrittr::extract(lines, .) %>%
                    stringr::str_match("^([0-9.]+)") %>%
                    unlist() %>%
                    magrittr::extract2(2) %>%
                    as.numeric()
                run_time <- lines %>%
                    stringr::str_detect("# Runtime: ") %>%
                    magrittr::extract(lines, .) %>%
                    stringr::str_extract("[0-9.]+$") %>%
                    as.numeric()
                run_time <- ifelse(length(run_time) == 0, NA, run_time)
                tibble::tibble(
                    run_accession = log_fname %>% stringr::str_replace(".log", ""),
                    n_reads = n_reads,
                    n_reads_aligned_0_times = n_reads_aligned_0_times,
                    n_reads_aligned_1_time = n_reads_aligned_1_time,
                    n_reads_aligned_mutiple_times = n_reads_aligned_mutiple_times,
                    overall_alignment_percent = overall_alingment_percent,
                    run_time = run_time)
            }, error = function(e) {
                cat("  ERROR reading log file '", log_fname, "'\n", sep = "")
                tibble::tibble()
            })
	})
}
