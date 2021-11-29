

library(tidyverse)
library(CalCEN)

# download alphafold predictions
system("
cd raw_data/
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/UP000000559_237561_CANAL.tar
untar UP000000559_237561_CANAL.tar
rm -rf *cif.gz
ls *pdb.gz | xargs gunzip -@
")



pdb_fnames <- list.files("raw_data/UP000000559_237561_CANAL/", pattern = "*pdb", full.names=TRUE)

pdb_fnames <- data.frame(pdb_fname = pdb_fnames) %>%
    dplyr::mutate(
        file_exists = pdb_fname %>%
            basename() %>%
            stringr::str_replace("[.]pdb$", "") %>%
            paste0("intermediate_data/af_deep_fri/", ., ".tsv") %>%
            file.exists())

pdb_fnames <- pdb_fnames %>%
    dplyr::mutate(
        uniprot_accession = pdb_fname %>%
            stringr::str_extract("AF-[^-]+") %>%
            stringr::str_replace("AF-", ""))

pdb_fnames %>%
    dplyr::rowwise() %>%
    dplyr::do({
        pdb_fname <- .$pdb_fname[1]
        tag <- basename(pdb_fname) %>% stringr::str_replace("[.]pdb$", "")
        output_fname <- paste0("intermediate_data/af_deep_fri/", tag, ".tsv")
        if (file.exists(output_fname)) {
            cat("Prediction output exsits for '", tag, "', skipping...\n", sep = "")
            predictions <- data.frame()
        } else {
            tryCatch({
                predictions <- CalCEN::get_deepfri_predictions(
                    pdb_fname = pdb_fname,
                    tag = tag,
                    max_retries = 2000,
                    output_path = "intermediate_data/af_deep_fri")
            }, error = function(e){
                cat("Failed to get predictions for '", tag, "' saving to intermediate_data/af_deep_fri/failed.tsv ...")            
                data.frame(pdb_fname = pdb_fname) %>%
                    readr::write_tsv(
                        "intermediate_data/af_deep_fri/failed.tsv",
                        append = TRUE)
                predictions <- data.frame()
            })
        }
        predictions
    })


######################
# gather predictions #
######################

# many uniprot_accession to one primary_cgd_id
uniprot_to_cgd <- readr::read_tsv("raw_data/uniprot-database_(type_CGD).tab") %>%
    dplyr::transmute(
        uniprot_accession = Entry,
        primary_cgd_id = `Cross-reference (CGD)` %>% stringr::str_replace(";.*", ""))

load("intermediate_data/chromosome_features.Rdata")
load("raw_data/ca_genes.Rdata")
deep_fri_predictions <- data.frame(feature_name = ca_genes) %>%
    dplyr::mutate(
        index = dplyr::row_number()) %>%
    dplyr::left_join(
        chromosome_features %>%
        dplyr::select(feature_name, primary_cgd_id),
        by = c("feature_name")) %>%
    dplyr::left_join(
        uniprot_to_cgd,
        by = "primary_cgd_id") %>%
    dplyr::left_join(
        pdb_fnames,
        by = "uniprot_accession")



deep_fri_predictions <- deep_fri_predictions %>%
    dplyr::rowwise() %>%
    dplyr::do({
        metadata <- .
        if (!is.na(metadata$file_exists[1])) {
            pdb_fname <- metadata$pdb_fname[1]
            tag <- basename(pdb_fname) %>% stringr::str_replace("[.]pdb$", "")
            predictions <- paste0("intermediate_data/af_deep_fri/", tag, ".tsv") %>%
                readr::read_tsv(
                    file = .,
                    col_types = readr::cols(
                        structure_tag = readr::col_character(),
                        prediction_type = readr::col_character(),
                        aspect = readr::col_character(),
                        go_term = readr::col_character(),
                        go_term_name = readr::col_character(),
                        go_term_score = readr::col_double())) %>%
                dplyr::mutate(feature_name = metadata$feature_name[1])
        } else {
            data.frame()
        }
    })
                
deep_fri_predictions %>%
    readr::write_tsv("product/deep_fri_predictions_20210904.tsv")

load("intermediate_data/ca_go_annotations.Rdata")
annotations <- ca_go_annotations %>%
    EGAD::filter_network_cols(min = 20, max = 1000)

predictions <- deep_fri_predictions %>%
    dplyr::filter(prediction_type == "cnn") %>%
    dplyr::select(feature_name, go_term, go_term_score) %>%
    tidyr::pivot_wider(
        id_cols = feature_name,
        names_from = go_term,
        values_from = go_term_score,
        values_fill = 0) %>%
    tibble::column_to_rownames(var = "feature_name") %>%
    EGAD::filter_network_cols(min = 20, max = 1000)

#evalaute_annotation_predictions <- function(
#    annotations,
#    predictions){

    common_genes <- intersect(rownames(annotations), rownames(predictions))
    common_terms <- intersect(colnames(annotations), colnames(predictions))
    annotations <- annotations[
        match(common_genes, rownames(annotations)),
        match(common_terms, colnames(annotations))]
    predictions <- predictions[
        match(common_genes, rownames(predictions)),
        match(common_terms, colnames(predictions))]
    assertthat::assert_that(
        all(rownames(annotations) == rownames(predictions)),
        msg = "matching gene names")
    assertthat::assert_that(
        all(colnames(annotations) == colnames(predictions)),
        msg = "matching term names")

    n_genes <- length(common_genes)
    n_terms <- length(common_terms) 
    np <- colSums(annotations)
    nn <- dim(annotations)[1] - colSums(annotations)
    predictions <- apply(abs(predictions), 2, rank, na.last = "keep", 
        ties.method = "average")
    predictions[negatives] <- 0
    p <- apply(predictions, 2, sum, na.rm = TRUE)
    rocN <- (p/np - (np + 1)/2)/nn
}

auc_scores <- data.frame(
    go_term = names(rocN),
    auroc = rocN) %>%
    dplyr::left_join(
        deep_fri_predictions %>%
        dplyr::distinct(go_term, go_term_name),
        by = "go_term")

