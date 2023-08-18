library(plyr)
library(dplyr)
library(magrittr)
library(stringr)
library(readr)
library(reshape2)
library(EGAD)
library(Hotelling)
library(huge)

parameters <- CalCEN::load_parameters()

runs <- readr::read_tsv(
    "product/runs_20210803.tsv",
    show_col_types = FALSE)

load("intermediate_data/h99_transcript_annotations.Rdata")
genes <- h99_transcript_annotations$cnag_id

load("intermediate_data/estimated_expression.Rdata")

exprs <- reshape2::acast(
    data = estimated_expression,
    formula = cnag_id ~ run_accession,
    value.var = "FPKM") |>
    magrittr::extract(genes, )

#############################
# spearman rank correlation #
#############################

coexp_network_full_spearman <- EGAD::build_coexp_network(
  exprs = exprs,
  gene.list = genes)


save(
    coexp_network_full_spearman,
    file = paste0(
        "product/coexp_network_full_spearman_", CalCEN::date_code(), ".Rdata"))

coexp_network_full_spearman |>
  as.data.frame() |>
    readr::write_tsv(
        file = paste0(
            "product/coexp_network_full_spearman_", CalCEN::date_code(), ".tsv"))

coexp_full_spearman <- coexp_network_full_spearman |>
  data.frame() |>
  tibble::rownames_to_column("cnag_id_1") |>
    tidyr::gather(
        key = "cnag_id_2",
        value = "score",
        -cnag_id_1) |>
    dplyr::mutate(
        cnag_id_2 = cnag_id_2 |> stringr::str_replace("[.]", "-"))

save(
    coexp_full_spearman,
    file = "intermediate_data/coexp_full_spearman.Rdata")

#####################
# intra study coexp #
#####################

coexp_intra_study_coexp_network_spearman <- estimated_expression |>
    dplyr::semi_join(runs, by = c("run_accession")) |>
    dplyr::group_by(study_accession) |>
    dplyr::group_map(function(study_expression, study) {
        exprs <- reshape2::acast(
            data = study_expression,
            formula = cnag_id ~ run_accession,
            value.var = "FPKM") |>
            magrittr::extract(genes, )
        cat(
            "building network for study: ",
            study$study_accession[1],
            " with ", ncol(exprs), " samples ...\n", sep = "")
        network <- EGAD::build_coexp_network(
            exprs = exprs,
            gene.list = genes)
        ifelse(is.na(network), .5, network)
    }) |>
    purrr::reduce(`+`)


coexp_intra_study_coexp_network_spearman <-
    coexp_intra_study_coexp_network_spearman /
    length(unique(runs$study_accession))

save(
    coexp_intra_study_coexp_network_spearman,
    file = paste0(
        "product/coexp_intra_study_coexp_network_spearman_",
        CalCEN::date_code(), ".Rdata"))

coexp_intra_study_coexp_network_spearman |>
  as.data.frame() |>
    readr::write_tsv(
        file = paste0(
            "product/coexp_intra_study_coexp_network_spearman_",
            CalCEN::date_code(), ".tsv"))


coexp_intra_study_coexp_spearman <- coexp_intra_study_coexp_network_spearman |>
    data.frame() |>
    tibble::rownames_to_column("cnag_id_1") |>
    tidyr::gather(
        key = "cnag_id_2",
        value = "score",
        -cnag_id_1) |>
    dplyr::mutate(
        cnag_id_2 = cnag_id_2 |> stringr::str_replace("[.]", "-"))

save(
    coexp_intra_study_coexp_spearman,
    file = "intermediate_data/coexp_intra_study_coexp_spearman.Rdata")


###################################
# Pearson Correlation coefficient #
###################################

coexp_network_full_pearson <- EGAD::build_coexp_network(
    exprs = exprs,
    gene.list = genes,
    method = "pearson",
    flag = FALSE)

save(
    coexp_network_full_pearson,
    file = paste0("product/coexp_network_full_pearson_",
        CalCEN::date_code(), ".Rdata"))
coexp_network_full_pearson |>
  as.data.frame() |>
    readr::write_tsv(
        file = paste0(
            "product/coexp_network_full_pearson_", CalCEN::date_code(), ".tsv"))

coexp_full_pearson <- coexp_network_full_pearson |>
    data.frame() |>
    tibble::rownames_to_column("cnag_id_1") |>
    tidyr::gather(
        key = "cnag_id_2",
        value = "score",
        -cnag_id_1) |>
    dplyr::mutate(
        cnag_id_2 = cnag_id_2 |> stringr::str_replace("[.]", "-"))

save(
    coexp_full_pearson, file = "intermediate_data/coexp_full_pearson.Rdata")

######################
# Direct correlation #
######################

# more numerically stable version of Hotelling::clr
# https://stackoverflow.com/q/2602583/198401
clr <- function(data) {
    log_gms <- apply(data, 1, function(x) {mean(log(x))})
    log(data) - log_gms
}

nlambda <- 60
coexp_full_direct <- (exprs + 1) |>
    t() |>
    clr() |>
    huge::huge.npn() |>
    huge::huge(method = "mb", nlambda = nlambda)

save(coexp_full_direct, file = "intermediate_data/coexp_full_direct.Rdata")

coexp_full_direct_bootstrap_60 <- coexp_full_direct |>
     huge::huge.select(criterion = "stars", stars.thresh = 0.05)

for (i in c(1:nlambda)) {
    coexp_full_direct_bootstrap_60$path[[i]]@Dimnames <-
        list(as.character(genes), as.character(genes))
    coexp_full_direct_bootstrap_60$beta[[i]]@Dimnames <-
        list(as.character(genes), as.character(genes))
    coexp_full_direct_bootstrap_60$merge[[i]]@Dimnames <-
        list(as.character(genes), as.character(genes))
}

save(
    coexp_full_direct_bootstrap_60,
    file = "intermediate_data/coexp_full_direct_bootstrap_60.Rdata")
