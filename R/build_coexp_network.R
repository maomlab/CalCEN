#' Build Co-exp network for each study and average
#'
#' @param estimated_expression data.frame with columns
#'        study_accession,
#'        run_accession,
#'        gene_id,
#'        FPKM
#' @param genes vector of gene labels
#' @param verbose verbose output (default: FALSE)
#' @return n_gene by n_gene matrix with values with the average
#'     spearman rank correlation across studies of of the FPKM values
#'     for each run in the study
#' @export
build_coexp_intra_study_spearman_network <- function(
     estimated_expression,
     genes,
     verbose = FALSE) {

    if (
        !("study_accesion" %in% names(estimated_expression)) |
        !("run_accession" %in% names(estimated_expression)) |
        !("gene_id" %in% names(estimated_expression)) |
        !("FPKM" %in% names(estimated_expression))) {

        stop(
            "estimated expression must have columns ",
            "['study_accession', 'run_accession', 'gene_id', 'FPKM']\n",
            "instead it has columns ",
            "['", paste(names(estimated_expression), collapse = "','"), "']\n",
            sep = "")
    }

    integrated_network <- matrix(
        data = 0,
        nrow = length(genes),
        ncol = length(genes),
        dimnames = list(genes, genes))

    estimated_expression |>
        dplyr::group_by(study_accession) |>
        dplyr::group_map(function(study_expression, study) {
            exprs <- reshape2::acast(
                data = study_expression,
                formula = gene_id ~ run_accession,
                value.var = "FPKM") |>
                magrittr::extract(genes, )
            if (verbose) {
                cat(
                    "building network for study: ",
                    study$study_accession[1],
                    " with ", ncol(exprs), " samples ...\n", sep = "")
            }
            network <- EGAD::build_coexp_network(
                exprs = exprs,
                gene.list = genes)
            ifelse(is.na(network), .5, network)
            integrated_network <<- integrated_network + network
            NULL
        })
    integrated_network / length(unique(estimated_expression$study_accession))
}

#' Build Direct Co-Expression network
#'
#' @description Use the graphical lasso to estimate a sparse
#' co-expression matrix huge package
#'
#' @param estimated_expression data.frame with columns gene_id
#'     run_accession
#' @param genes vector of gene labels
#' @param verbose verbose output (default: FALSE)
#' @return result of huge::huge.select In particular
#'     .$merge[<nlambda>] is a n_gene by n_gene co-expression matrix
#' @export
build_coexp_direct_network <- function(
     estimated_expression,
     genes,
     nlambda = 30,
     verbose = FALSE) {

    exprs <- reshape2::acast(
        data = estimated_expression,
        formula = gene_id ~ run_accession,
        value.var = "FPKM") |>
        magrittr::extract(genes, )

    # more numerically stable version of Hotelling::clr
    # https://stackoverflow.com/q/2602583/198401
    clr <- function(data) {
        log_gms <- apply(data, 1, function(x){mean(log(x))})
        log(data) - log_gms
    }

    coexp_direct <- (exprs + 1) |>
        t() |>
        clr() |>
        huge::huge.npn() |>
        huge::huge(method = "mb", nlambda = nlambda)

    coexp_direct <- full_direct |>
        huge::huge.select(criterion = "stars", stars.thresh = 0.05)

    for (i in c(1:nlambda)) {
        coexp_direct$path[[i]]@Dimnames <-
            list(as.character(genes), as.character(genes))
        coexp_direct$beta[[i]]@Dimnames <-
            list(as.character(genes), as.character(genes))
        coexp_direct$merge[[i]]@Dimnames <-
            list(as.character(genes), as.character(genes))
    }
    coexp_direct
}
