#' Use the DeepFRI webserver for gene function prediction
#'
#'  @description Using the https://beta.deepfri.flatironinstitute.org/
#'     webserver Use Structure-based gene function predictions using
#'     graph convolutional networks (Gligorijevic, et al., 2019)
#'
#'  If you use for more than few predictions please coordinate with
#'  the Bonneau lab to not abuse their resources
#'
#' @param pdb_fname proteins structure .pdb file
#' @param tag used to tracking the prediction on the DeepFRI webserver and for
#'    the output path if not using use the basename without .pdb extension.
#'    (default: NULL)
#' @param output_path folder where the output predictions should be written
#'    if null, then don't save output to file. (default: NULL)
#' @param max_retries (default: 120)
#' @param verbose use verbose output (default: TRUE)
#' @return data.frame with columns
#'     structure_tag: <tag>
#'     prediction_type: [cnn, gcn] for sequence and structure based predictions
#'     go_term: e.g. GO:0044237
#'     go_term_name: e.g. "cellular metabolic process"
#'     go_termscore: prediction score
#'  @usage
#'    On disk: <tag>.pdb
#'    In R:
#'
#'     predictions <- get_deepfri_predictions(
#'        pdb_fname = "<tag>.pdb",
#'        output_path = "intermediate_data")
#'
#' If you use please cite:
#'
#'      Gligorijevic, Vladimir and Renfrew, P. Douglas and Kosciolek,
#'      Tomasz and Leman, Julia Koehler and Cho, Kyunghyun and
#'      Vatanen, Tommi and Berenberg, Daniel and Taylor, Bryn and
#'      Fisk, Ian M. and Xavier, Ramnik J. and Knight, Rob and
#'      Bonneau, Richard Structure-Based Function Prediction using
#'      Graph Convolutional Networks 2019,
#'      https://www.biorxiv.org/content/early/2019/10/04/786236
#'
#'
#' @export
get_deepfri_predictions <- function(
  pdb_fname,
  tag = NULL,
  output_path = NULL,
  max_retries = 120,
  verbose = TRUE) {

  if (is.null(tag)) {
    tag <- pdb_fname |> basename() |> stringr::str_replace("[.]pdb$", "")
  }

  if (verbose) {
    cat(
      "Getting DeepFRI gene-function prediction for pdb '", pdb_fname, "'\n",
      "  tag '", tag, "'\n",
      "  saving results to '", output_path, "/", tag, ".tsv'\n", sep = "")
  }

  if (!file.exists(pdb_fname)) {
    cat("  WARNING: pdb '", pdb_fname, "' does not exist\n", sep = "")
  }

  if (!dir.exists(output_path)) {
      cat("  creating output path '", output_path, "'\n", sep = "")
      dir.create(output_path)
  }

  if (verbose) {
    cat("  Submitting prediction for ", tag, " ...\n", sep = "")
  }
  submit_result <- httr::POST(
    url = "https://beta.api.deepfri.flatironinstitute.org/workspace/6PXMJP/predictions",
    httr::add_headers(
      content_type = "multipart/form-data"),
    body = list(
      file = httr::upload_file(
        path = pdb_fname,
        type = "application/octet-stream"),
      inputType = "structureFile",
      tags = c("tag")))
  if (submit_result$status_code != 200) {
    cat("WARNING status of submit request is '", submit_request$status_code, "' is not 200\n", sep = "")
  }
  prediction_name <- httr::content(submit_result)$predictions[[1]]$name
  if (verbose) {
    cat("  Prediction name for tag '", tag, "': ", prediction_name, "\n", sep = "")
  }

  # wait till the results are ready
  retrieve_state <- "submitted"
  n_retries <- 0
  while(retrieve_state %in% c("submitted", "enqueued")) {
    retrieve_result <- httr::GET(
      url = paste0(
        "https://beta.api.deepfri.flatironinstitute.org/workspace/6PXMJP/predictions/",
        prediction_name))

    if (retrieve_result$status_code != 200) {
      cat("WARNING status of retrieve request is '", retrieve_result$status_code, "' is not 200\n", sep = "")
    }
    prediction <- httr::content(retrieve_result)$prediction

    retrieve_state <- ifelse(is.null(prediction$state), "failed", prediction$state)

    n_retries <- n_retries + 1
    if (retrieve_state == "enqueued") {
      if (n_retries < max_retries) {
        if (verbose) {
          cat("  Retrieve retry #", n_retries, "\n", sep = "")
        }
        Sys.sleep(5)
      } else {
        cat("Retried for > 10 minutes, failing\n")
        retrieve_state <<- "failed"
        stop("Failed")
      }
    } else {
      cat("  Retrieve state: ", retrieve_state, "\n", sep = "")
    }
  }

  if (verbose) {
    cat("  Gathering predictions for tag ", tag, "\n", sep = "")
  }

  function_predictions <- dplyr::bind_rows(
    prediction$data$A$cnn_bp$predictions |>
      purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
      dplyr::mutate(prediction_type = "cnn", aspect = "BP", .before = 1),
    prediction$data$A$cnn_cc$predictions |>
     purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
     dplyr::mutate(prediction_type = "cnn", aspect = "CC", .before = 1),
    prediction$data$A$cnn_ec$predictions |>
      purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
      dplyr::mutate(prediction_type = "cnn", aspect = "EC", .before = 1),
    prediction$data$A$cnn_mf$predictions |>
      purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
      dplyr::mutate(prediction_type = "cnn", aspect = "MF", .before = 1),
    prediction$data$A$gcn_bp$predictions |>
      purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
      dplyr::mutate(prediction_type = "gcn", aspect = "BP", .before = 1),
    prediction$data$A$gcn_cc$predictions |>
      purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
      dplyr::mutate(prediction_type = "gcn", aspect = "CC", .before = 1),
    prediction$data$A$gnc_ec$predictions |>
      purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
      dplyr::mutate(prediction_type = "gcn", aspect = "EC", .before = 1),
    prediction$data$A$gnc_mf$predictions |>
      purrr::map_dfr(~.[c("go_term", "go_term_name", "go_term_score")]) |>
      dplyr::mutate(prediction_type = "gcn", aspect = "MF", .before = 1))

    if (length(function_predictions) > 1) {
        function_predictions <- function_predictions |>
            dplyr::mutate(structure_tag = tag, .before = 1)
        if (!is.null(output_path)) {
            if (verbose) {
                cat("  Saving ", nrow(function_predictions),
                    " predictions for ", "tag ", tag, " to ", output_path,
                    "/", tag, ".tsv\n", sep = "")
            }
            function_predictions |>
                readr::write_tsv(
                    file = paste0(output_path, "/", tag, ".tsv"))
        }
        return(function_predictions)
    } else {
        if (verbose) {
            cat("  No predictions found.\n")
        }
        return(data.frame())
    }
}
