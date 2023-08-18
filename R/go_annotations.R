#' Propagate GO term annotations
#'
#' @description GO terms are organized into a directed acyclic graph where an
#' annotation fo a more specific term implies the annotation of the
#' more general.
#'
#' @param go_annotations data.frame with column go_id and object
#'   identifiers. A warning will be issued if it ia columns term,
#'   ontology, or definition as these are added from the
#'   \pkg{GO.db}.
#'
#' @param verbose (Default: TRUE)
#'
#' @return data.frame with columns of input <go_annotations>
#'   where the annotations are propagated up the hierarchy
#' @export
propagate_go_annotations <- function(
    go_annotations,
    verbose = TRUE) {

    if (
        "term" %in% names(go_annotations) ||
        "ontology" %in% names(go_annotations) ||
        "definition" %in% names(go_annotations)) {
        warning(
            paste0(
                "the input go_annotations data.frame shouldn't have ",
                "columns 'term', 'ontology' or 'definition', as these ",
                "are retrieved from the GO terms themselves"))
    }

    if (verbose) {
        cat("Retriving all GO terms ...\n")
    }

    go_terms <- GO.db::GO_dbconn() |>
        dplyr::tbl("go_term") |>
        dplyr::collect(n = Inf)
    if (verbose) {
        cat("Got ", nrow(go_terms), " go terms from Go.db\n", sep = "")
    }

    if (verbose) {
        cat("Getting a map of each go term to each of its parents\n")
    }

    go_parents <- GO.db::GO_dbconn() |>
        dplyr::tbl("go_cc_parents") |>
        dplyr::collect(n = Inf) |>
        dplyr::left_join(
            go_terms |>
            dplyr::select(`_id`, go_id),
            by = "_id") |>
        dplyr::left_join(
            go_terms |>
            dplyr::select(
                `_parent_id` = `_id`,
                parent_go_id = go_id),
            by = "_parent_id") |>
        dplyr::select(-`_id`, -`_parent_id`) |>
        dplyr::filter(parent_go_id != "all")

    go_ancestors <- go_parents |>
        dplyr::select(-relationship_type) |>
        dplyr::rename(ancestor_go_id = parent_go_id)
    n_ancestors <- nrow(go_ancestors)
    while (TRUE) {
        prev_n_ancestors <- n_ancestors
        go_great_ancestors <- go_ancestors |>
            dplyr::inner_join(
                go_parents |> dplyr::transmute(
                    ancestor_go_id = go_id,
                    great_ancestor_go_id = parent_go_id),
                by = "ancestor_go_id",
                relationship = "many-to-many") |>
            dplyr::transmute(
                go_id,
                ancestor_go_id = great_ancestor_go_id) |>
            dplyr::distinct(go_id, ancestor_go_id)

        go_ancestors <- dplyr::bind_rows(
            go_ancestors,
            go_great_ancestors) |>
            dplyr::distinct(go_id, ancestor_go_id)

        n_ancestors <- nrow(go_ancestors)
        if (prev_n_ancestors == n_ancestors) {
            break
        } else {
            cat(
                "adding ", n_ancestors - prev_n_ancestors,
                " additional annotations\n", sep = "")
        }
    }

    if (verbose) {
        cat(
            "Propagating provided ", nrow(go_annotations),
            " annotations to each ancestor.\n", sep = "")
    }

    go_annotations_propagated <- dplyr::bind_rows(
        go_annotations,
        go_annotations |>
            dplyr::inner_join(go_ancestors, by = "go_id") |>
            dplyr::mutate(go_id = ancestor_go_id) |>
            dplyr::select(-ancestor_go_id))

    if (verbose) {
        cat(
            "Found an additional ",
            nrow(go_annotations_propagated) - nrow(go_annotations),
            " annotations\n", sep = "")
    }

    go_annotations_propagated <- go_annotations_propagated |>
        dplyr::left_join(
            go_terms |>
            dplyr::select(-`_id`),
            by = "go_id") |>
        dplyr::filter(!is.na(ontology)) |>
        dplyr::distinct()
}
