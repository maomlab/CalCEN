
#' Propagate GO term annotations
#'
#' GO terms are organized into a directed acyclic graph where an
#' annotation fo a more specific term implies the annotation of the
#' more general.
#'
#' GO terms are not 
#' 
#' @param go_annotations data.frame with column go_id and object identifiers
#' @param verbose (Default: TRUE)
#'   
#' @return data.frame with columns of input <go_annotations>
#'   where the annotations are propagated up the hierarchy
#' @export
propagate_go_annotations <- function(
    go_annotations,
    verbose = TRUE) {

    if (verbose) {
        cat("Retriving all GO terms ...\n")
    }
    
    go_terms <- GO.db::GO_dbconn() %>%
        dplyr::tbl("go_term") %>%
        dplyr::collect(n = Inf)

    if (verbose) {
        cat("Getting a map of each go term to each of its parents\n")
    }
    go_parents <- GO.db::GO_dbconn() %>%
        dplyr::tbl("go_cc_parents") %>%
        dplyr::collect(n = Inf) %>%
        dplyr::left_join(
            go_terms %>%
            dplyr::select(`_id`, go_id),
            by = "_id") %>%
        dplyr::left_join(
            go_terms %>%
            dplyr::select(
                `_parent_id` = `_id`,
                parent_go_id = go_id),
            by = "_parent_id") %>%
        dplyr::select(-`_id`, -`_parent_id`) %>%
        dplyr::filter(parent_go_id != "all")
        

    if (verbose) {
        cat("Propagating provided ", nrow(go_annotations), " annotations to each parent.\n", sep = "")
    }

    go_annotations_propagated <- dplyr::bind_rows(
        go_annotations,
        go_annotations %>%
            dplyr::inner_join(go_parents, by = "go_id") %>%
	    dplyr::mutate(go_id = parent_go_id) %>%
            dplyr::select(-parent_go_id, -relationship_type))
    
    if (verbose) {
        cat(
            "Found an additional ",
            nrow(go_annotations_propagated) - nrow(go_annotations),
            " annotations\n", sep = "")
    }
    
    go_annotations_propagated <- go_annotations_propagated %>%
	dplyr::left_join(
            go_terms %>%
            dplyr::select(-`_id`),
            by = "go_id") %>%
	dplyr::filter(!is.na(ontology)) %>%
        dplyr::distinct()
}
