
library(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
library(monocle3, quietly = TRUE, warn.conflicts = FALSE)
library(CalCEN)


runs <- readr::read_tsv("product/runs_20210803.tsv")
load("intermediate_data/estimated_expression.Rdata")
#load("coexp_intra_study_coexp_spearman.Rdata")

estimated_expression <- estimated_expression %>%
    dplyr::mutate(FPKM_normed = log(FPKM + 1))

if (! dir.exists(paths = "product/figures/coexp_embedding")) {
    cat("Creating 'product/figures/coexp_embedding' ...\n", sep = "")
    dir.create(path = "product/figures/coexp_embedding", recursive = TRUE)
}



#####################
# Embed experiments #
#####################

# translating to monocle terminology
# "gene" <- feature_name a.k.a gene id
# "cell" <- run_accession

exprs <- reshape2::acast(
    data = estimated_expression,
    formula = cnag_id ~ run_accession,
    value.var = "FPKM_normed")

gene_metadata <- data.frame(gene_short_name = exprs %>% colnames())

cell_metadata <-  data.frame(run_accession = exprs %>% colnames) %>%
    dplyr::left_join(
        runs,
        by = "run_accession") %>%
    dplyr::mutate(
        study_accession = as.factor(study_accession))

row.names(gene_metadata) <- gene_metadata$gene_short_name
names(exprs) <- gene_metadata$gene_short_name
row.names(cell_metadata) <- cell_metadata$run_accession

coexp_cds <- monocle3::new_cell_data_set(
    expression_data = exprs,
    gene_metadata = gene_metadata,
    cell_metadata = cell_metadata)

coexp_cds <- coexp_cds %>% monocle3::preprocess_cds(
    num_dims = 100)

coexp_cds <- coexp_cds %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        umap.min_dis = .5,
        umap.n_neighbors = 30L,
        verbose = TRUE)

coexp_cds <- coexp_cds %>%
    monocle3::cluster_cells(
        k = 30,
        num_iter = 10,
        resolution = .1,
        verbose = TRUE)

# id like to label the graph by study accession...
clusters <- cell_metadata$study_accession
names(clusters) <- cell_metadata$run_accession
coexp_cds@clusters[["UMAP"]]$clusters <- clusters

coexp_cds %>%
    monocle3::plot_cells(
        show_trajectory_graph = FALSE,
        cell_size = 0.8)

ggplot2::ggsave(
    filename = paste0(
        "product/figures/coexp_embedding/UMAP_experiments_", CalCEN::date_code(), ".pdf"),
    width = 6,
    height = 6)


clusters <- coexp_cds@clusters[["UMAP"]]$clusters %>%
    data.frame(
        run_accession = exprs %>% colnames(),
        cluster_label = .)

runs <- runs %>%
    dplyr::left_join(
        clusters,
        by = c("run_accession"))

runs %>%
    dplyr::count(study_accession, cluster_label) %>%
    tidyr::pivot_wider(
        id_cols = "study_accession",
        names_from = "cluster_label",
        values_from = "n")


###############
# embed genes #
###############
# translating to monocle terminology
# "gene" expression_data columns <- rna-seq run accession
# "cell" expression_data rows    <- gene (chromosome feature)

exprs <- reshape2::acast(
    data = estimated_expression,
    formula = cnag_id ~ run_accession,
    value.var = "FPKM") %>%
    t()

gene_metadata <- data.frame(run_accession = exprs %>% rownames) %>%
    dplyr::left_join(
        runs %>%
        dplyr::mutate(
            gene_short_name = run_accession),
        by = "run_accession") %>%
    dplyr::mutate(
        study_accession = as.factor(study_accession))
row.names(gene_metadata) <- gene_metadata$gene_short_name

cell_metadata <- data.frame(feature_name = exprs %>% colnames())
row.names(exprs) <- gene_metadata$gene_short_name
row.names(cell_metadata) <- cell_metadata$feature_name

coexp_cds <- monocle3::new_cell_data_set(
    expression_data = exprs,
    gene_metadata = gene_metadata,
    cell_metadata = cell_metadata)

coexp_cds <- coexp_cds %>% monocle3::preprocess_cds(
    num_dims = 500)

coexp_cds <- coexp_cds %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        #umap.min_dis = 0.01,
        #spread = .3,
        a = 50,
        b = 0.5,
        umap.n_neighbors = 30L,
        verbose = TRUE,
        n_epochs = 2000,
        negative_sample_rate = 50,
        repulsion_strength = 3)


coexp_cds <- coexp_cds %>%
    monocle3::reduce_dimension(
        preprocess_method = "PCA",
        verbose = TRUE,
        a = 30,
        b = .8)




coexp_cds <- coexp_cds %>%
    monocle3::cluster_cells(
        k = 30,
        num_iter = 15,
        resolution = .001,
        verbose = TRUE)

save(coexp_cds, file = paste0("intermediate_data/coexp_cds_by_gene_", CalCEN::date_code(), ".Rdata"))


monocle3::plot_cells(cds = coexp_cds,
    show_trajectory_graph = FALSE,
    cell_size = 0.5) +
    ggplot2::coord_fixed() +
    ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
    filename = paste0("product/figures/coexp_embedding/UMAP_genes_", CalCEN::date_code(), ".pdf"),
    width = 15,
    height = 6)

# write out gene embedding and clusters
gene_clusters <- coexp_cds@clusters[["UMAP"]]$clusters %>%
    data.frame(
        cnag_id = exprs %>% colnames(),
        cluster_label = .)

load("intermediate_data/h99_transcript_annotations.Rdata")
gene_clusters <- gene_clusters %>%
    dplyr::left_join(
        h99_transcript_annotations,
        by = "cnag_id")

umap_coordinates <- coexp_cds %>% reducedDim("UMAP") %>%
    as.data.frame() %>%
    dplyr::rename(
        UMAP_1 = 1,
        UMAP_2 = 2)

gene_clusters <- gene_clusters %>%
    dplyr::bind_cols(umap_coordinates)

gene_clusters %>% readr::write_tsv(
    paste0("product/figures/coexp_embedding/UMAP_genes_cluster_labels_", CalCEN::date_code(), ".tsv"))




### plot with old clusters ###
old_gene_clusters <- readr::read_tsv(
    file = paste0(
        "product/figures/coexp_embedding/UMAP_genes_cluster_labels_", CalCEN::date_code(), ".tsv"))

data.frame(
    old_clustering = as.character(old_gene_clusters$cluster_label),
    new_clustering = as.character(gene_clusters$cluster_label)) %>%
    dplyr::count(old_clustering, new_clustering) %>%
    tidyr::pivot_wider(
        names_from = new_clustering,
        values_from = n)



coexp_cds@clusters[["UMAP"]]$clusters <- old_gene_clusters$cluster_label %>% factor()

coexp_cds %>%
    monocle3::plot_cells(
        show_trajectory_graph = FALSE,
        cell_size = 0.8)

ggplot2::ggsave(
    filename = paste0("product/figures/coexp_embedding/UMAP_genes_old_clustering", CalCEN::date_code(), ".pdf"),
    width = 6,
    height = 6)



##### Consensus Clustering #####

n_clusterings <- 2

clusterings <- data.frame(clustering = 1:n_clusterings) %>%
    dplyr::rowwise() %>%
    dplyr::do({
        data <- .
        cat("Computing embedding ", data$clustering, "\n", sep = "")
        coexp_cds <- coexp_cds %>%
            monocle3::reduce_dimension(
                preprocess_method = "PCA",
                init = "random",
                #umap.min_dis = 0.01,
                #spread = .3,
                a = 50,
                b = 0.5,
                umap.n_neighbors = 30L,
                verbose = TRUE,
                n_epochs = 200, #2000,
                negative_sample_rate = 50,
                repulsion_strength = 3)
        coexp_cds <- coexp_cds %>%
            monocle3::cluster_cells(
                k = 30,
                num_iter = 10,
                resolution = .001,
                verbose = TRUE)
        gene_clusters <- coexp_cds@clusters[["UMAP"]]$clusters %>%
            data.frame(
                clustering = data$clustering,
                feature_name = exprs %>% colnames(),
                cluster_label = .)
        gene_clusters
    })

clustering <- clusterings %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_label = as.character(cluster_label)) %>%
    tidyr::pivot_wider(
        id_cols = feature_name,
        names_from = clustering,
        values_from = cluster_label)
# these are all the same


library(Dune)
merger <- Dune::Dune(
    clusMat = clusterings,
    verbose = TRUE)
