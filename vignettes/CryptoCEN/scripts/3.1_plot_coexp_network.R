
library(plyr)
library(dplyr)
library(stringr)
library(readr)
library(seriation)
library(EGAD)
library(CalCEN)

load("intermediate_data/estimated_expression.Rdata")
load("product/coexp_intra_study_coexp_network_spearman_20210812.Rdata")


#############################
# Plot co-expression matrix #
#############################

irlba::irlba(
    A = network,
    nv = min(pca_dim, min(dim(network) - 1)),
    center = network |> colMeans2(),
    scale = network |> colVars() |> sqrt())

irlba_res <- sparse_prcomp_irlba(
    Matrix::t(FM),
    n = min(num_dim, min(dim(FM)) - 1),
    center = scaling,
    scale. = scaling)
preproc_res <- irlba_res$x
row.names(preproc_res) <- colnames(cds)
irlba_rotation <- irlba_res$rotation
row.names(irlba_rotation) <- rownames(FM)
cds@preprocess_aux$gene_loadings <- irlba_rotation
cds@preprocess_aux$prop_var_expl <- irlba_res$sdev^2/sum(irlba_res$sdev ^ 2)


coexp_intra_study_coexp_network_spearman |>
    CalCEN::plot_network_heatmap(
        output_fname = paste0(
            "product/figures/coexp_network_heatmap_TSP_raster_", CalCEN::date_code(), ".png"),
        gene_order = "TSP",
        verbose = TRUE)

coexp_intra_study_coexp_network_spearman |>
    CalCEN::plot_network_heatmap(
        output_fname = paste0(
            "product/figures/coexp_network_heatmap_OLO_raster_", CalCEN::date_code(), ".png"),
        gene_order = "OLO",
        verbose = TRUE)


### OLO ##
gene_order <- seriation::seriate(
    x = dist(1-CalCEN_network_full),
    method = "OLO") |>
    seriation::get_order()

heatmap <- as.raster(CalCEN_network_full[gene_order, gene_order])
png(
    "product/figures/CalCEN_heatmap_raster1_20201130.png",
    width = length(gene_order),
    height = length(gene_order))
plot.new()
rasterImage(
    image = heatmap,
    xleft = 0,
    ybottom = 0,
    xright = 1,
    ytop = 1,
    interpolate = FALSE)
dev.off()

### TSP ##
gene_order <- seriation::seriate(
    x = dist(1 - CalCEN_network_full),
    method = "TSP") |>
    seriation::get_order()

heatmap <- as.raster(CalCEN_network_full[gene_order, gene_order])
png(
    "product/figures/CalCEN_heatmap_raster_TSP_20201112.png",
    width = length(gene_order),
    height = length(gene_order))
plot.new()
rasterImage(
    image = heatmap,
    xleft = 0,
    ybottom = 0,
    xright = 1,
    ytop = 1,
    interpolate = FALSE)
dev.off()

###########################
# Plot expression heatmap #
###########################


estimated_expression <- estimated_expression |>
    dplyr::semi_join(ca_runs_final, by = c("run_accession"))

exprs <- reshape2::acast(
  data = estimated_expression,
  formula = gene_id ~ run_accession,
  value.var = "FPKM")

exprs <- log(1 + exprs)

exprs_ranked <- exprs |>
    rank(na.last = "keep", ties.method = "average") |>
    matrix(nrow = 6226, ncol = 853)
exprs_ranked <- exprs_ranked / max(exprs_ranked)

exprs_gene_dist <- dist(1 - exprs_ranked)
exprs_run_dist <- dist(1 - t(exprs_ranked))

gene_order <- seriation::seriate(
    x = exprs_gene_dist,
    method = "OLO") |>
    seriation::get_order()

run_order <- seriation::seriate(
    x = exprs_run_dist,
    method = "OLO") |>
    seriation::get_order()

heatmap <- as.raster(exprs_ranked[gene_order, run_order])

png(
    "product/figures/ca_expression_heatmap_raster_20201203.png",
    width = length(run_order),
    height = length(gene_order))
plot.new()
rasterImage(
    image = heatmap,
    xleft = 0,
    ybottom = 0,
    xright = 1,
    ytop = 1,
    interpolate = FALSE)
dev.off()


# label expression matrix with "Rank Expression" with values [0-1]
#   (black to white)
# label co-exprssion matrix with "Rank Co-expression" with values [0-1]
#   (black to white)
scale_bar <- seq(0, 1, length.out = 11) |>
    rep(3) |>
    matrix(nrow = 3, ncol = 11, byrow = TRUE)
png(
    "product/figures/heatmap_raster_scale_bar_20201203.png",
    width = 11 * 50,
    height = 3 * 50)
plot.new()
rasterImage(
    image = scale_bar,
    xleft = 0,
    ybottom = 0,
    xright = 1,
    ytop = 1,
    interpolate = FALSE)
dev.off()
