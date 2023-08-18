


#' Plot network heatmap
#'
#' Create a grayscale heatmap for large matrix
#'
#' @param network matrix
#' @param output_fname (should be .png)
#' @param gene_order sort the genes using the traveling salesman path using
#'   the seriation package
#' @param verbose (Default: FALSE)
#'
#' @export
plot_network_heatmap <- function(
    network,
    output_fname,
    gene_order = "TSP",
    verbose = FALSE) {

    if (verbose) {
        cat(
            "Plotting network of size",
            " (", dim(network)[1], ", ", dim(network)[2], ") ",
            "to '", output_fname, "'\n", sep = "")
    }
    
    if (!dir.exists(dirname(output_fname))) {
        cat(
            "Output path '", dirname(output_fname), "' does not exist, ",
            "creating ...\n", sep = "")
        dir.create(dirname(output_fname))
    }

    if (verbose) {
        cat("Generating ", gene_order, " gene order...\n", sep = "")
    }
    gene_order <- seriation::seriate(
        x = dist(1 - network),
        method = gene_order) |>
        seriation::get_order()

    if (verbose) {
        cat("Writing out plot...\n")
    }
    
    heatmap <- as.raster(network[gene_order, gene_order])
    png(
        output_fname,
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
}
