### network-building ##############################################################################
library(minet);
library(Seurat);
library(igraph);
library(parallel);

# load in technical variance calculation function:
source('/u/home/l/lixinzhe/project-github/WQE/code/function/technical-variance-estimate.R')

# NSCP-ARACNE
nscp.aracne <- function(
    expression,
    batch.effect = NULL,
    technical.var = TRUE,
    count = TRUE,
    eps.cutoff = 0,
    n.bins = 10
    ){
        # module for handling the expression matrix:
        if (count) {
            # create the count matrix: 
            sim.seurat <- CreateSeuratObject(
                counts = t(expression),
                project = "CreateSeuratObject",
                assay = "RNA",
                names.field = 1,
                names.delim = "_",
                min.cells = 3,
                min.features = 10,
                meta.data = NULL
                );
            
            # normalization:
            sim.seurat <- NormalizeData(
                sim.seurat,
                normalization.method = "LogNormalize",
                scale.factor = 1e6
                );

            # extract out the normalized simulation:
            normalized.sim <- t(as.matrix(sim.seurat@assays$RNA@data));
            }
            
        else {
            normalized.sim <- expression;
            }
        
        # module for estimating technical variance:
        if (technical.var) {
            technical.var <- compute.technical.variance(normalized.sim, span = 0.3, degree = 2);
            }

        else {
            technical.var <- rep(1, ncol(normalized.sim));
            }

        # module for building MIM with and without batch effect:
        if (is.null(batch.effect)) {
            normalized.mim <- build.mim(
                dataset = data.frame(normalized.sim),
                estimator = 'mi.mm',
                disc = 'equalwidth',
                nbins = n.bins
                );

            # shift by technical variance before entering ARACNE:
            shifted.mim <- sapply(
                X = seq(1, length(technical.var)),
                FUN = function(column) normalized.mim[, column] / technical.var[column]
                );

            # rename the rownames and column names:
            rownames(shifted.mim) <- rownames(normalized.mim);
            colnames(shifted.mim) <- colnames(normalized.mim);
            }
        else {
            normalized.sim <- cbind(
                normalized.sim,
                batch.effect[row.names(normalized.sim), , drop = FALSE]
                );
            normalized.mim <- build.mim(
                dataset = data.frame(normalized.sim),
                estimator = 'mi.mm',
                disc = 'equalwidth',
                nbins = n.bins
                );

            # shift by technical variance:
            for (gene in seq(1, length(technical.var))) {
                normalized.mim[, gene] <- normalized.mim[, gene] / technical.var[gene];
                }
            shifted.mim <- normalized.mim;
            }

        # module for building network:
        # compute result:
        aracne.network <- aracne(shifted.mim, eps = eps.cutoff);
        return(aracne.network);
        }