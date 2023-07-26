### batch-positive-simulation.R ###################################################################
# purpose: positive simulation with batch effects (batch have off diagnoal covariances)

### PREAMBLE ######################################################################################
# load in libraries:
library(extraDistr);
library(MASS);
library(Matrix)

# specify paths:
out.dir <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/simulated-data/'
function.dir <- '/u/home/l/lixinzhe/project-github/WQE/code/function/';
r.save <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/R-save/'

# set up progress bar function:
source(paste0(function.dir, 'progress-bar.R'));

### SIMULATION ####################################################################################
# set up simulation parameters:
# global:
n.cell <- 500;
n.gene <- 100;
seeds <- seq(1, 100);

# about batch effect:
n.batch.cell <- c(100, 100, 100, 100, 100);
n.batch <- length(n.batch.cell);
min.batch.mean <- 0;
max.batch.mean <- 10;
batch.sd <- 0.5;
batch.precision <- 3;

# about covariance:
n.covariation.values <- 100;
covariation.mean <- 0;
covariation.variance <- 9;
covariation.values <- rnorm(
    n = n.covariation.values,
    mean = covariation.mean,
    sd = sqrt(covariation.variance)
    );

# about global correlation features:
combination.features <- 10;
n.combine <- 5;
combine.weight.mean <- 0;
combine.weight.variance <- 1;

# about Xij:
min.expected.Xij <- 0;
max.expected.Xij <- 10;
min.variance.Xij <- 0;
max.variance.Xij <- 10;

# about technical variance:
min.tech.var <- 0;
max.tech.var <- 10;
tech.mean <- 0;

# about library size:
min.library.size <- 1;
max.library.size <- 100000;

### SIMULATE ######################################################################################
for (seed in seeds) {

    # generate batch effect mean:
    batch.effect.ita <- runif(n = n.batch, min = min.batch.mean, max = max.batch.mean);

    # generated batch effect variance:
    batch.effect.gamma <- abs(batch.effect.ita + rnorm(n = n.batch, mean = 0, sd = batch.sd));

    # generate technical variance:
    tech.var <- runif(n = n.gene, min = min.tech.var, max = max.tech.var);

    # generate library size for each cell:
    library.size <- sample(
        x = min.library.size:max.library.size,
        size = n.cell
        );

    # generate Xij with multivariate normal:
    batch.xij <- vector('list', length = n.batch);
    covariance.matrices <- vector('list', length = n.batch);
    covariance.indexes <- vector('list', length = n.batch);
    linear.combination <- NULL;
    for (batch in seq(1, n.batch)) {
        # generate genes mean and variance within the batch:
        mu.batch <- rnorm(n = n.gene, mean = batch.effect.ita[batch], sd = batch.precision);
        sigma.batch <- rtnorm(
            n = n.gene,
            mean = batch.effect.gamma[batch],
            sd = batch.precision,
            a = 0,
            b = Inf
            );

        # generate covariance matrix:
        covariance.matrix <- matrix(0, ncol = n.gene, nrow = n.gene);
        diag(covariance.matrix) <- sigma.batch;
        rownames(covariance.matrix) <- seq(1, nrow(covariance.matrix));
        colnames(covariance.matrix) <- seq(1, ncol(covariance.matrix));

        # generate off diagonal of covariance matrix for each batch:
        # simulate random gene-gene covariations:
        covary.gene.pair <- replicate(
            n = n.covariation.values,
            expr = sample(
                x = 1:n.gene,
                size = 2,
                replace = FALSE
                )
            );
        # for each covariation values, we will fill it into the covariance matrix:
        for (pair in seq(1, ncol(covary.gene.pair))) {
            index <- covary.gene.pair[, pair];
            covariance.matrix[index[1], index[2]] <- covariation.values[pair];
            covariance.matrix[index[2], index[1]] <- covariation.values[pair];
            }
        covariance.matrix <- as.matrix(nearPD(covariance.matrix, keepDiag = TRUE, maxit = 1e6)$mat);

        # save the covariance matrix:
        rownames(covariance.matrix) <- colnames(covariance.matrix) <- paste0('gene', seq(1, n.gene));
        colnames(covary.gene.pair) <- paste0('pair', seq(1, n.covariation.values));
        rownames(covary.gene.pair) <- paste0('gene', seq(1, 2));
        covariance.matrices[[batch]] <- covariance.matrix;
        covariance.indexes[[batch]] <- covary.gene.pair;
        names(covariance.matrices) <- names(covariance.indexes) <- paste0('batch', seq(1, n.batch));
        
        # generate Xij:
        batch.xij[[batch]] <- mvrnorm(
            n = n.batch.cell[[batch]],
            mu = mu.batch,
            Sigma = covariance.matrices[[batch]]
            );
        }

    # piece the expression matrix together:
    xij <- Reduce('rbind', batch.xij);
    stopifnot(nrow(xij) == n.cell);
    stopifnot(ncol(xij) == n.gene);

    # Add in additional features that are linear combination of other features:
    for (feature in seq(1, combination.features)){
        feature.to.combine <- sample(1:n.gene, size = n.combine);
        weight.to.combine <- rnorm(
            n = n.combine,
            mean = combine.weight.mean,
            sd = sqrt(combine.weight.variance)
            );
        combined.features <- sweep(xij[, feature.to.combine], 2, weight.to.combine, '*');
        combined.features <- rowSums(combined.features);
        xij <- cbind(xij, combined.features);

        # save the weights and the genes that are combined:
        names(weight.to.combine) <- paste0('gene', feature.to.combine);
        linear.combination[[feature]] <- weight.to.combine;
        }
    names(linear.combination) <- paste0('combined', seq(1, combination.features));
    baseline.expression <- xij;

    # Create count:
    count.expression <- baseline.expression;
    non.zero.baseline <- baseline.expression + abs(min(baseline.expression));
    probability.gene <- t(
        apply(
            non.zero.baseline,
            1,
            FUN = function(cell) cell / sum(cell)
            )
        );

    for (cell in seq(1, n.cell)) {
        count.expression[cell, ] <- rmvhyper(
            nn = 1,
            n = round(max.library.size * 1e11 * probability.gene[cell, ]),
            k = library.size[cell]
            );
        }

    # save the count matrix based on seeds:
    colnames(count.expression) <- c(
        paste0('gene', seq(1, n.gene)),
        paste0('combined', seq(1, combination.features))
        );
    rownames(count.expression) <- paste0('cell', seq(1, n.cell));

    # output the simulated data:
    write.table(
        count.expression,
        file = paste0(out.dir, 'seed-', seed, '-batch-positive-expression.csv'),
        sep = ',',
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
        );
    
    # write out the covariance matrix:
    saveRDS(
        covariance.matrices,
        file = paste0(out.dir, 'seed-', seed, '-batch-positive-covariance-matrix.rds')
        );
    saveRDS(
        covariance.indexes,
        file = paste0(out.dir, 'seed-', seed, '-batch-positive-covariance-index.rds')
        );

    # write out the linear combination:
    saveRDS(
        linear.combination,
        paste0(out.dir, 'seed-', seed, '-batch-positive-linear-combination.rds')
        );

    # set up and print progress bar:
    progress_bar(current.step = seed, total.step = length(seeds));
    }
