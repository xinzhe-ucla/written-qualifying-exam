### noisy-simulation.R ###########################################################################
# purpose: simulate expression data that is noisy with batch effects

### PREAMBLE ######################################################################################
# load in library:
library(extraDistr);
library(MASS);

# specify paths:
out.dir <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/simulated-data/'
function.dir <- '/u/home/l/lixinzhe/project-github/WQE/code/function/';
r.save <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/R-save/'

# set up progress bar function:
source(paste0(function.dir, 'progress-bar.R'));

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
        cov.matrix <- matrix(0, ncol = n.gene, nrow = n.gene);
        diag(cov.matrix) <- sigma.batch;

        # generate Xij:
        batch.xij[[batch]] <- mvrnorm(
            n = n.batch.cell[[batch]],
            mu = mu.batch,
            Sigma = cov.matrix
            );
        }

    # piece the expression matrix together:
    xij <- Reduce('rbind', batch.xij);
    stopifnot(nrow(xij) == n.cell);
    stopifnot(ncol(xij) == n.gene);

    # Add in the technical noise:
    baseline.expression <- xij;
    # for (gene in seq(1,ncol(xij))) {
    #     technical.noise <- rnorm(n.cell, mean = tech.mean, sd = sqrt(tech.var[gene]));
    #     baseline.expression[,gene] <- xij[, gene] + technical.noise;
    #     }

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
    
    # output the simulated null data:
    write.table(
        xij,
        file = paste0(out.dir, 'seed-', seed, 'noisy-null-expression.csv'),
        sep = ',',
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
        );

    # set up and print progress bar:
    progress_bar(current.step = seed, total.step = length(seeds));
    }
