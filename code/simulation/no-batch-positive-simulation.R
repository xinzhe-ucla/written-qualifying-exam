### no-batch-positive-simulation.R ################################################################
# purpose: simulation of positive signal when there is no batches:

### PREAMBLE ######################################################################################
# load in library:
library(extraDistr);
library(MASS);
library(Matrix);

# Specify paths:
out.dir <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/simulated-data/'
function.dir <- '/u/home/l/lixinzhe/project-github/WQE/code/function/';
r.save <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/R-save/'

# set up progress bar function:
source(paste0(function.dir, 'progress-bar.R'));

# set up simulation parameters:
n.cell <- 500;
n.batch <- 10;
min.expected.Xij <- 0;
max.expected.Xij <- 10;
variance.deviation <- 1;
min.library.size <- 1;
max.library.size <- 100000;
n.gene <- 100;
seeds <- seq(1, 100);
n.covariation.values <- 100;
covariation.mean <- 0;
covariation.variance <- 9;
covariation.values <- rnorm(
    n = n.covariation.values,
    mean = covariation.mean,
    sd = sqrt(covariation.variance)
    );

### SIMULATION ####################################################################################
for (seed in seeds) {
    # set seed for the range of seeds:
    set.seed(seed);

    # first generate baseline expression:
    expression.mean <- runif(
        n = n.gene,
        min = min.expected.Xij,
        max = max.expected.Xij
        );

    expression.variance <- rtnorm(
        n = n.gene,
        mean = expression.mean,
        sd = variance.deviation,
        a = 0,
        b = Inf
        );

    library.size <- sample(
        x = min.library.size:max.library.size,
        size = n.cell
        );

    # next simulate baseline expression level:
    # use a covariance matrix to simulate:
    covariance.matrix <- matrix(0, nrow = n.gene, ncol = n.gene);
    diag(covariance.matrix) <- expression.variance;
    rownames(covariance.matrix) <- seq(1, nrow(covariance.matrix));
    colnames(covariance.matrix) <- seq(1, ncol(covariance.matrix));

    # simulate random gene-gene covariations:
    covary.gene.pair <- replicate(
        n = n.covariation.values,
        expr = sample(
            x = 1:n.gene,
            size = 2,
            replace = FALSE)
            );
    # for each covariation values, we will fill it into the covariance matrix:
    for (pair in seq(1, ncol(covary.gene.pair))) {
        index <- covary.gene.pair[, pair];
        covariance.matrix[index[1], index[2]] <- covariation.values[pair];
        covariance.matrix[index[2], index[1]] <- covariation.values[pair];
        }
    
    covariance.matrix <- as.matrix(nearPD(covariance.matrix, keepDiag = TRUE, maxit = 1e6)$mat);
    baseline.expression <- mvrnorm(
        n = n.cell,
        mu = expression.mean,
        Sigma = covariance.matrix
        );

    colnames(baseline.expression) <- paste0('gene', seq(1, n.gene));
    rownames(baseline.expression) <- paste0('cell', seq(1, n.cell));

    # directly simulate count using multinomial distribution without replacement:
    count.expression <- baseline.expression;
    non.zero.baseline <- baseline.expression + abs(min(baseline.expression));
    probability.gene <- t(
        apply(
            non.zero.baseline,
            1,
            FUN = function(cell) cell / sum(cell)
            )
        );
    # check genes for each cells can sum to 1:
    stopifnot(all(signif(rowSums(probability.gene), 10) == 1))

    # simulate count from the result:
    for (cell in seq(1, n.cell)) {
        count.expression[cell, ] <- rmvhyper(
            nn = 1,
            n = round(max.library.size * 1e11 * probability.gene[cell, ]),
            k = library.size[cell]
            );
        }

    # output the simulated null data:
    write.table(
        count.expression,
        file = paste0(out.dir, 'seed-', seed, '-positive-batchless-expression.csv'),
        sep = ',',
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
        );

    # also output the covariance matrix:
    write.table(
        covariance.matrix,
        file = paste0(out.dir, 'seed-', seed, '-positive-batchless-covariance.csv'),
        sep = ',',
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
        );

    # write out the covarying gene pairs as well:
    write.table(
        covary.gene.pair,
        file = paste0(out.dir, 'seed-', seed, '-positive-batchless-gene-pairs.csv'),
        sep = ',',
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
        );

    # set up and print progress bar:
    progress_bar(current.step = seed, total.step = length(seeds));
    }

#