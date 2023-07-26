### simulation-null-correlation-null-noise.R ######################################################
# purpose: simulate data that doesn't have any correlation structure or added noise
# to show that there is no pattern that the minet detect

### PREAMBLE ######################################################################################
# load in library:
library(extraDistr);
library(MASS);

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

    # nex simulate baseline expression level:
    baseline.expression <- matrix(NA, ncol = n.gene, nrow = n.cell);
    for (gene in seq(1, n.gene)) {
        baseline.expression[, gene] <- rnorm(
            n = n.cell,
            mean = expression.mean[gene],
            sd = sqrt(expression.variance[gene])
            );
        }
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
        file = paste0(out.dir, 'seed-', seed, 'sim-null-expression.csv'),
        sep = ',',
        quote = FALSE,
        row.names = TRUE,
        col.names = TRUE
        );

    # set up and print progress bar:
    progress_bar(current.step = seed, total.step = length(seeds));
}
