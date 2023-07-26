### noisy-simulation-analysis.R ###################################################################
# purpose: run both of the network building for the noisy signal:

### PREABLE #######################################################################################
# load in libraries:
library(minet);
library(Seurat);
library(igraph);
library(parallel);

# write in paths:
simulation.path <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/simulated-data/'
function.dir <- '/u/home/l/lixinzhe/project-github/WQE/code/function/';
r.save <- '/u/home/l/lixinzhe/project-pasaniuc/gym/WQE/R-save/'
system.date <- Sys.Date();

# load in functions:
source(paste0(function.dir, 'progress-bar.R'));
source(paste0(function.dir, 'technical-variance-estimate.R'));
source(paste0(function.dir, 'network-building.R'));
source(paste0(function.dir, 'power-calculation.R'));

# hyperparameters:
seeds <- seq(1, 100);
threads <- 8;
eps.cutoff <- 0; # greater the less edges pruned
n.bin <- 10;
covariation.mean <- 0;
covariation.variance <- 9;
bin.breaks <- qnorm(
    p = seq(0, 1, length.out = n.bin + 1),
    mean = covariation.mean,
    sd = sqrt(covariation.variance)
    );

# load in data:
simulated.count <- vector('list', length = length(seeds));
covariance.matrix <- vector('list', length = length(seeds));
covary.pair <- vector('list', length = length(seeds));
for (seed in seeds){
    simulated.count[[seed]] <- read.table(
        file = paste0(simulation.path, 'seed-', seed, '-positive-batchless-expression.csv'),
        sep = ',',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );

    covariance.matrix[[seed]] <- read.table(
        file = paste0(simulation.path, 'seed-', seed, '-positive-batchless-covariance.csv'),
        sep = ',',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );
    
    covary.pair[[seed]] <- read.table(
        file = paste0(simulation.path, 'seed-', seed, '-positive-batchless-gene-pairs.csv'),
        sep = ',',
        stringsAsFactors = FALSE
        );

    }

### DEPLOY NSCP-ARACNE ############################################################################
nscp.aracne.result <- mclapply(
    simulated.count,
    FUN = function(count) nscp.aracne(
        expression = count,
        batch.effect = NULL,
        technical.var = TRUE,
        count = TRUE,
        eps.cutoff = 0
        ),
    mc.cores = threads
    );

aracne.result <- mclapply(
    simulated.count,
    FUN = function(count) nscp.aracne(
        expression = count,
        batch.effect = NULL,
        technical.var = FALSE,
        count = TRUE,
        eps.cutoff = 0
        ),
    mc.cores = threads
    );

### ANALYSIS - POWER #############################################################################
# calculate power for nscp:
nscp.power <- power.calculation(
    network = nscp.aracne.result,
    latent = covariance.matrix,
    index = covary.pair,
    bin.breaks = bin.breaks
    );

aracne.power <- power.calculation(
    network = aracne.result,
    latent = covariance.matrix,
    index = covary.pair,
    bin.breaks = bin.breaks
    );

# 