### batch-positive-simulation-analysis.R ##########################################################
# purpose: analyze batch positive simulation analysis and check proportion of edges:

### PREAMBLE ######################################################################################
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

# hyperparameters:
seeds <- seq(1, 100);
threads <- 1;
eps.cutoff <- 0; # greater the less edges pruned

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
linear.combination <- vector('list', length = length(seeds));

for (seed in seeds){
    simulated.count[[seed]] <- read.table(
        file = paste0(simulation.path, 'seed-', seed, '-batch-positive-expression.csv'),
        sep = ',',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );

    covariance.matrix[[seed]] <- readRDS(
        file = paste0(simulation.path, 'seed-', seed, '-batch-positive-covariance-matrix.csv')
        );
    
    covary.pair[[seed]] <- readRDS(
        file = paste0(simulation.path, 'seed-', seed, '-batch-positive-covariance-index.csv')
        );

    linear.combination[[seed]] <- readRDS(
        file = paste0(simulation.path, 'seed-', seed, '-batch-positive-linear-combination.csv')
        );
    }

### ANALYSIS ######################################################################################
# for each of the simulated data, find out the power across batches:

