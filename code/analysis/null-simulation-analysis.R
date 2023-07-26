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

# hyperparameters:
seeds <- seq(1, 100);
threads <- 8;
percentile.cutoff <- 0.95; # greater the more edges pruned
eps.cutoff <- 0; # greater the less edges pruned

# load in data:
simulated.count <- vector('list', length = length(seeds));
for (seed in seeds){
    simulated.count[[seed]] <- read.table(
        file = paste0(simulation.path, 'seed-', seed, 'sim-null-expression.csv'),
        sep = ',',
        header = TRUE,
        row.names = 1,
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

# compute result:
aracne.shifted <- aracne.original <- vector('list', length = length(aracne.result));
shifted.zero.counts <- original.zero.counts <- NULL;
for (mim in seq(1, length(aracne.result))) {
    shifted.zero.count <- sum(nscp.aracne.result[[mim]] == 0);
    original.zero.count <- sum(aracne.result[[mim]] == 0);
    shifted.zero.counts <- c(shifted.zero.counts, shifted.zero.count);
    original.zero.counts <- c(original.zero.counts, original.zero.count);
    }
cat('mean of shifted zeros:', mean(shifted.zero.count), '\n')
cat('mean of original zeros:', mean(original.zero.count), '\n')

# save result:
saveRDS(aracne.result, file = paste0(r.save, system.date, '-aracne-result-under-null.rds'));
saveRDS(nscp.aracne.result, file = paste0(r.save, system.date, '-nscp-aracne-result-under-null.rds'));
