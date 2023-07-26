### noisy-simulation-analysis.R ##################################################################
# purpose: run network checking for noisy simulation

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
percentile.cutoff <- 0.95; # greater the more edges pruned
eps.cutoff <- 0; # greater the less edges pruned

# load in data:
simulated.count <- vector('list', length = length(seeds));
for (seed in seeds){
    simulated.count[[seed]] <- read.table(
        file = paste0(simulation.path, 'seed-', seed, 'noisy-null-expression.csv'),
        sep = ',',
        header = TRUE,
        row.names = 1,
        stringsAsFactors = FALSE
        );
    }

# add in the batch effect:
batch.effect.label <- c(
    rep(1, 100),
    rep(2, 100),
    rep(3, 100),
    rep(4, 100),
    rep(5, 100)
    );
cov.matrix <- data.frame(batch.effect.label);
rownames(cov.matrix) <- seq(1, 500);

# Run ARACNE:
nscp.aracne.result <- mclapply(
    simulated.count,
    FUN = function(count) nscp.aracne(
        expression = count,
        batch.effect = cov.matrix,
        technical.var = FALSE,
        count = FALSE,
        eps.cutoff = 0
        ),
    mc.cores = threads
    );

# run aracne results:
aracne.result <- mclapply(
    simulated.count,
    FUN = function(count) nscp.aracne(
        expression = count,
        batch.effect = NULL,
        technical.var = FALSE,
        count = FALSE,
        eps.cutoff = 0
        ),
    mc.cores = threads
    );

# compute result:
shifted.edge.counts <- batch.edge.counts <- original.edge.counts <- NULL;
for (mim in seq(1, length(aracne.result))) {
    shifted.edge.count <- sum(nscp.aracne.result[[mim]] != 0);
    batch.edge.count <- sum(nscp.aracne.result[[mim]][101, ] != 0);
    original.edge.count <- sum(aracne.result[[mim]] != 0);

    shifted.edge.counts <- c(shifted.edge.counts, shifted.edge.count);
    batch.edge.counts <- c(batch.edge.counts, batch.edge.count);
    original.edge.counts <- c(original.edge.counts, original.edge.count);
    }

cat('ARACNE edge counts:', mean(original.edge.counts), '\n')
cat('batch edge counts:', mean(batch.edge.counts), '\n')
cat('NSCP edge counts:', mean(shifted.edge.counts), '\n')

# save result:
saveRDS(aracne.result, file = paste0(r.save, system.date, '-aracne-result-under-null-with-noise.rds'));
saveRDS(nscp.aracne.result, file = paste0(r.save, system.date, '-nscp-aracne-result-under-null-with-noise.rds'));
