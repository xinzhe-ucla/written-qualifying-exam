### power calculation.R ###########################################################################
# purpose: calculate power for each quantile

# INPUT:
# network: list of aracne styled networks 
# latent:  latent covariance matrix used to generate the data
# index: index of covariance matrix (2 * pairs of pertubed edges)
# quantile.bins: number of bins you quantiling covariance and computing power
# bin.breaks: break points for binning the covariance matrix
# return: results that you want, support 'power' and 'all'

# OUTPUT:
power.calculation <- function(network, latent, index, bin.breaks, return = 'power') {
    # create three empty place holder:
    latent.perturbation <- matrix(NA, nrow = length(seeds), ncol = ncol(index[[1]]));
    nscp.edge <- latent.perturbation;
    seeds <- seq(1, length(network));
    n.bin <- length(bin.breaks) - 1;
    nscp.classification <- matrix(NA, nrow = length(seeds), ncol = n.bin);
    rownames(nscp.classification) <- paste0('simulation', seeds);
    colnames(nscp.classification) <- paste0('bin', seq(1, n.bin));
    perturbed.pairs.bin <- nscp.classification;

    for (simulation in seeds) {
        # obtain the gene pairs that has been altered:
        for (pair in seq(1, ncol(index[[simulation]]))) {
            # grab out the first and the second gene:
            gene1 <- index[[simulation]][, pair][1];
            gene2 <- index[[simulation]][, pair][2];

            # grab out the covariation strength:
            latent.covariance <- latent[[simulation]][gene1, gene2];

            # create a vector that enclose this information:
            latent.perturbation[simulation, pair] <- latent.covariance;

            # grab out the pruned edge:
            nscp.edge[simulation, pair] <- network[[simulation]][gene1, gene2];
            }

        # bin the real covariation values into 10 different quantiles:
        binned.simulation <- seq(1, ncol(latent.perturbation));
        names(binned.simulation) <- cut(
            latent.perturbation[simulation, ],
            breaks = bin.breaks,
            include.lowest = TRUE,
            labels = paste0('bin', seq(1, n.bin))
            );

        # find the number of pruned edges that still formed:
        for (bin in paste0('bin', seq(1, n.bin))){
            # grab out the bin index:
            bin.index <- binned.simulation[names(binned.simulation) == bin];
            
            # also grab out the number of gene pairs that are contained in each bin for latent perturbation:
            perturbed.pairs.bin[simulation, bin] <- table(names(binned.simulation))[bin];

            # fill in the edges classification matrix:
            nscp.classification[simulation, bin] <- sum(nscp.edge[simulation, bin.index] != 0);            
            }
        }
        nscp.classification[is.na(nscp.classification)] <- 0;
        perturbed.pairs.bin[is.na(perturbed.pairs.bin)] <- 0;

        # calculate power:
        all.sim.nscp <- colSums(nscp.classification);
        all.sim.latent <- colSums(perturbed.pairs.bin);
        power <- sapply(seq(1, n.bin), FUN = function(ix) all.sim.nscp[ix] / all.sim.latent[ix]);

        # return result:
        if (return == 'power') {
            return(power);
            }

        if (return == 'all') {
            collector <- list(nscp.classification, perturbed.pairs.bin, power);
            names(collector) <- c('nscp.classification', 'perturbed.pairs.bin', 'power')
            return(collector);
        }

    }