### technical-variance-estimate.R #################################################################
# purpose: estimate technical variance estimate using methods

### COMPUTATION ###################################################################################
compute.technical.variance <- function(expression, span = 0.3, degree = 2, loggedCPM = TRUE) {
    
    if (loggedCPM) {
        log.cpm.variance <- apply(expression, 2, var);
        log.cpm.mean <- apply(expression, 2, mean);
        cpm.expression <- exp(expression) - 1;
        log.expression <- expression;
        stopifnot(signif(rowSums(cpm.expression), 10) == 1e6);
        }

    else {
        stop('log cpm expected! \n')
        }

    # compute the gene mean and variance at original space:
    cpm.mean <- apply(cpm.expression, 2, mean);
    cpm.variance <- apply(cpm.expression, 2, var);

    # check for log scaled mean/variance vs original scale
    stopifnot(cpm.variance > log.cpm.variance);

    # construct local polynomial fit using original scale cpm variance and mean:
    local.polynomial.model <- loess(
        formula = log10(cpm.variance) ~ log10(cpm.mean),
        span = span,
        degree = degree
        );

    # extract the predicted variance:
    predicted.log.variance <- local.polynomial.model$fit;
    predicted.variance <- 10^(predicted.log.variance);

    # compute the %technical variance:
    # defined as 10^fited value / original scale variance:
    tech.percentage <- predicted.variance / cpm.variance ;

    # compute the technical variance:
    # defined as log scaled variance * percentage:
    technical.noise <- log.cpm.variance * tech.percentage;
    return(tech.percentage);
    }
