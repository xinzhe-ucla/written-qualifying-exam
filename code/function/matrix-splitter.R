# function:
matrix.spliter <- function(data, fold.number) {
    # PURPOSE:
    # split the data by columns into different fold (return indexing)

    # input:
    # data: matrix that we want to split
    # fold.number: number of fold to split the data into

    # output:
    # collector: index in each fold in a nested list

    # inititate collector for return:
    collector <- vector('list', length = fold.number);

    # generate index for the columns based on number of folds required:
    fold.index <- floor(seq(from = 0, to = ncol(data), length.out = fold.number + 1));

    # split data into requried folds:
    for (fold in seq(1, fold.number)){
        # grab the columns that will be in this fold:
        features.in.fold <- colnames(data)[(fold.index[fold] + 1):(fold.index[fold + 1])];

        # output the splitted data:
        collector[[fold]] <- features.in.fold;
        }
    names(collector) <- paste0('fold', seq(1, fold.number));

    # add in check to make sure that splitting is sensible:
    stopifnot(identical(Reduce('c', collector), colnames(data)));
    return(collector);
    }