### FUNCTION LOADING ##############################################################################
# progress bar:
# write a function that prints out progress bar:
progress_bar <- function(current.step, total.step) {
    # first get the width of the progress bar according to the screen size:
    width <- options("width")$width / 2; # for split screen

    # calculate the percentage of the progress bar:
    percentage <- current.step / total.step;

    # find the width of the program bar to print:
    current.progress <- round(percentage * width);
    left <- width - current.progress;

    cat(
        # print the opening bracket:
        '\r', '| ',
        # print the * as current progress:
        paste(rep('*', current.progress), collapse = ''),
        # print the ' ' as blank for the left over progress:
        paste(rep(' ', left), collapse = ''), '|'
        );

    if (current.step == total.step) {
        cat('\n Prgress Complete! \n');
        }
    }
