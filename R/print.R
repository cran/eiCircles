#' 	Print a summary of an output of the BPF function
#'
#' @description Print method for objects obtained with the BPF function.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param x An object output of the **BPF** function.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param margins A `TRUE/FALSE` argument informing if the margins of the transition matrix should be displayed. Default, `TRUE`.
#' @param digits Integer indicating the number of decimal places to be shown. Default, 2.
#'
#' @return
#' {No return value, called for side effects.}
#'
#' @export
#' @method print BPF
#' @examples
#' votes1 <- structure(list(P1 = c(16L, 4L, 13L, 6L, 1L, 16L, 6L, 17L, 48L, 14L),
#'                          P2 = c(8L, 3L, 0L, 5L, 1L, 4L, 7L, 6L, 28L, 8L),
#'                          P3 = c(38L, 11L, 11L, 3L, 13L, 39L, 14L, 34L, 280L, 84L),
#'                          P4 = c(66L, 5L, 18L, 39L, 30L, 57L, 35L, 65L, 180L, 78L),
#'                          P5 = c(14L, 0L, 5L, 2L, 4L, 21L, 6L, 11L, 54L, 9L),
#'                          P6 = c(8L, 2L, 5L, 3L, 0L, 7L, 7L, 11L, 45L, 17L),
#'                          P7 = c(7L, 3L, 5L, 2L, 3L, 17L, 7L, 13L, 40L, 8L)),
#'                          row.names = c(NA, 10L), class = "data.frame")
#' votes2 <- structure(list(C1 = c(2L, 1L, 2L, 2L, 0L, 4L, 0L, 4L, 19L, 14L),
#'                          C2 = c(7L, 3L, 1L, 7L, 2L, 5L, 3L, 10L, 21L, 6L),
#'                          C3 = c(78L, 7L, 28L, 42L, 28L, 84L, 49L, 85L, 260L, 100L),
#'                          C4 = c(56L, 14L, 20L, 7L, 19L, 54L, 22L, 50L, 330L, 91L),
#'                          C5 = c(14L, 3L, 6L, 2L, 3L, 14L, 8L, 8L, 45L, 7L)),
#'                          row.names = c(NA, 10L), class = "data.frame")
#' example <- BPF(votes1, votes2, local = "none")
#' print(example, digits = 1, margins = TRUE)
#'

#' @export
print.BPF  <- function(x,
                         ...,
                         margins = TRUE,
                         digits = 2)
{

  print.summary.BPF(x = summary.BPF(object = x),
                        margins = margins,
                        digits = digits)

}

#' 	Print a summary of a summary.BPF object
#'
#' @description Print method for `summary.BPFC` objects
#' @inheritParams print.BPF
#' @param x An `summary.BPF` class object.
#' @return
#' {No return value, called for side effects.}
#' @method print summary.BPF
#' @export
print.summary.BPF  <- function(x,
                                 ...,
                                 margins = TRUE,
                                 digits = 2)
{

    tabla <- format(round(x$prop.matrix*100, digits), nsmall = digits)
    tabla <- apply(tabla, 2, as.character)
    rownames(tabla) <- rownames(x$prop.matrix)

    if (margins){
      nr <- nrow(tabla)
      tabla <- rbind(tabla, format(round(x$col.margins[1L:ncol(tabla)]*100, digits), nsmall = digits))
      tabla <- cbind(tabla, c(format(round(x$row.margins[1L:nr]*100, digits), justify = "right", nsmall = digits), ""))
    }

    cat("Estimated row-standardized transfer matrix \n")
    print(as.table(tabla))

}
