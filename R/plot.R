#' 	Graphical representation of a RxC ecological inference (vote transfer) matrix
#'
#' @description Plot method for objects obtained with BPF.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @param x An object output of the **BPF** function.
#' @param margins A `TRUE/FALSE` argument informing whether the margins of the matrix should be displayed. Default, `TRUE`.
#' @param digits Integer indicating the number of decimal places to be shown. Default, 2.
#' @param row.names Names to be used for the rows of the matrix.
#' @param col.names Names to be used for the columns of the matrix.
#' @param size.numbers A reference number indicating the average font size to be used for the transfer numbers. Default, 6.
#' @param size.labels A number indicating the font size to be used for labels. Default, 4.
#' @param size.margins A number indicating the font size to be used for margin numbers. Default, 6.
#' @param colour.cells Background base colour for cells.
#' @param colour.grid Colour to be used for grid lines.
#' @param alpha A \[0,1\] number of colour transparency.
#' @param which A vector of integers informing the units for which the aggregate transfer matrix should be plotted. Default, `NULL`: the global matrix is shown.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param show.plot A `TRUE/FALSE` value indicating whether the plot should be displayed as a side-effect. By default, `TRUE`.
#'
#' @return
#' Invisibly returns the (ggplot) description of the plot, which is a list with components that contain the plot itself, the data, information about the scales, panels etc.
#'
#' @note ggplot2 is needed to be installed for this function to work.
#'
# @import ggplot2
#'
#' @export
#' @method plot BPF
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
#' example <- BPF(votes1, votes2)
#' p <- plot(example, show.plot = FALSE)
#' p
#'
plot.BPF <- function(x,
                     margins = TRUE,
                     digits = 2,
                     row.names = NULL,
                     col.names = NULL,
                     size.numbers = 6,
                     size.labels = 4,
                     size.margins = 6,
                     colour.cells = "darkolivegreen3",
                     colour.grid = "floralwhite",
                     alpha = 0.5,
                     which = NULL,
                     ...,
                     show.plot = TRUE){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }

  prop <- x$TM*100
  votos <- x$TM.votes

  n.fil <- nrow(prop)
  n.col <- ncol(prop)

  if (!is.null(which)){
    if(is.null(x$TR.votes.units)){
      stop("Unit matrices are not available in the input object, please set 'which = NULL'")
    }
    if (!all(which %in% 1L:dim(x$TR.votes.units)[3L])){
      stop("The 'which' argument that you are using is not valid. Please check it.")
    }
    x$TM.votes <- apply(x$TM.votes.units[, , which], c(1, 2), sum)
    votos <- x$TM.votes[1L:n.fil, 1L:n.col]
    prop <- (x$TM.votes/rowSums(x$TM.votes)*100)[1L:n.fil, 1L:n.col]
  }

  votos.filas <- rowSums(x$TM.votes)[1L:n.fil]
  votos.columnas <- colSums(x$TM.votes)[1L:n.col]
  prop.filas <- format(round(rowSums(x$TM.votes)[1L:n.fil]/sum(votos.filas)*100,
                             digits), n.small = digits)
  prop.columnas <- format(round(colSums(x$TM.votes)[1L:n.col]/sum(votos.columnas)*100,
                                digits), n.small = digits)

  prop2 <- as.vector(prop)

  ## base de trabajo
  bbdd <- cbind(expand.grid(n.fil:1L, 1L:n.col), prop2,
                format(round(prop2, digits), n.small = digits))
  bbdd <- as.data.frame(bbdd)

  names(bbdd) <- c("y", "x", "coefficient", "label")
  bbdd$color <- paste0("gray", round((100 - round(bbdd$coefficient))/1.5))

  ## Tamanyos numeros
  factor.size <- log(votos/sum(votos)*100 + 1L)
  factor.size <- factor.size/max(max(factor.size)) + 0.5
  bbdd$size <- as.vector(factor.size*size.numbers)

  ## Se añaden marginales
  if (margins){
    suma.fila <- data.frame(y = n.fil:1L, x = n.col + 1L, coefficient = 0,
                            label = prop.filas, color = "gray27",
                            size = size.margins)
    suma.columna <- data.frame(y = 0, x = 1L:n.col, coefficient = 0,
                               label = prop.columnas, color = "gray27",
                               size = size.margins)
    bbdd <- rbind(bbdd, suma.fila, suma.columna)
  }

  ## Se añaden nombres
  if (is.null(row.names)){
    row.names <- rownames(prop)
  }
  nombres.fila <- data.frame(y = n.fil:1L, x = 0, coefficient = 0,
                             label = row.names, color = "gray27", size = size.labels)

  if (is.null(col.names)){
    col.names <- colnames(prop)
  }
  nombres.columna <- data.frame(y = n.fil + 1L, x = 1L:n.col, coefficient = 0,
                                label = col.names, color = "gray27", size = size.labels)

  bbdd <- rbind(bbdd, nombres.fila, nombres.columna)

  p <- ggplot2::ggplot(bbdd, ggplot2::aes(x = !!quote(x), y = !!quote(y))) +
    ggplot2::geom_tile(ggplot2::aes(fill = !!quote(coefficient)),
                       color = colour.grid) +
    ggplot2::scale_fill_continuous(high = scales::alpha(colour = colour.cells, alpha = alpha),
                                   low = "white", trans = "sqrt") +
    ggplot2::geom_text(ggplot2::aes(label = !!quote(label)),
                       size = bbdd$size, colour = bbdd$color) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  if (show.plot) print(p)
  return(p)
}
