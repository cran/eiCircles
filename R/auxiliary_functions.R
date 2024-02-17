## Auxiliary functions

## function tests_inputs ##
tests_inputs <- function(argg){

  # Tests numeric data
  x <- as.matrix(argg$X)
  y <- as.matrix(argg$Y)
  if (nrow(x) != nrow(y))
    stop('The number of units (rows) is different in arguments "X" and "Y".')
  if (min(x,y) < 0L) stop('Negative values are not allowed in arguments "X" and "Y".')

  # Test census.changes
  if (!(argg$census.changes %in% c("adjust1", "adjust2")))
    stop("The value set for argument 'census.changes' is invalid. Only 'adjust1' and 'adjust2' are allowed.")

  # Test local
  if (!(argg$local %in% c("none", "IPF", "lik")))
    stop("The value set for argument 'local' is invalid. Only 'none', 'IPF' and 'lik' are allowed.")

  # Test stability.par
  if (argg$stability.par < 0)
    stop("The argument 'stability.par' must be non-negative.")

  # Test confidence
  if(argg$confidence >= 1L | argg$confidence <= 0L)
    stop('Only values in the interval (0, 1) are allowed for the "confidence" argument.')

  # Test max.iter
  if(argg$max.iter < 0 | abs(argg$max.iter - round(argg$max.iter)) > 0)
    stop('Only positive integers are allowed for the "max.iter" argument.')

  # Test covariates
  covars <- argg$covariates
  if(!is.null(argg$covariates)){
    if (length(covars) != 2L) {
      stop('The object used as argument "covariates" does not have the proper structure. It should have two components.')
    }
    covars[[2L]] <- to_matrix(covars[[2L]])
    if (nrow(covars[[1L]]) != nrow(x)) {
      stop('The object used as argument "covariates" does not have the proper structure. The number of units for which covariates are defined is different to the number of units (rows) in arguments "X" and "Y".')
    }
    if (ncol(covars[[2L]]) != 3L) {
      stop('The object used as argument "covariates" does not have the proper structure. The matrix of metadata "meta" should have three columns.')
    }
    if (nrow(covars[[2L]]) < 1L) {
      stop('The object used as argument "covariates" does not have the proper structure. The matrix of metadata "meta" should have at least a row.')
    }
    n.covariate <- text2number(colnames(covars[[1L]]), covars[[2L]][, 3L])
    if (sum(is.na(n.covariate)) > 0L){
      stop('The object used as argument "covariates" does not have the proper structure. At least a name declared in "meta" is not in "covar".')
    } else {
      covars[[2L]][, 3L] <- n.covariate
    }
    n.rows <- text2number(colnames(x), covars[[2L]][, 1L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "covariates" does not have the proper structure. At least a name declared in "meta" is not in "X".')
    } else {
      covars[[2L]][, 1L] <- n.rows
    }
    n.cols <- text2number(colnames(y), covars[[2L]][, 2L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "covariates" does not have the proper structure. At least a name declared in "meta" is not in "Y".')
    } else {
      covars[[2L]][, 2L] <- n.cols
    }
    covars[[2L]] <- to_matrix(suppressWarnings(apply(covars[[2L]], 2, as.numeric)))
    invalid <- test_n(ncol(x), covars[[2L]][, 1L]) + test_n(ncol(y), covars[[2L]][, 2L]) +
      test_n(ncol(covars[[1L]]), covars[[2L]][, 3L])
    if (invalid > 0L)
      stop('The object used as argument "covariates" is not properly defined. At least a number included in "meta" is invalid.')
  }


  # Test null.cells
  null.cells <- to_matrix(argg$null.cells)
  if(!is.null(null.cells)){
    if (ncol(null.cells) != 2L) {
      stop('The object used as argument "null.cells" does not have the proper structure. It should have two columns.')
    }
    if (nrow(null.cells) < 1L) {
      stop('The object used as argument "null.cells" does not have the proper structure.')
    }
    n.rows <- text2number(colnames(x), null.cells[, 1L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "null.cells" is not properly defined. At least a name declared is not in "X".')
    } else {
      null.cells[, 1L] <- n.rows
    }
    n.cols <- text2number(colnames(y), null.cells[, 2L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "null.cells" is not properly defined. At least a name declared is not in "Y".')
    } else {
      null.cells[, 2L] <- n.cols
    }
    null.cells <- to_matrix(suppressWarnings(apply(null.cells, 2L, as.numeric)))
    invalid <- test_n(ncol(x), null.cells[, 1L]) + test_n(ncol(y) - 1L, null.cells[, 2L])
    if (invalid > 0L)
      stop('The object used as argument "null.cells" is not properly defined. At least a number included is invalid.')
  }


  # Test row.cells.relationships
  row.cells.relationships <- to_matrix(argg$row.cells.relationships)
  if(!is.null(row.cells.relationships)){
    if (ncol(row.cells.relationships) != 4L) {
      stop('The object used as argument "row.cells.relationships" does not have the proper structure. It should have four columns.')
    }
    if (nrow(row.cells.relationships) < 1L) {
      stop('The object used as argument "row.cells.relationships" does not have the proper structure.')
    }
    n.rows <- text2number(colnames(x), row.cells.relationships[, 1L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "row.cells.relationships" is not properly defined. At least a name declared is not in "X".')
    } else {
      row.cells.relationships[, 1L] <- n.rows
    }
    n.cols <- text2number(colnames(y), row.cells.relationships[, 2L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "row.cells.relationships" is not properly defined. At least a name declared is not in "Y".')
    } else {
      row.cells.relationships[, 2L] <- n.cols
    }
    n.cols <- text2number(colnames(y), row.cells.relationships[, 3L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "row.cells.relationships" is not properly defined. At least a name declared is not in "Y".')
    } else {
      row.cells.relationships[, 3L] <- n.cols
    }
    row.cells.relationships <- to_matrix(suppressWarnings(apply(row.cells.relationships, 2L, as.numeric)))
    invalid <- test_n(ncol(x), row.cells.relationships[, 1L]) + test_n(ncol(y) - 1L, row.cells.relationships[, 2L]) +
      test_n(ncol(y) - 1L, row.cells.relationships[, 3L]) + test_p_real(row.cells.relationships[, 4L])
    if (invalid > 0L)
      stop('The object used as argument "row.cells.relationships" is not properly defined. At least a number included is invalid.')
  }


  # Test row.cells.relationships.C
  row.cells.relationships.C <- to_matrix(argg$row.cells.relationships.C)
  if(!is.null(row.cells.relationships.C)){
    if (ncol(row.cells.relationships.C) != 3L) {
      stop('The object used as argument "row.cells.relationships.C" does not have the proper structure. It should have three columns.')
    }
    if (nrow(row.cells.relationships.C) < 1L) {
      stop('The object used as argument "row.cells.relationships.C" does not have the proper structure.')
    }
    n.rows <- text2number(colnames(x), row.cells.relationships.C[, 1L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "row.cells.relationships.C" is not properly defined. At least a name declared is not in "X".')
    } else {
      row.cells.relationships.C[, 1L] <- n.rows
    }
    n.cols <- text2number(colnames(y), row.cells.relationships.C[, 2L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "row.cells.relationships.C" is not properly defined. At least a name declared is not in "Y".')
    } else {
      row.cells.relationships.C[, 2L] <- n.cols
    }
    row.cells.relationships.C <- to_matrix(suppressWarnings(apply(row.cells.relationships.C, 2, as.numeric)))
    invalid <- test_n(ncol(x), row.cells.relationships.C[, 1L]) + test_n(ncol(y) - 1L, row.cells.relationships.C[, 2L]) +
      test_p_real(row.cells.relationships.C[, 3L])
    if (invalid > 0L)
      stop('The object used as argument "row.cells.relationships.C" is not properly defined. At least a number included is invalid.')
  }


  # Test pair.cells.relationships
  pair.cells.relationships  <- to_matrix(argg$pair.cells.relationships)
  if(!is.null(pair.cells.relationships)){
    if (ncol(pair.cells.relationships) != 7L) {
      stop('The object used as argument "pair.cells.relationships" does not have the proper structure. It should have seven columns.')
    }
    if (nrow(pair.cells.relationships) < 1L) {
      stop('The object used as argument "pair.cells.relationships" does not have the proper structure.')
    }
    n.rows <- text2number(colnames(x), pair.cells.relationships[, 1L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "pair.cells.relationships" is not properly defined. At least a name declared is not in "X".')
    } else {
      pair.cells.relationships[, 1L] <- n.rows
    }
    n.cols <- text2number(colnames(y), pair.cells.relationships[, 2L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "pair.cells.relationships" is not properly defined. At least a name declared is not in "Y".')
    } else {
      pair.cells.relationships[, 2L] <- n.cols
    }
    n.cols <- text2number(colnames(y), pair.cells.relationships[, 3L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "pair.cells.relationships" is not properly defined. At least a name declared is not in "Y".')
    } else {
      pair.cells.relationships[, 3L] <- n.cols
    }
    n.rows <- text2number(colnames(x), pair.cells.relationships[, 4L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "pair.cells.relationships" is not properly defined At least a name declared is not in "X".')
    } else {
      pair.cells.relationships[, 4L] <- n.rows
    }
    n.cols <- text2number(colnames(y), pair.cells.relationships[, 5L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "pair.cells.relationships" is not properly defined. At least a name declared is not in "Y".')
    } else {
      pair.cells.relationships[, 5L] <- n.cols
    }
    n.cols <- text2number(colnames(y), pair.cells.relationships[, 6L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "pair.cells.relationships" is not properly defined. At least a name declared is not in "Y".')
    } else {
      pair.cells.relationships[, 6L] <- n.cols
    }
    pair.cells.relationships <- to_matrix(suppressWarnings(apply(pair.cells.relationships, 2, as.numeric)))
    invalid <- test_n(ncol(x), pair.cells.relationships[, 1L]) + test_n(ncol(y), pair.cells.relationships[, 2L]) +
      test_n(ncol(y), pair.cells.relationships[, 3L]) + test_n(ncol(x), pair.cells.relationships[, 4L]) +
      test_n(ncol(y), pair.cells.relationships[, 5L]) + test_n(ncol(y), pair.cells.relationships[, 6L]) +
      test_p_real(pair.cells.relationships[, 7L])
    if (invalid > 0L)
      stop('The object used as argument "pair.cells.relationships" is not properly defined. At least a number included is invalid.')
  }

  # Test cells.fixed.logit
  cells.fixed.logit <- to_matrix(argg$cells.fixed.logit)
  if(!is.null(cells.fixed.logit)){
    if (ncol(cells.fixed.logit) != 3L) {
      stop('The object used as argument "cells.fixed.logit" does not have the proper structure. It should have three columns.')
    }
    if (nrow(cells.fixed.logit) < 1L) {
      stop('The object used as argument "cells.fixed.logit" does not have the proper structure.')
    }
    n.rows <- text2number(colnames(x), cells.fixed.logit[, 1L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "cells.fixed.logit" is not properly defined. At least a name declared is not in "X".')
    } else {
      cells.fixed.logit[, 1L] <- n.rows
    }
    n.cols <- text2number(colnames(y), cells.fixed.logit[, 2L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "cells.fixed.logit" is not properly defined. At least a name declared is not in "Y".')
    } else {
      cells.fixed.logit[, 2L] <- n.cols
    }
    cells.fixed.logit <- to_matrix(suppressWarnings(apply(cells.fixed.logit, 2, as.numeric)))
    invalid <- test_n(ncol(x), cells.fixed.logit[, 1L]) + test_n(ncol(y), cells.fixed.logit[, 2L]) +
      test_real(row.cells.relationships.C[, 3L])
    if (invalid > 0L)
      stop('The object used as argument "cells.fixed.logit" is not properly defined. At least a number included is invalid.')

  }


  # Test dispersion.rows
  dispersion.rows <- to_matrix(argg$dispersion.rows)
  if(!is.null(dispersion.rows)){
    if (ncol(dispersion.rows) != 2L) {
      stop('The object used as argument "dispersion.rows" does not have the proper structure. It should have two columns.')
    }
    if (nrow(dispersion.rows) < 1L) {
      stop('The object used as argument "dispersion.rows" does not have the proper structure.')
    }
    n.rows <- text2number(colnames(x), dispersion.rows[, 1L])
    if (sum(is.na(n.rows)) > 0L){
      stop('The object used as argument "dispersion.rows" is not properly defined. At least a name declared is not in "X".')
    } else {
      dispersion.rows[, 1L] <- n.rows
    }
    n.cols <- text2number(colnames(x), dispersion.rows[, 2L])
    if (sum(is.na(n.cols)) > 0L){
      stop('The object used as argument "dispersion.rows" is not properly defined. At least a name declared is not in "Y".')
    } else {
      dispersion.rows[, 2L] <- n.cols
    }
    dispersion.rows <- to_matrix(suppressWarnings(apply(dispersion.rows, 2, as.numeric)))
    invalid <- test_n(ncol(x), dispersion.rows[, 1L]) + test_n(ncol(x), dispersion.rows[, 2L])
    if (invalid > 0L)
      stop('The object used as argument "dispersion.rows" is not properly defined. At least a number included is invalid.')
    dispersion.rows <- t(apply(dispersion.rows, 1L, sort))
    dispersion.rows <- dispersion.rows[(dispersion.rows[, 2L] - dispersion.rows[, 1L]) != 0L, ]
  }

  return(list("x" = x, "y" = y, "covars" = covars,
              "null.cells" = null.cells,
              "row.cells.relationships" = row.cells.relationships,
              "row.cells.relationships.C" = row.cells.relationships.C,
              "pair.cells.relationships" = pair.cells.relationships,
              "cells.fixed.logit" = cells.fixed.logit,
              "dispersion.rows" = dispersion.rows)
  )
}

##------------------------------------------------------------------------------

## function to_matrix ##
to_matrix <- function(x){
  if(is.null(x)){
    return(x)
  } else {
    x <- as.matrix(x)
    if(ncol(x) == 1L) x <- t(x)
    return(x)
  }
}
##------------------------------------------------------------------------------

## function text2number ##
# This function transform in numbers the definition of variables or options introduced by the user.
text2number <- function(vector1, vector2){
  vector2.t <- suppressWarnings(as.numeric(vector2))
  vector2.t[is.na(vector2.t)] <- match(vector2, vector1)[is.na(vector2.t)]
  return(vector2.t)
}
##------------------------------------------------------------------------------

## functions test_rows, test_columns, test_numbers ##
test_n <- function(n, x){
  out <- sum(!(x %in% 1L:n))
  return(out)
}

test_real <- function(x){
  out <- sum(is.na(is.numeric(x)))
  return(out)
}

test_p_real <- function(x){
  out <- sum(is.na(x < 0L)) + sum(x < 0L & !is.na(x))
  return(out)
}
##------------------------------------------------------------------------------

## function preprocessing() ##
# This function selects stable units and adjust census changes
preprocessing <- function(X0, Y0, stable.units, stability.par, census.changes){
  if(stable.units){
    r <- abs(log(rowSums(Y0) / rowSums(X0)))
    units.selected <- which(r <= stability.par) # = in case zero be selected
    N <- X0[units.selected, ]
    Y <- Y0[units.selected, ]
    # plot(r) # Maybe this should be parametrized!!! and decided if plotted
  } else{
    units.selected <- 1:nrow(X0)
    N <- X0
    Y <- Y0
  }

  # Redimension and re-scale
  if(census.changes == "adjust1"){
    dn <- rowSums(N)
    dy <- rowSums(Y)
    N <- dy / dn * N
  } else if (census.changes == "adjust2"){
    dn <- rowSums(N)
    dy <- rowSums(Y)
    Y <- dn / dy * Y
  }

  output <- list("N" = N, "Y" = Y, "units.selected" = units.selected)
  return(output)
}
##------------------------------------------------------------------------------

## Function to_matrix_constraint() ##
# This function is an auxiliary function for converting the constraint objects into proper matrices
to_matrix_constraint <- function(x){
  output <- x
  if (is.null(x)){
    output <- matrix(0L, nrow = 0L, ncol = 0L)
  } else if (is.null(nrow(x))){
    output <- matrix(x, nrow = 1L)
  }
  return(output)
}
##------------------------------------------------------------------------------

## Function to_two_row_matrix_constraint() ##
# This function is an auxiliary function for converting type 5 constraints into a proper matrix where each constraint takes 4 cells on two consecutive rows of the matrix
to_two_row_matrix_constraint <- function(x){
  if (is.null(x)){
    output <- matrix(0, nrow = 0L, ncol = 0L)
  } else if (is.null(nrow(x))){
    output <- matrix(c(x[1L:3L], x[7L], x[4L:6L], 1L),
                     nrow = 2L, byrow = TRUE)
  } else {
    nf <- nrow(x)
    output <- matrix(1L, nrow = 2L*nf, ncol = 4L)
    for (ff in 1L:nf){
      output[2L*(ff - 1L) + 1L, 1L:3L] <- x[ff, 1L:3L]
      output[2L*(ff - 1L) + 1L, 4L] <- x[ff, 7L]
      output[2L*ff, 1L:3L] <- x[ff, 4L:6L]
    }
  }
  return(output)
}
##------------------------------------------------------------------------------

## Function Imod_matrix() ##
# This functions performs the work relating creating the Imod matrix to handle constraints and covariates
Imod_matrix <- function(I, covariates,
                        null.cells,
                        null.cells.C = NULL,
                        row.cells.relationships,
                        row.cells.relationships.C,
                        pair.cells.relationships,
                        cells.fixed.logit,
                        dispersion.rows,
                        dispersion.row = NULL,
                        dispersion.fixed = NULL){

  # Constraints and covariates

  ########################################################
  # Definition of the auxilar Imod matrix
  # constraints: conventions for Imod
  #     %- cod, row, col1, col2, log(a).
  #     %-> 0: no restriction           => [0 0 0 0 0]
  #     %-> 1: p(ij) = 0                => [1 i j 0 0]
  #     %-> 2: p(iJ) = 0                => [2 i j 0 0] -> j=new ref cat # NOT ACTIVE
  #     %-> 3: p(ik)=a*p(ij)            => [3 i j k log(a)]
  #     %-> 4: p(ij) = a p(iJ)          => [4 i j 0 log(a)]
  #     %-> 5: p(hj)/p(hk) = a p(im)/p(in)
  #         %                           => [5 h j k log(a)
  #         %                           => [5 i m n 0]
  #     %-> 6: la(i) = la(h)            => [6 0 i h 0]
  #     %-> 7: la(i)=0                  => [7 i 0 0 0] # NOT INCLUDED
  #     %-> 8: covariate dicotomicche   => [8 i j 1 0 #1=code
  #         %  o continue che agiscono      8 h k 1 0 #2=row
  #         %  su certi sottogruppi di      8 l m 1 0 #3=col
  #         %  caselle, per "group"         ......... #4=group
  #         %  si intende la colonna        8 r s 2 0
  #         %  della matrice C delle        8 t v 2 0
  #         %  covariate                    .........
  #     %-> 9: fix logit of cell i, j => [9 i j 0 logit ]
  #     %->10: fix la(i)            => [10 i 0 0 valor] # NOT INCLUDED


  # Definition of the auxiliary matrix Imod
  Imod <- NULL

  # 1: cells constrained to be 0
  # Ize <- null.cells
  null.cells <- to_matrix_constraint(null.cells)
  ir <- nrow(null.cells)
  if (ir > 0L) {
    Imod <- rbind(
      Imod,
      cbind(rep(1L, ir), null.cells, matrix(0L, nrow = ir, ncol = 2L))
    )
  }

  # 2: Reference option constrained to be 0: NOT ACTIVE
  # Inz <- null.cells.C
  null.cells.C <- to_matrix_constraint(null.cells.C)
  ir <- nrow(null.cells.C)
  if (ir > 0L) {
    Imod <- rbind(
      Imod,
      cbind(rep(2L, ir), null.cells.C, matrix(0L, nrow = ir, ncol = 2L))
    )
  }

  # 3: constrain relationships between cells in same row
  Icacb <- to_matrix_constraint(row.cells.relationships)
  ir <- nrow(Icacb )
  if (ir > 0L) {
    Icacb[, 4L] <- log(Icacb[, 4L])
    Imod <- rbind(
      Imod,
      cbind(rep(3L, ir), Icacb)
    )
  }

  # 4: constrain relationships with last column reference
  Iac <- to_matrix_constraint(row.cells.relationships.C)
  ir <- nrow(Iac)
  if (ir > 0) {
    Iac <- cbind(Iac, log(Iac[, 3L]))
    Iac[, 3L] <- 0L
    Imod <- rbind(
      Imod,
      cbind(rep(4L, ir), Iac)
    )
  }

  # 5: constraints among 4 cells on two different rows
  Irarb <- to_two_row_matrix_constraint(pair.cells.relationships)
  ir <- nrow(Irarb)
  if (ir > 0) {
    Irarb[, 4L] <- log(Irarb[, 4L])
    Imod <- rbind(
      Imod,
      cbind(rep(5L, ir), Irarb)
    )
  }

  # 9: fixed the logit of selected cells
  IcaL <- to_matrix_constraint(cells.fixed.logit)
  ir <- nrow(IcaL)
  if (ir > 0L) {
    IcaL <- cbind(IcaL, IcaL[, 3L])
    IcaL[, 3L] <- 0L
    Imod <- rbind(
      Imod,
      cbind(rep(9L, ir), IcaL)
    )
  }

  # 6: constant over-dispersions
  Overd <- to_matrix_constraint(dispersion.rows)
  ir <- nrow(Overd)
  if (ir > 0L) {
    Overd <- cbind(rep(6L, ir), rep(0L, ir), Overd, rep(0L, ir))
    Imod <- rbind(
      Imod,
      Overd
    )
  }

  # Covariates
  Ico <- to_matrix_constraint(covariates[[2L]])
  ir <- nrow(Ico)
  if (ir > 0L) {
    Ico <- cbind(Ico, 1L:ir, matrix(0, nrow = ir, ncol = 1))
    Ico <- Ico[, -3L, drop = FALSE]
    Imod <- rbind(
      Imod,
      cbind(rep(8L, ir), Ico)
    )
  }

  return(Imod)

}
##------------------------------------------------------------------------------

#### Function model_building
model_building <- function(N, Y, Imod, start, cvar, cs, seed = NULL) {

  # Preliminaries
  I <- ncol(N)
  J <- ncol(Y)
  Jm <- ncol(Y) - 1L
  Ym <- Y[, 1L:Jm]
  IJm <- I * Jm
  off <- rep(0L, IJm + I)
  X <- diag(IJm + I)
  Z <- matrix(0L, nrow = IJm + I, ncol = 0L)
  ibet <- c()
  if(!is.null(seed)) set.seed(seed)

  if (is.null(start)) {
    # Initial estimates for beta
    La <- rep(1L, I)/100
    P <- kronecker(rep(1L, I), t(colSums(Y) / sum(Y))) +
      matrix(stats::runif(I*J, 0L, 1L), nrow = I, ncol = J) / 50
    P <- P * (P > 10^(-8)) + 10^(-8) * (P <= 10^(-8))
    beta <- log(P[, 1L:Jm]/P[, J])
    beta <- c(t(beta), log(La / (1L - La)))
    beta <- beta + stats::rnorm(IJm + I) / 50
  } else {
    beta <- start
  }

  if(!is.null(Imod)){
    # 1: a transition set to 0
    Im <- which(Imod[, 1L] == 1L)
    ni <- length(Im)
    if (ni > 0L){
      Im <- Imod[Im, , drop = FALSE]
      for (i in 1L:ni) {
        hh <- (Im[i, 2L] - 1L) * Jm + Im[i, 3L]
        ibet <- c(ibet, hh)
        off[hh] <- -20
      }
    }

    # 2: a reference set to 0 ## NOT ACTIVE
    Im <- which(Imod[, 1L] == 2L)
    ni <- length(Im)
    if (ni > 0L){
      Im <- Imod[Im, , drop = FALSE]
      for (i in 1L:ni) {
        h1 <- (Im[i, 2L] - 1L) * Jm
        h2 <- (h1 + 1L):(h1 + Jm)
        h1 <- h1 + Im[i, 3L]
        off[h2] <- off[h2] + 20L
        ibet <- c(ibet, h1)
      }
    }

    # 3:  p(ik) = a*p(ij)
    Im <- which(Imod[, 1L] == 3L)
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      for (i in 1L:ni) {
        h1 <- (Im[i, 2L] - 1L) * Jm + Im[i, 3L]
        h2 <- (Im[i, 2L] - 1L) * Jm + Im[i, 4L]
        X[, h1] <- X[, h1] + X[, h2]
        ibet <- c(ibet, h2)
        off[h2] <- Im[i, 5L]
      }
    }

    # 4:  p(ij) = a*p(iJ)
    Im <- which(Imod[, 1L] == 4L)
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      for (i in 1L:ni) {
        hh <- (Im[i, 2L] - 1L) * Jm + Im[i, 3L]
        ibet <- c(ibet, hh)
        off[hh] <- Im[i, 5L]
      }
    }

    # 5:  p(hj)/p(hk) = a*p(im)/p(in)
    Im <- which(Imod[, 1L] == 5L)
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      ni <- ni / 2
      hh <- 2L * (1L:ni)
      Im <- cbind(Im[hh - 1L, , drop = FALSE], Im[hh, , drop = FALSE])
      for (i in 1L:ni) {
        h1 <- (Im[i, 2L] - 1L) * Jm + Im[i, 3L]
        h2 <- (Im[i, 2L] - 1L) * Jm + Im[i, 4L]
        k1 <- (Im[i, 7L] - 1L) * Jm + Im[i, 8L]
        k2 <- (Im[i, 7L] - 1L) * Jm + Im[i, 9L]
        X[, h1] <- X[, h1] - X[, k2]
        X[, h2] <- X[, h2] + X[, k2]
        X[, k1] <- X[, k1] + X[, k2]
        ibet <- c(ibet, k2)
        off[k2] <- Im[i, 5L]
      }
    }

    # 6: equal overdispersion
    Im <- which(Imod[, 1L] == 6L)
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      for(i in 1:ni) {
        h1 <- IJm + Im[i, 3L]
        h2 <- IJm + Im[i, 4L]
        X[, h1] <- X[, h1, drop = FALSE] + X[, h2, drop = FALSE]
        ibet <- c(ibet, h2)
      }
    }

    # 7: la(i) set to 0 # NOT ACTIVE
    Im <- which(Imod[, 1L] == 7L)
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      for(i in 1:ni) {
        hh <- IJm + Im[i, 2L]
        ibet <- c(ibet, hh)
        off[hh] <- -20
      }
    }

    # 8: effects of covariates
    Im <- which(Imod[, 1L] == 8L, )
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      ng <- Im[ni, 4L]
      Z <- matrix(0L, nrow = IJm + I, ncol = ng)
      for (i in 1L:ng) {
        id1 <- c()
        if (Im[i, 3L] <= Jm) {
          # only one logit in the row
          id1 <- c(id1, (Im[i, 2L] - 1L) * Jm + Im[i, 3L])
        } else {
          # all logits in the row
          id1 <- c(id1, ((Im[i, 2L] - 1L) * Jm + 1L):(Im[i, 2L] * Jm))
        }
        Z[id1, i] <- 1L
      }
    }

    # 9: fix the logit of a transition
    Im <- which(Imod[, 1L] == 9L)
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      for(i in 1:ni) {
        hh <- (Im[i, 2] - 1) * Jm + Im[i, 3]
        off[hh] <- Im[i, 5]
      }
    }

    # 10: fix overdispersion # NOT ACTIVE
    Im <- which(Imod[, 1L] == 10)
    ni <- length(Im)
    if (ni > 0L) {
      Im <- Imod[Im, , drop = FALSE]
      for(i in 1:ni) {
        hh <- IJm + Im[i, 2L]
        ibet <- c(ibet, hh)
        off[hh] <- log(Im[i, 5L]) # This is in the log scale
      }
    }
  }

  if (!is.null(ibet)){
    ier <- as.numeric(names(table(ibet))[table(ibet) > 1])
    if (length(ier) > 0L){
      fila <- ifelse( ier/Jm - floor(ier/Jm) == 0L, floor(ier/Jm), floor(ier/Jm) + 1L)
      col <- ier - (fila - 1L)*Jm
      fila.p <- fila <= I
      fila.o <- col[fila > I]
      fila.p <- paste0("(", fila[fila.p], ", ", col[fila.p], ")")
      fila.o <- ifelse(fila > I, paste0(" Constraints corresponding to over-dispersions in which the row(s) ", fila.o, ' is/are involved could be not consistent.'), '')
      fila.p <- ifelse(fila <= I, paste0('Constraints in which the cell(s) ', fila.p, ' of the transition matrix is/are involved could be not consistent.'), '')
      message <- paste0(fila.p, fila.o)
      warning(message)
    }
    # Adapt beta and design matrix
    # remove constrained elements
    X <- X[, -ibet]
    #k <- nrow(Z)
    z <- ncol(Z)
    if (is.null(start)) {
      beta <- beta[-ibet]
      if (cvar > 0L) {
        beta <- c(beta, rep(0L, z))
      }
    }
  } else {
    if (is.null(start) & cvar > 0L) {
      z <- ncol(Z)
      beta <- c(beta, rep(0L, z))
    }
  }

  if (cvar == 0L) {
    Z <- matrix(0, nrow = nrow(X), ncol = 0)
  }

  output <- list("X1" = X, "X2" = Z, "off" = off, "cs" = cs, "beta" = beta)
  return(output)
}
##------------------------------------------------------------------------------

## Function votlan()##
votlan <- function(N, Y, C, Imod, start, cs, seed, mit, tol, verbose, save.beta) {
  # N: results in election 1
  # Y: results ion election 2
  # C: matrix of covariates
  # Imod: desing matrix, output of Imod_function
  # start: matrix of initial estimates for the beta vector
  # mit: maximum number of iterations?

  #  -  MLE of vote transitions (Forcina, Marcehtti, CLDAG 2010)
  #     rows of N and Y are electoral results at time t-1 and t
  #     if model is constant across units, beta can be omitted and
  #     initial estimates of P and La provided; last col of P is
  #     used as reference; log-oversdipersion; C contains covariates
  #     or has 0 columns
  #  -  calls  ->  fiscor

  # Preliminaries
  # J <- ncol(Y)
  # I <- ncol(N)
  # Jm <- ncol(Y) - 1L
  # Ym <- Y[, 1L:Jm]
  # IJm <- I * Jm

  # Adapt beta to the model
  Marg <- model_building(N = N, Y = Y, Imod = Imod, start = start, cvar = ncol(C), cs = cs)
  # beta <- Marg[[5L]]
  if(nrow(Y) > 1L){
    Darg <- list(N = N, Y = Y[, 1L:(ncol(Y) - 1L)], C = C)
  } else {
    Darg <- list(N = N, Y = Y[, 1L:(ncol(Y) - 1L)], C = C)
    Darg[[2]] <- matrix(Darg[[2]], nrow = 1L, ncol = length(Darg[[2]]))
  }

  # Fisher scoring
  if (ncol(Darg[[3L]]) > 0) { # with covariates
    f.res <- fiscoVT(Darg = Darg, Marg = Marg, maxit = mit,
                     tol = tol, verbose = verbose, save.beta = save.beta)
  } else { # without covariates
    f.res <- fiscoxLC(Darg = Darg, Marg = Marg, maxit = mit,
                      tol = tol, verbose = verbose, save.beta = save.beta)
  }

  beta <- f.res$beta
  V <- f.res$V

  # Compute estimates
  res <- TraSeM(be = beta, V = V, Darg = Darg, Marg = Marg)
  V[V < 10^-20] <- 10^-20

  return(list("Q" = res$Q, "beta" = beta, "La" = res$La, "Vp" = res$Vp,
              "sbe" = sqrt(diag(V)), "madis" = res$madis, "lk" = f.res$lk, "iter" = f.res$iter,
              "V" = V))
}
##------------------------------------------------------------------------------

##------------------------------------------------------------------------------

## Function votlan_unit()##
votlan_unit <- function(N, Y, C, Imod, start, cs, seed, mit, tol, verbose, save.beta) {
  # N: results in election 1
  # Y: results ion election 2
  # C: matrix of covariates
  # Imod: desing matrix, output of Imod_function
  # start: matrix of initial estimates for the beta vector
  # mit: maximum number of iterations?

  #  -  MLE of vote transitions (Forcina, Marcehtti, CLDAG 2010)
  #     rows of N and Y are electoral results at time t-1 and t
  #     if model is constant across units, beta can be omitted and
  #     initial estimates of P and La provided; last col of P is
  #     used as reference; log-oversdipersion; C contains covariates
  #     or has 0 columns
  #  -  calls  ->  fiscor

  # Preliminaries
  # J <- ncol(Y)
  # I <- ncol(N)
  # Jm <- ncol(Y) - 1L
  # Ym <- Y[, 1L:Jm]
  # IJm <- I * Jm

  # Adapt beta to the model
  Marg <- model_building(N = N, Y = Y, Imod = Imod, start = start, cvar = ncol(C), cs = cs)
  # beta <- Marg[[5L]]
  if(nrow(Y) > 1L){
    Darg <- list(N = N, Y = Y[, 1L:(ncol(Y) - 1L)], C = C)
  } else {
    Darg <- list(N = N, Y = Y[, 1L:(ncol(Y) - 1L)], C = C)
    Darg[[2]] <- matrix(Darg[[2]], nrow = 1L, ncol = length(Darg[[2]]))
  }

  # Fisher scoring
  if (ncol(Darg[[3L]]) > 0) { # with covariates
    f.res <- fiscoVT(Darg = Darg, Marg = Marg, maxit = mit,
                     tol = tol, verbose = verbose, save.beta = save.beta)
  } else { # without covariates PENDIENTE
    f.res <- fiscoxLC(Darg = Darg, Marg = Marg, maxit = mit,
                      tol = tol, verbose = verbose, save.beta = save.beta)
  }

  beta <- f.res$beta
  V <- f.res$V

  # Compute estimates
  res <- TraSeM_unit(be = beta, V = V, Darg = Darg, Marg = Marg)
  res$Q[res$Q < 10^-15] <- 10^-15

  return(list("Q" = res$Q))
}
##------------------------------------------------------------------------------


## Function fiscoVT##
# Performs Fisher scoring when covariates are available
fiscoVT <- function(Darg, Marg, maxit = 100, tol = 0.0001, o = 5.8, verbose = TRUE, save.beta = TRUE) {
  # Darg: data arguments
  # Marg: matrix argument
  # flik: likelihood function
  # Finds mle by Fisher scoring, given data in Darg, retrieved from GM
  # modified version without steepest ascent
  # tol <- 0.0005  # For testing convergence
  # o <- 5.8  # Adjust step length

  # Preliminares
  # Initial variables
  g <- rep(0L, 5L)
  s <- 1L
  ico <- 1L
  beta <- Marg[[5L]]
  q <- length(beta)

  # Inicialization
  it <- 1L
  test <- TRUE
  # vl0 <- rep(Inf, 4L)

  # Iterate
  while (test & it <= maxit) {
    sw <- 3
    Marg[[5]] <- beta
    # First call
    if (ico > 0) {
      res1 <- liderlac(Darg = Darg, Marg = Marg, sw = sw)
      lk <- res1$lk
      d1 <- res1$d1
      D <- res1$d2
      d1r <- d1[abs(beta) < 12]
      tgra <- mean(abs(d1r))

      # check singularity
      eigs <- eigen(D)
      dL <- eigs$values
      mL <- min(dL)
      if (mL > -10^(-10)) {
        dL <- dL + 10^(-9)*rep(1L, length(beta))
      } else {
        stop(paste("Singular matrix encountered. The process stopped before converging.\n",
                   "If you use `save.beta = TRUE` when calling the function, you could restart the process\n",
                   "using as initial values the beta vector saved in `beta.RData`."))
      }
      Di <- eigs$vectors %*% diag(dL^(-1)) %*% t(eigs$vectors)
      dlt <- Di %*% d1
    } else {
      o <- o / 2
      dlt <- dlt + stats::rnorm(q) / 200
    } # steepest ascent

    # compute step direction
    delta <- o * dlt / (sqrt(sum(dlt^2)) + ico + 1)
    if (any(is.complex(delta))) {
      stop(paste("Complex numbers in delta. The process stopped before converging.\n",
                 "If you use `save.beta = TRUE` when calling the function, you could restart the process\n",
                 "using as initial values the beta vector saved in `beta.RData`."))
    }

    # loglik and derivative
    g[1] <- lk
    g[4] <- sum(d1 * dlt)

    # Update beta
    b <- beta + delta
    sw <- 2

    # Second call
    Marg[[5]] <- b
    res2 <- liderlac(Darg = Darg, Marg = Marg, sw = sw)
    d1 <- res2$d1
    d1r <- d1[abs(b) < 12]
    tgra <- min(tgra, mean(abs(d1r))) # test convergence

    # Log-likelihood and derivative
    g[2] <- res2$lk
    g[5] <- sum(d1 * dlt)

    # adjust step length
    if (g[4] > g[5] & g[5] > 0) {
      o <- o * (1 + 2 * (g[5] / g[4]))
    } else if (g[5] < 0 & abs(g[5]) > g[4]) {
      o <- o / 6
    }
    if (o < 0.01) {
      o <- 0.01
    }
    if (o > 8) {
      o <- 8
    }

    # Compute new step length
    s <- cubic(g, 1)
    if (s > 8) {
      s <- 8
    } else if (s < 0.01) {
      s <- 0.01
    }
    # Second guess and update
    sw <- 1
    b3 <- beta + s * delta
    Marg[[5]] <- b3
    res3 <- liderlac(Darg = Darg, Marg = Marg, sw = sw)
    g[3] <- res3$lk

    # check convergence
    if (tgra < tol) {
      test <- FALSE
    }

    # second guess the best
    if (g[3] >= g[1] & g[3] >= g[2]) {
      beta <- b3
      sw <- 3
      ico <- 3
      # first guess the best
    } else if (g[2] > g[3] & g[2] >= g[1]) {
      beta <- b
      sw <- 3
      ico <- 2
      # else steepest ascent or bisect
    } else if (g[1] > g[2] & g[1] > g[3]) {
      # last attempt
      s <- s / 3
      sw <- 1
      b4 <- beta + s * delta
      Marg[[5]] <- b4
      lk <- liderlac(Darg = Darg, Marg = Marg, sw = sw)$lk
      if (lk > g[1]) {
        g[3] <- lk
        beta <- b4
        ico <- 4
      } else {
        ico <- 1
      }
    }

    # if(max(abs(vl0 - c(s, o, ico, tgra))) < tol*10^-4) test <- FALSE
    # vl0 <- c(s, o, ico, tgra)

    # display LL and derivatives
    if(verbose){
      cat("%%%%%\n")
      cat(it, "\n")
      cat(-g[1:3] / 100, "\n")
      cat(g[4:5], "\n")
      cat(s, o, ico, tgra / 100, "\n")
    }
    it <- it + 1L

    if (save.beta) {
       save(beta, file = "beta.RData")
    }

  }  # End while

  # information matrix
  V <- Di

  # print results
  if(verbose){
    cat("%%%%%\n")
    cat(it, -g[1:3], "\n")
    cat(g[4:5], "\n")
    cat(s, o, tgra, "\n")
  }

  # Output
  return(list("beta" = beta, "V" = V, "lk" = g[1L], "iter" = it))
}
##------------------------------------------------------------------------------


## Function liderlac ##
liderlac <- function(Darg, Marg, sw = 3){
  # Marg is the 4 first components of the output of model_building
  # beta is the fifth component of the output of model_building
  # Darg: data N, Ym, C and ncol(C)
  # sw: determines what computes,
  #    1: just likelihood,
  #    2: likelihood and first derivatives,
  #    3: likelihood, first derivatives and matrix of expected values of second derivatives
  #    4: likelihood, first derivatives and empirical information matrix (faster)

  # Preliminaries
  X1 <- Marg[[1L]]
  X2 <- Marg[[2L]]
  csi <- Marg[[4L]]
  beta <- Marg[[5L]]
  N <- Darg[[1L]]
  Ym <- Darg[[2L]]
  C <- Darg[[3L]]
  I <- ncol(N)
  K <- nrow(N)
  Jm <- ncol(Ym)
  IJm <- I * Jm
  d1 <- d2 <- S <- NULL

  # Compute P and Theta
  Q_index_La <- linklac(Darg = Darg, Marg = Marg)
  Q <- Q_index_La$Q
  index <- Q_index_La$index
  La <- Q_index_La$lambda

  # Initialize
  lk <- 0
  sib <- length(beta)
  sid <- IJm + I

  if (sw > 1) {
    d1 <- rep(0L, sib)
  }

  if (sw > 2) {
    d2 <- matrix(0L, nrow = sib, ncol = sib)
  }

  # if (sw > 3) {
  #   S <- d2
  # }

  #if (ncol(C) == 0L) {
  #  C <- matrix(0L, nrow = K, ncol = 1L)
  #}

  # T <- '.-+*'
  U <- as.matrix(rep(1L, IJm + I))

  # Clicking time
  # QK <- ceiling(K * (1:60) / 60)

  # Cycle across units
  for (k in 1L:K) {

    n <- as.matrix(N[k, ])
    y <- as.matrix(Ym[k, ])
    csn <- csi * (n > 0) * (n - 1L) / (n + (n == 0))
    Ck <- U %*% C[k, ]
    X <- cbind(X1, Ck * X2)
    P <- Q[index[, k], ]
    nb <- as.matrix(n * (1L + La * csn))

    # Variance matrix
    Vy <- diag(as.vector(t(P) %*% nb)) - t(P) %*% diag(as.vector(nb)) %*% P
    e <- y - t(P) %*% n
    H <- solve(Vy)

    # Twice likelihood
    lk <- lk - (log(det(Vy)) + t(e) %*% H %*% e)

    # First derivatives
    if (sw > 1) {
      He <- H %*% e
      Hd <- H - He %*% t(He)
      dH <- diag(Hd)
      u <- as.matrix(rep(0L, sid))
      V <- matrix(0L, nrow = sid, ncol = sid)

      for (i in 1L:I) {
        ni <- n[i]
        nbi <- nb[i]
        pi <- P[i, ]
        Oi <- diag(pi) - pi %*% t(pi)
        Lai <- csn[i] * La[i] * (1 - La[i])
        # components wrt P
        rp <- ((i - 1L) * Jm + 1L):(i * Jm)
        u[rp] <- Oi %*% (nbi * (-dH / 2 + Hd %*% pi) + ni * He)
        u[IJm + i] <- -ni * Lai * sum(diag(Hd %*% Oi)) / 2

        # Second derivatives
        if (sw > 2) {
          H2 <- H^2 / 2
          hi <- H %*% pi
          for (l in 1L:I) {
            nl <- n[l]
            nbl <- nb[l]
            pl <- P[l, ]
            Lal <- csn[l] * La[l] * (1 - La[l])
            Ol <- diag(pl) - pl %*% t(pl)
            hl <- H %*% pl
            jn <- ((l - 1L) * Jm + 1L):(l * Jm)
            if (l <= i) {
              # wrt to beta_i, beta_l'
              A <- nbi * nbl * (H2 - diag(as.vector(hl))%*%H - t(diag(as.vector(hi))%*%H) +
                                  hl %*% t(hi) + H * as.vector(t(pl) %*% H %*% pi))
              A <- A + ni * nl * H
              V[rp, jn] <- Oi %*% A %*% Ol
              #  wrt to La_i, La_l
              a <- ni * nl * Lai * Lal / 2
              a <- a * sum(H %*% Ol %*% H %*% Oi)
              V[IJm + i, IJm + l] <- a
            }
            #  wrt to La_i,beta_l'
            a <- ni * nbl * Lai / 2
            B <- H %*% Oi %*% H
            V[IJm + i, jn] <- as.vector(a * (Ol %*% (diag(B) - 2 * B %*% pl)))
          } # end for d2_units
        } # end sw > 2
      } # end for d1_units
    } # end sw > 1

    # Cumulate
    if (sw > 1) {
      si <- t(X) %*% u
      d1 <- d1 + si
    }

    if (sw > 2) {
      V1 <- V2 <- V
      V1[!lower.tri(V, diag = TRUE)] <- 0
      V2[!lower.tri(V, diag = FALSE)] <- 0
      V <- V1 + t(V2)
      d2 <- d2 + t(X) %*% V %*% X
      d2 <- (d2 + t(d2)) / 2
    }

    # if (sw > 3) {
    #    S <- S + si %*% t(si)
    # }
  }

  lk <- lk / 2

  return(list("lk"= lk, "d1" = d1, "d2" = d2, "S"= S))

}
##------------------------------------------------------------------------------

## Function cubic ##
cubic <- function(g, o = 1, s = NULL, L = NULL) {
  # Computes step length by fitting a cubic polynomial
  # Preliminaries
  g <- g[-3]

  # Computation of parameters
  X <- matrix(c(1, 0, 0, 0,           # lk(0)
                1, o, o^2, o^3,       # lk(1)
                0, 1, 0, 0,           # d1(0)
                0, 1, 2*o, 3*o^2),    # d1(1)
              nrow = 4, byrow = TRUE)

  if (!is.null(s) & !is.null(L)) {
    A <- cbind(rep(1, nrow(s)), s, s^2, s^3)
    g <- c(g, L)
    X <- rbind(X, A)
  }

  be <- solve(t(X) %*% X) %*% (t(X) %*% as.matrix(g))

  # Finding the maximum
  c <- be[2]
  b <- 2 * be[3]
  a <- 3 * be[4]

  dis <- b^2 - 4 * a * c

  if (dis > 0) {
    dis <- sqrt(dis)
    xn <- (-b - dis) / (2 * a)
    xp <- (-b + dis) / (2 * a)
    Hn <- b + 2 * a * xn
    Hp <- b + 2 * a * xp

    if (is.na(xn) | is.na(Hn)){
      s <- 0.1
    } else if (xn > 0 & Hn < 0) {
      s <- xn
    } else if (xp > 0 & Hp < 0) {
      s <- xp
    } else {
      s <- 0.1
    }
  } else {
    s <- -1 # steepest ascent
  }

  return(s)
}
##------------------------------------------------------------------------------

## Function linklac ##
linklac <- function(Darg, Marg){
  # This function is used by liderlac

  # Preliminaries
  X1 <- Marg[[1L]]
  X2 <- Marg[[2L]]
  off <- Marg[[3L]]
  beta <- Marg[[5L]]
  C <- as.matrix(Darg[[3L]])
  I <- ncol(Darg[[1L]])
  Jm <- ncol(Darg[[2L]])
  K <- nrow(Darg[[3L]])
  IJm <- I * Jm

  # Initialize variables
  u <- as.matrix(rep(1L, IJm + I))
  P <- matrix(0L, nrow = I * K, ncol = Jm)
  index <- matrix(1L:(I * K), nrow = I)

  # Cycle across units
  for (k in 1L:K) {
    Ck <- u %*% C[k, ]
    X <- cbind(X1, Ck * X2)
    ett <- X %*% beta + off
    ett <- ett * (abs(ett) < 680) + 680 * sign(ett) * (abs(ett) >= 680)
    et <- matrix(exp(ett[1L:IJm]), ncol = Jm, byrow = TRUE)
    # et <- t(matrix(exp(ett[1L:IJm]), ncol = Jm, byrow = FALSE))

    # Compute P
    P1 <- diag(1L / (1L + rowSums(et))) %*% et
    P[index[, k], ] <- P1
  }

  # Compute Lambda
  etl <- ett[(IJm + 1L):(IJm + I)]
  etl <- etl * (etl < 12) + 12 * (etl >= 12)
  La <- exp(etl) / (1L + exp(etl))

  return(list("Q"= P, "index" = index, "lambda" = La))
}
##------------------------------------------------------------------------------

## Function TraSeM()##
TraSeM <- function(be, V, Darg, Marg) {
  # Preliminaries
  X1 <- Marg[[1L]]
  X2 <- Marg[[2L]]
  off <- Marg[[3L]]
  csi <- Marg[[4L]]
  N <- Darg[[1L]]
  Ym <- Darg[[2L]]
  C <- Darg[[3L]]
  K <- nrow(N)
  I <- ncol(N)
  Jm <- ncol(Ym)
  sT <- I * Jm
  IJm <- 1L:sT
  nt <- t(as.matrix(colSums(N)))
  R <- t(diag(1/as.vector(nt))%*% t(N))
  Im <- t(matrix(IJm, nrow = Jm))
  Q <- rep(0L, I * Jm)
  madis <- rep(0, K)
  be <- be * (abs(be) < 640) + 640 * sign(be) * (abs(be) >= 640)

  # Betas for transitions
  X <- cbind(X1, X2)
  isx <- which(colSums(X[IJm, ]) > 0)
  btra <- be[isx]
  Om <- matrix(0, nrow = sT, ncol = sT)
  sx <- length(btra)
  D <- matrix(0, nrow = sT, ncol = sx)
  V <- V[isx, isx]

  # Compute Lambda
  et <- X %*% be + off
  etL <- et[(sT + 1L):(sT + I)]
  La <- exp(etL) / (1 + exp(etL))

  # Compute P and Se
  if (ncol(C) == 0) { # without covariates
    X <- X1[IJm, isx]
    et <- X %*% btra + off[IJm]
    Et <- t(matrix(exp(et), nrow = Jm, byrow = FALSE))
    P <- diag(1/(1 + rowSums(Et))) %*% Et
    for (i in 1:I) {
      pik <- P[i, ]
      im <- Im[i, ]
      Om[im, im] <- diag(pik) - pik %*% t(pik) # Omega(pik)
    }
    for (k in 1:K) {
      # Mahalanobis
      n <- as.matrix(N[k, ])
      y <- as.matrix(Ym[k, ])
      csn <- csi * (n > 1) * (n - 1) / (n + (n == 0))
      nb <- n * (1 + La * csn)
      A <- kronecker(as.matrix(nt), diag(Jm)) %*% Om %*% X
      Vy <- diag(as.vector(t(P) %*% nb)) - t(P) %*% diag(as.vector(nb)) %*% P + A %*% V %*% t(A)
      e <- y - t(P) %*% n
      madis[k] <- t(e) %*% solve(Vy) %*% e
    }
    D <- Om %*% X
    Q <- P

  } else { # with covariates

    for (k in 1:K) {
      # Compute P
      ck <- C[k, ]
      if (ncol(C) == 1L) ck <- as.matrix(ck)
      X <- cbind(X1, X2 %*% diag(ck))[IJm, isx]
      et <- X %*% btra + off[IJm]
      Et <- t(matrix(exp(et), nrow = Jm, byrow = FALSE))
      Pk <- diag(1 / (1 + rowSums(Et))) %*% Et

      # Variances of transitions
      for (i in 1:I) {
        pik <- Pk[i, ]
        im <- Im[i, ]
        Om[im, im] <- Omega(pik)
      } # end for i

      A <- kronecker(diag(R[k, ]), diag(Jm))
      Q <- Q + A %*% as.vector(t(Pk))
      A <- A %*% Om %*% X
      D <- D + A

      # Mahalanobis
      n <- as.matrix(N[k, ])
      y <- as.matrix(Ym[k, ])
      csn <- csi * (n > 1) * (n - 1) / (n + (n == 0))
      nb <- n * (1 + La * csn)
      A <- kronecker(as.matrix(nt), diag(Jm)) %*% Om %*% X
      Vy <- diag(as.vector(t(Pk) %*% nb)) - t(Pk) %*% diag(as.vector(nb)) %*% Pk + A %*% V %*% t(A)
      e <- y - t(Pk) %*% n
      madis[k] <- t(e) %*% solve(Vy) %*% e
    } # end for k
    Q <- t(matrix(Q, nrow = Jm, byrow = FALSE))
  }

  # Se transitions
  Vp <- D %*% V %*% t(D)
  Q <- cbind(Q, 1 - rowSums(Q))
  Q[Q < 10^-15] <- 10^-15
  TT <- kronecker(diag(I), rbind(diag(Jm), rep(1, Jm)))
  Vp <- TT %*% Vp %*% t(TT)
  Vp[Vp < 10^-10] <- 10^-10
  Vp <- sqrt(t(matrix(diag(Vp), nrow = Jm + 1L, byrow = FALSE)))
  return(list("Q" = Q, "Vp" = Vp, "La" = La, "madis" = madis))
}
##------------------------------------------------------------------------------

## Function TraSeM_unit()##
TraSeM_unit <- function(be, V, Darg, Marg) {
  # Preliminaries
  X1 <- Marg[[1L]]
  X2 <- Marg[[2L]]
  off <- Marg[[3L]]
  csi <- Marg[[4L]]
  N <- Darg[[1L]]
  Ym <- Darg[[2L]]
  C <- Darg[[3L]]
  K <- nrow(N)
  I <- ncol(N)
  Jm <- ncol(Ym)
  sT <- I * Jm
  IJm <- 1L:sT
  nt <- t(as.matrix(colSums(N)))
  R <- t(diag(1/as.vector(nt))%*% t(N))
  Im <- t(matrix(IJm, nrow = Jm))
  Q <- rep(0L, I * Jm)
  be <- be * (abs(be) < 640) + 640 * sign(be) * (abs(be) >= 640)

  # Betas for transitions
  X <- cbind(X1, X2)
  isx <- which(colSums(X[IJm, ]) > 0)
  btra <- be[isx]
  Om <- matrix(0, nrow = sT, ncol = sT)
  sx <- length(btra)
  D <- matrix(0, nrow = sT, ncol = sx)
  V <- V[isx, isx]

  # Compute Lambda
  et <- X %*% be + off
  etL <- et[(sT + 1L):(sT + I)]
  La <- exp(etL) / (1 + exp(etL))

  # Compute P and Se
  if (ncol(C) == 0) { # without covariates
    X <- X1[IJm, isx]
    et <- X %*% btra + off[IJm]
    Et <- t(matrix(exp(et), nrow = Jm, byrow = FALSE))
    P <- diag(1/(1 + rowSums(Et))) %*% Et
    for (i in 1:I) {
      pik <- P[i, ]
      im <- Im[i, ]
      Om[im, im] <- diag(pik) - pik %*% t(pik) # Omega(pik)
    }
    for (k in 1:K) {
      # Mahalanobis
      n <- as.matrix(N[k, ])
      y <- as.matrix(Ym[k, ])
      csn <- csi * (n > 1) * (n - 1) / (n + (n == 0))
      nb <- n * (1 + La * csn)
      A <- kronecker(as.matrix(nt), diag(Jm)) %*% Om %*% X
      Vy <- diag(as.vector(t(P) %*% nb)) - t(P) %*% diag(as.vector(nb)) %*% P + A %*% V %*% t(A)
      e <- y - t(P) %*% n
    }
    D <- Om %*% X
    Q <- P

  } else { # with covariates

    for (k in 1:K) {
      # Compute P
      ck <- C[k, ]
      if (ncol(C) == 1L) ck <- as.matrix(ck)
      X <- cbind(X1, X2 %*% diag(ck))[IJm, isx]
      et <- X %*% btra + off[IJm]
      Et <- t(matrix(exp(et), nrow = Jm, byrow = FALSE))
      Pk <- diag(1 / (1 + rowSums(Et))) %*% Et

      # Variances of transitions
      for (i in 1:I) {
        pik <- Pk[i, ]
        im <- Im[i, ]
        Om[im, im] <- Omega(pik)
      } # end for i

      A <- kronecker(diag(R[k, ]), diag(Jm))
      Q <- Q + A %*% as.vector(t(Pk))
      A <- A %*% Om %*% X
      D <- D + A

      # Mahalanobis
      n <- as.matrix(N[k, ])
      y <- as.matrix(Ym[k, ])
      csn <- csi * (n > 1) * (n - 1) / (n + (n == 0))
      nb <- n * (1 + La * csn)
      A <- kronecker(as.matrix(nt), diag(Jm)) %*% Om %*% X
      Vy <- diag(as.vector(t(Pk) %*% nb)) - t(Pk) %*% diag(as.vector(nb)) %*% Pk + A %*% V %*% t(A)
      e <- y - t(Pk) %*% n
    } # end for k
    Q <- t(matrix(Q, nrow = Jm, byrow = FALSE))
  }

  # Se transitions
  Q <- cbind(Q, 1 - rowSums(Q))
  return(list("Q" = Q))
}
##------------------------------------------------------------------------------


## Function Omega##
Omega <- function(x){
  x <- as.matrix(x)
  size <- dim(x)
  if (size[1L] < size[2L]) x <- t(x)
  y <- diag(as.vector(x)) - x %*% t(x)
  return(y)
}
##------------------------------------------------------------------------------

## Function fiscoxLC ##
fiscoxLC <- function(Darg, Marg, maxit = 100, tol = 0.0002, o = 5, verbose = TRUE, save.beta = FALSE) {
  # Darg: data arguments
  # Marg: matrix argument
  # flik: likelihood function
  # Finds mle by Fisher scoring,
  # Finds mle by Fisher scoring, given data in Darg
  #  Output
  #       s: step length
  #       g(1,2,3): log-lik at the 3 basic steps
  #       g(4,5): derivative of g(1,2)
  #       a(1,2): coefficients of quadratic or linear approximation
  #       a(3): discriminant in quadratic approximation


  # Preliminaries
  # Initial variables
  g <- a <- rep(0L, 5L)
  s <- 0.1
  beta <- Marg[[5L]]
  q <- length(beta)
  sain <- 0

  # Initialize
  it <- 1
  test <- TRUE
  q <- length(beta)

  # Iterate
  while (test & it <= maxit) {
    bind <- 0
    sw <- 3
    Marg[[5]] <- beta
    res1 <- liderla(Darg = Darg, Marg = Marg,  sw = sw)
    lk <- res1$lk
    d1 <- res1$d1
    D <- res1$d2
    tgra <- max(abs(d1))

    # check singularity
    if (sain == 0) {
      eigs <- eigen(D)
      dL <- eigs$values
      mL <- min(dL)
      if (mL < 10^(-6)) {
        dL <- dL + 10^(-6)*rep(1L, length(beta))
        if (mL < 10^(-5) & verbose) print(mL)
      }
      Di <- eigs$vectors %*% diag(dL^(-1)) %*% t(eigs$vectors)
      dlt <- Di %*% d1
    }

    if (sain > 0) { # steepest ascent
      dlt <- d1
    }

    # compute step direction
    delta <- o * dlt / (sqrt(sum(dlt^2)) + sain + 1)

    # Reduce step length if last cycle went too far
    if (s < 0.1) {
      delta <- s * delta
    }

    if (any(!is.na(delta) & !is.complex(delta) == 0)) {
      stop(paste("Step direction contains complex numbers. The process stopped before converging.\n",
                 "If you use `save.beta = TRUE` when calling the function, you could restart the process\n",
                 "using as initial values the beta vector saved in `beta.RData`."))
    }

    # Loglik and der
    g[1] <- 2*lk
    g[4] <- sum(d1 * dlt)

    # Update beta
    b <- beta + delta
    sw <- 2

    # Second call
    Marg[[5]] <- b
    res2 <- liderla(Darg = Darg, Marg = Marg,  sw = sw)
    lk <- res2$lk
    d1 <- res2$d1

    tgra <- min(tgra, max(abs(d1))) #  test convergence

    # Loglik and der
    g[2] <- 2*lk
    g[5] <- sum(d1 * dlt)

    # adjust step length
    if (g[4] > g[5] & g[5] > 0) { #  step too short
      o <- o * (1 + 10 * (g[5] / g[4]) ^ 0.5)
    } else if (g[5] < 0 & abs(g[5]) > g[4]) { # step too long
      o <- o / 1.6
    }
    if (o < 0.01) {
      o <- 0.01
    }
    if (o > 10) {
      o <- 10
    }

    # Compute new step length
    a[1] <- 3 * (g[4] + g[5]) - 6 * (g[2] - g[1])
    a[4] <- g[5] - g[4]
    a[2] <- a[4] - a[1]
    a[3] <- a[2]^2 - 4 * a[1] * g[4]
    a[5] <- abs(a[1])

    s <- fit_ag(a, g)

    if (s > 2) s <- 2
    if (s < 0.1) s <- 0.1

    # Second guess and update
    rep <- TRUE

    # Repeated bisection
    while (rep) {
      rep <- FALSE
      sw <- 1
      b3 <- beta + s * delta
      Marg[[5]] <- b3
      res3 <- liderla(Darg = Darg, Marg = Marg,  sw = sw)
      lk <- 2*res3$lk

      g[3] <- lk

      # display
      if(verbose){
        cat("%%%%%\n")
        cat(it, "\n")
        cat(-g[1:3] / 100, "\n")
        cat(g[4:5], "\n")
        cat(s, o, tgra / 1000, "\n")
      }

      # Check convergence
      if (tgra < tol) {
        test <- FALSE
      }

      # Second guess best
      if (g[3] >= g[1] & g[3] >= g[2]) {
        beta <- b3
        sain <- 0
      }
      # First guess best
      if (g[2] > g[3] & g[2] >= g[1]) {
        beta <- b
        sain <- 0
      }
      # else steepest ascent or bisect
      if (g[1] > g[2] & g[1] > g[3]) {
        if (bind < 1) {
          delta <- delta / 10
          bind <- bind + 1
          sain <- 0
          rep <- TRUE
          if(verbose) cat("Bisection.\n")
        } else {
          sw <- 2
          sain <- sain + 1
          if(verbose) cat("Steepest ascent.\n")
        }
      }
    } # end bisection

    # Save results
    if (any(is.na(beta))) {
      stop(paste("NA produced. The process stopped before converging.\n",
                 "If you use `save.beta = TRUE` when calling the function, you could restart the process\n",
                 "using as initial values the beta vector saved in `beta.RData`."))
    } else if (save.beta) {
       save(beta, file = "beta.RData")
    }

    it <- it + 1L
  } # end while

  # Information matrix
  V <- Di

  if(verbose){
    # Print results
    cat("%%%%%\n")
    cat(it, "\n")
    cat(-g[1:3] / 100, "\n")
    cat(g[4:5], "\n")
    cat(s, o, tgra / 1000, "\n")
  }

  return(list("beta" = beta, "V" = V, "lk" = lk/2, "iter" = it))
}
##------------------------------------------------------------------------------

## Function fit_ag ##
fit_ag <- function(a, g) {
  if (a[3] > 0 & a[5] > 0.1) {
    s <- (-a[2] - sqrt(a[3])) / (2 * a[1])

    if (g[5] < 0) {
      s <- (s - 2 * g[4] / a[4]) / 3
    }
  }

  if (a[3] <= 0 | a[5] <= 0.1) {
    if (a[4] < 0) {
      s <- -g[4] / a[4]
    }

    if (a[4] >= 0) {
      s <- 0.5
    }
  } else {
    s <- 0.1
  }

  return(s)
}
##------------------------------------------------------------------------------

## Function  liderla ##
liderla <- function(Darg, Marg, sw) {
  # Preliminaries
  X <- Marg[[1]]
  cs <- Marg[[4]]
  beta <- Marg[[5]]
  N <- Darg[[1]]
  Ym <- Darg[[2]]
  C <- Darg[[3]]
  I <- ncol(N)
  K <- nrow(N)
  Jm <- ncol(Ym)
  IJm <- I * Jm
  d1 <- d2 <- NULL

  # Compute P and Lambda
  res <- linkla(Darg = Darg, Marg = Marg)
  P <- res$P
  La <- res$La

  # Initialize
  lk <- 0
  sid <- IJm + I

  if (sw > 1) {
    d1 <- matrix(0, nrow = sid, ncol = 1)
  }

  if (sw > 2) {
    d2 <- matrix(0, nrow = sid, ncol = sid)
  }

  # Cycle across units
  # QK <- ceiling(K * (1:60) / 60)

  for (k in 1:K) {
    n <- as.matrix(N[k, ])
    y <- as.matrix(Ym[k, ])
    csn <- cs * (n > 0) * (n - 1L) / (n + (n == 0))
    nb <- as.matrix(n * (1L + La * csn))

    # Variance matrix
    Vy <- diag(as.vector(t(P) %*% nb)) - t(P) %*% diag(as.vector(nb)) %*% P
    e <- y - t(P) %*% n
    H <- solve(Vy)

    # Twice likelihood
    lk <- lk - (log(det(Vy)) + t(e) %*% H %*% e)

    # First derivatives
    if (sw > 1) {
      He <- H %*% e
      Hd <- H - He %*% t(He)
      dH <- diag(Hd)
      u <- as.matrix(rep(0L, sid))
      V <- matrix(0L, nrow = sid, ncol = sid)

      for (i in 1:I) {
        ni <- n[i]
        nbi <- nb[i]
        pi <- P[i, ]
        Oi <- diag(pi) - pi %*% t(pi)
        Lai <- csn[i] * La[i] * (1 - La[i])
        #  components wrt P
        rp <- ((i - 1L) * Jm + 1L):(i * Jm)
        u[rp] <- Oi %*% (nbi * (-dH / 2 + Hd %*% pi) + ni * He)
        u[IJm + i] <- -ni * Lai * sum(diag(Hd %*% Oi)) / 2

        # Second derivatives
        if (sw > 2) {
          H2 <- H * H / 2
          hi <- H %*% pi

          for (l in 1:I) {
            nl <- n[l]
            nbl <- nb[l]
            pl <- P[l, ]
            Lal <- csn[l] * La[l] * (1 - La[l])
            Ol <- diag(pl) - pl %*% t(pl)
            hl <- H %*% pl
            jn <- ((l - 1L) * Jm + 1L):(l * Jm)

            if (l <= i) {
              #  wrt to beta_i, beta_l'
              A <- (H2 - diag(as.vector(hl))%*%H - t(diag(as.vector(hi))%*%H)) +
                hl %*% t(hi) + H * as.vector(t(pl) %*% H %*% pi)
              A <- nbi * nbl * A + ni * nl * H
              V[rp, jn] <- Oi %*% A %*% Ol
              # wrt to La_i, La_l
              a <- ni * nl * Lai * Lal / 2
              a <- a * sum(H %*% Ol %*% H %*% Oi)
              V[IJm + i, IJm + l] <- a
            }

            a <- ni * nbl * Lai / 2
            B <- H %*% Oi %*% H
            V[IJm + i, jn] <- as.vector(a * (Ol %*% (diag(B) - 2 * B %*% pl)))
          }
        }
      }
    }

    # Cumulate
    if (sw > 1) {
      d1 <- d1 + u
    }

    if (sw > 2) {
      d2 <- d2 + V
    }
  }

  lk <- lk/2

  if (sw > 1) {
    d1 <- t(X) %*% d1
  }

  if (sw > 2) {
    V1 <- V2 <- d2
    V1[!lower.tri(V, diag = TRUE)] <- 0
    V2[!lower.tri(V, diag = FALSE)] <- 0
    V <- V1 + t(V2)
    d2 <- t(X) %*% V %*% X
    d2 <- (d2 + t(d2)) / 2
  }

  return(list("lk" = lk, "d1" = d1, "d2" = d2))
}
##------------------------------------------------------------------------------

## Function linkla ##
linkla <- function(Darg, Marg, eps = 1e-11) {
  # This function is used by liderla
  # Preliminaries
  X1 <- Marg[[1L]]
  off <- Marg[[3L]]
  b <- Marg[[5L]]
  I <- ncol(Darg[[1L]])
  Jm <- ncol(Darg[[2L]])
  IJm <- I * Jm

  et <- X1 %*% b + off
  etl <- et[(IJm + 1L):(IJm + I)]
  et <- et[1L:IJm]
  et <- matrix(exp(et), ncol = Jm, byrow = TRUE)

  # Compute P
  o <- matrix(1, nrow = 1, ncol = Jm)
  P <- et / (1 + rowSums(et) %*% o)
  P[P < eps] <- eps
  p <- rowSums(P)
  p[p > 1 - eps] <- 1 - eps
  p <- p / rowSums(P)
  P <- diag(p) %*% P

  # Compute Lambda
  La <- exp(etl) / (1 + exp(etl))

  return(list(P = P, La = La))
}
##------------------------------------------------------------------------------

# Function IPF##
IPF <- function(matriz, vector.columna, vector.fila, precision = 0.000001){
  nc <- length(vector.columna)
  nf <- length(vector.fila)
  vector.fila1 <- rowSums(matriz)
  R1 <- diag(as.vector(vector.fila)) %*% (diag(as.vector(1/vector.fila1)))
  R1[is.nan(R1)] <- 0L
  R1[is.infinite(R1)] <- 1L
  # R1 <- diag(as.vector(vector.fila)) %*% MASS::ginv(diag(as.vector(vector.fila1)))
  X1 <- R1 %*% matriz
  vector.columna1 <- colSums(X1)
  S1 <- diag(as.vector(vector.columna)) * (diag(as.vector(1/vector.columna1)))
  S1[is.nan(S1)] <- 0L
  S1[is.infinite(S1)] <- 1L
  # S1 <-  diag(as.vector(vector.columna)) * ginv(diag(as.vector(vector.columna1)))
  X2 <- X1 %*% S1
  dif <- max(abs(X2 - matriz))

  while (dif > precision){
    matriz <- X2
    vector.fila1 <- rowSums(matriz)
    R1 <- diag(as.vector(vector.fila)) %*% (diag(as.vector(1/vector.fila1)))
    R1[is.nan(R1)] <- 0L
    R1[is.infinite(R1)] <- 1L
    #   R1 <- diag(as.vector(vector.fila)) %*% ginv(diag(as.vector(vector.fila1)))
    X1 <- R1 %*% matriz
    vector.columna1 <- colSums(X1)
    S1 <- diag(as.vector(vector.columna)) * (diag(as.vector(1/vector.columna1)))
    S1[is.nan(S1)] <- 0L
    S1[is.infinite(S1)] <- 1L
    #   S1 <-  diag(as.vector(vector.columna)) * ginv(diag(as.vector(vector.columna1)))
    X2 <- X1 %*% S1
    dif <- max((abs(X2 - matriz)))
  }
  X2 <- X2*(vector.fila/rowSums(X2))
  X2[is.nan(X2)] <- 0L
  X2[is.infinite(X2)] <- 1L
  X2 <- t(t(X2)*(vector.columna/colSums(X2)))
  X2[is.nan(X2)] <- 0L
  X2[is.infinite(X2)] <- 1L

  return(X2)
}
##------------------------------------------------------------------------------

## Function local_units ##
## Function to adjust local margins using IPF
local_units <- function(X, Y, TM, tol){
  I <- ncol(X)
  J <- ncol(Y)
  K <- nrow(X)
  TR.units <- TR.votes.units <- array(NA, c(I, J, K))
  for (ii in 1L:K){
    # Adjustment of initial underlying probabilities
    TR.votes.units[, , ii] <- IPF(TM, Y[ii, ], X[ii, ], precision = tol)
    TR.votes.units[is.nan(TR.votes.units)] <- 0L
    TR.units[, , ii] <- TR.votes.units[, , ii]/rowSums(TR.votes.units[, , ii])
    TR.units[X[ii, ] == 0L, , ii] <- 0L
  }
  # Composition/aggregation matrix
  TR.votes <- apply(TR.votes.units, c(1L, 2L), sum)
  TR <- TR.votes/rowSums(TR.votes)

  return(list("TR" = TR, "TR.votes" = TR.votes, "TR.units" = TR.units,
              "TR.votes.units" = TR.votes.units)
  )
}
##------------------------------------------------------------------------------

## Function unit_lk ##
## Function to estimate each unit table using the likelihood assumed by the model
## employing as initial values the global estimates

unit_lk <- function(N,
                    Y,
                    start,
                    dispersion.rows,
                    cs = 50,
                    seed = NULL,
                    max.iter = 100,
                    tol = 0.0001
){

  I <- ncol(N) # Number of rows of the transfer matrix
  J <- ncol(Y) # Number of columns of the transfer matrix
  K <- nrow(N) # Number of polling units
  Imod.mat <- Imod_matrix(I = I, covariates = NULL,
                          null.cells = NULL,
                          null.cells.C = NULL,
                          row.cells.relationships = NULL,
                          row.cells.relationships.C = NULL,
                          pair.cells.relationships = NULL,
                          cells.fixed.logit = NULL,
                          dispersion.rows = dispersion.rows,
                          dispersion.row = NULL,
                          dispersion.fixed = NULL)

  TR.units <- TR.votes.units <- array(0L, c(I, J, K))

   for (kk in 1L:K){
      vector.row <- N[kk, , drop = FALSE]
      vector.col <- Y[kk, , drop = FALSE]

      TR.units[, , kk] <- votlan_unit(N = vector.row, Y = vector.col,
                                      C = matrix(0, nrow = K, ncol = 0),
                                      Imod = Imod.mat, start = start, cs = cs,
                                      seed = seed, mit = max.iter, tol = tol,
                                      verbose = FALSE, save.beta = FALSE)$Q

      # suma <- sum(abs(vector.row%*%TR.units[, , kk] - vector.col))
      TR.votes.units[, , kk] <- TR.units[, , kk]*as.vector(vector.row)
      TR.votes.units[, , kk] <- IPF(TR.votes.units[, , kk], Y[kk, ],
                                    N[kk, ], precision = tol)
      TR.votes.units[is.nan(TR.votes.units)] <- 0L
      TR.units[, , kk] <- TR.votes.units[, , kk]/rowSums(TR.votes.units[, , kk])
      TR.units[N[kk, ] == 0L, , kk] <- 0L
    }

    TR.votes <- apply(TR.votes.units, c(1L, 2L), sum)
    TR <- TR.votes/rowSums(TR.votes)

  return(list("TR" = TR, "TR.votes" = TR.votes, "TR.units" = TR.units,
              "TR.votes.units" = TR.votes.units))
}
##------------------------------------------------------------------------------


