#' Simulate RxC Tables from Overdispersed-Multinomial Models
#'
#' @description  Generates at random a set of RxC tables with the joint distribution of voters in two elections according to the model proposed in Forcina et al. (2012), as extension of Brown and Payne (1986), under the assumption that local units are homogeneous (no covariates). Results in the first election may be provided by the user or generated at random according to the overdispersed multinomial model.
#'
#' @author Antonio Forcina, \email{forcinarosara@@gmail.com}
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#'
#' @references Brown, P. and Payne, C. (1986). Aggregate data, ecological regression and voting transitions. *Journal of the American Statistical Association*, 81, 453–460. \doi{10.1080/01621459.1986.10478290}
#' @references Forcina, A., Gnaldi, M. and Bracalente, B. (2012). A revised Brown and Payne model of voting behaviour applied to the 2009 elections in Italy. *Statistical Methods & Applications*, 21, 109–119. \doi{10.1007/s10260-011-0184-x}
#'
#' @param n.units Either a positive integer number, `K`, indicating the number of polling units to be simulated, or
#'                a `KxR` data.frame of a matrix with the number of votes gained in election 1
#'                for each of the `R` options in each of the `K` units. If `n.units` is a matrix (data.frame) of
#'                counts (votes) the values of arguments `prop1` and `theta1` are ommitted.
#'
#' @param TM A row-standardized RxC matrix with the underlying global transition
#'           probabilities of the simulated elections. If the matrix is not row-standardized,
#'           it is internally row-standardized by the function.
#'
#' @param prop1 A vector of length R with the initial assumed probabilities of voting (to be simulated)
#'              for each of the R competing options in the first election. If the provided vector
#'              is not a set of probabilities (i.e., a vector of positive numbers adding to 1),
#'              it is internally standardized by the function.
#'
#' @param polling.sizes Either a vector of two components with two positive integer
#'                      numbers indicating the minimum and maximum number of voters
#'                      for each unit or a vector of length `n.units` of positive integer
#'                      numbers informing about the number of voters in each unit. When
#'                      `polling.sizes` is a vector of length two, a number of voters is
#'                      randomly assigned for each unit using a uniform distribution with
#'                      parameters the minimum and maximum values included in `polling.sizes`.
#'
#' @param theta1 A number between 0 and 1 used as the overdispersion parameter.
#'               This parameter is employed by the underlying Dirichlet distribution,
#'               in conjunction with `prop1`, to randomly generate vectors of probabilities
#'               for each unit. These vectors are then used to simulate the results
#'               of the first election. The smaller the value of this parameter,
#'               the closer the unit-level marginal distributions for the first election
#'               are to `prop1`. Default, `0.1`.
#'
#' @param theta2 Either a single number between 0 and 1 or a vector of length `nrow(TM)`
#'               containing numbers between 0 and 1. The values in `theta2` serve as
#'               overdispersion parameters and are used alongside the row-probability
#'               vectors in `TM` within the underlying Dirichlet distributions.
#'               These distributions are employed to generate probability vectors
#'               for each combination of unit, cluster, and row, which are then used
#'               to simulate vote transfers from the first to the second election.
#'               If `theta2` is a vector, each row is assigned a distinct overdispersion
#'               parameter based on its corresponding value. Default, `0.1`.
#'
#' @param cs A positive number indicating the average number of cluster size. Default, `50`.
#'
#' @param noise Either a single number between 0 and 1 or a vector of length `nrow(TM)`
#'              containing numbers between 0 and 1. These numbers account for the
#'              proportion of causal voters of each origin party (row).
#'              These numbers are used to introduce more variability, compared
#'              to the BPF model, into the simulations.
#'              If `noise > 0`, a 100*noise percentage of votes of each row of each unit
#'              are randomly assigned among the column parties.
#'              Default, `0`.
#'
#' @param simplify A TRUE/FALSE argument indicating whether the simulated RxCxK array of
#'                 counts by polling unit should be rearranged as a matrix of order Kx(RC).
#'                 Default, FALSE.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.

#'
#' @return
#' A list with the following components
#'
#'  \item{votes1}{ A matrix of order KxR with the results simulated in each polling unit for the first election.}
#'  \item{votes2}{ A matrix of order KxC with the results simulated in each polling unit for the second election..}
#'  \item{TM.global}{ A matrix of order RxC with the actual simulated global transfer matrix of counts.}
#'  \item{TM.units}{ An array of order RxCxK with the simulated transfer matrices of votes by polling unit. If
#'                 `simplify = TRUE` the simulated transfer matrices of votes are returned organized in a Kx(RC) matrix.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'
#' @export
#'
#' @importFrom stats rmultinom rgamma
#'
#' @family simulators for ecological inference overdispersed-multinomial models
#'
#' @examples
#' TMg <- matrix(c(0.6, 0.1, 0.3, 0.1, 0.7, 0.2, 0.1, 0.1, 0.8),
#'              byrow = TRUE, nrow = 3)
#' example <- simula_BPF(n.units = 100, TM = TMg, prop1 = c(0.3, 0.3, 0.4),
#'                       polling.sizes = c(750, 850))


simula_BPF <- function(n.units,
                       TM,
                       prop1,
                       polling.sizes,
                       theta1 = 0.1,
                       theta2 = 0.1,
                       cs = 50,
                       noise = 0,
                       simplify = FALSE,
                       ...
){

  inputs <- c(as.list(environment()), list(...))
  arggs <- tests_inputs_simula_BF(inputs)
  TM <- arggs$TM
  prop1 <- arggs$prop1
  sizes <- arggs$sizes

  if(is.null(dim(n.units))){
    K <- n.units
    votes1 <- matrix(NA, nrow = K, ncol = nrow(TM))
    for (kk in 1L:K){
      votes1[kk, ] <- rcmult(sizes[kk], prop1, theta1)
    }
  } else {
    votes1 <- as.matrix(n.units)
    K <- nrow(votes1)
  }

  names1 <- paste0("R", 1L:nrow(TM))
  names2 <- paste0("C", 1L:ncol(TM))
  names3 <- paste0("unit", 1L:K)


  # Results in second election
  output0 <- SimElection2(P = TM, th2 = theta2, cs = cs,
                          X = votes1, noise = noise)

#  if (noise > 0){
#    output0 <- add_noise(init = output0, noise = noise, prop1 = prop1)
#    for (kk in 1L:K){
#      votes1[kk, ] <- rowSums(output0$TM.u[, , kk])
#    }
#  }

  votes2 <- output0$Y
  TM.units <- output0$TM.u

  # Improving outputs presentation
  rownames(votes1) <- rownames(votes2) <- rownames(output0$TM.u2) <- names3
  colnames(votes1) <- names1
  colnames(votes2) <- names2
  dimnames(TM.units) <- list(names1, names2, names3)

  TM.global <- apply(TM.units, c(1L, 2L), sum)

  if(simplify){
    TM.units <- output0$TM.u2
    names4 <- paste0(rep(names1, each = nrow(TM)), rep(names2, ncol(TM)))
    colnames(TM.units) <- names4
  }

  output <- list("votes1" = votes1, "votes2" = votes2,
                 "TM.global" = TM.global, "TM.units" = TM.units,
                 "inputs" = inputs)

  class(output) <- "simula_BPF"
  return(output)
}


######################################
# Auxiliary functions

### tests_inputs_simula_BF
tests_inputs_simula_BF <- function(arggs){

  if(is.null(dim(arggs$n.units))){
    if(!(arggs$n.units > 0 & (arggs$n.units - round(arggs$n.units) == 0)))
       stop("Argument 'n.units' must be either a positive integer or a matrix of non-negative integers.")
  } else {
    matriz <- as.matrix(arggs$n.units)
    arggs$prop1 <- rep(1, ncol(matriz))
    arggs$theta1 <- 0.1
    arggs$polling.sizes <- c(750, 850)
    arggs$n.units <- nrow(arggs$n.units)
    if (!all(matriz >= 0 & matriz == floor(matriz)))
       stop("Argument 'n.units' must be either a positive integer or a matrix of non-negative integers.")
    if (!all(rowSums(matriz) > 0) )
      stop("At least a row in argument 'n.units' has zeroes all its entries.")
  }

  if (!(all(arggs$TM >= 0)))
    stop("Non-negative values are allowed in argument 'TM'.")
  TM <- arggs$TM/rowSums(arggs$TM)

  if(any(is.na(rowSums(TM))))
    stop("At least of row in 'TM' contains NA's.")

  if (length(arggs$prop1) < 2L)
    stop("The argument 'prop1' must have at least length 2.")

  if (!(all(arggs$prop1 > 0)))
    stop("Only positive values are allowed in argument 'prop1'.")
  prop1 <- arggs$prop1/sum(arggs$prop1)

  if (!(all(arggs$polling.sizes > 0)))
    stop("Only positive values are allowed in argument 'polling.sizes'.")

  if (!(length(arggs$polling.sizes) == 2 | length(arggs$polling.sizes) == arggs$n.units))
    stop("The argument 'polling.sizes' must be a vector of length 2 or 'n.units'.")
  if(length(arggs$polling.sizes) == arggs$n.units){
    sizes <- arggs$polling.sizes
  } else {
    sizes <- sample(arggs$polling.sizes[1L]:arggs$polling.sizes[2L],
                    size = arggs$n.units, replace = TRUE)
  }

  if (!(arggs$theta1 >= 0 & arggs$theta1 < 1))
    stop("The argument 'theta1' must be between 0 and 1.")

  if (!(length(arggs$theta2) == 1 | length(arggs$theta2) == nrow(TM)))
    stop("The argument 'theta2' must be of length 1 or R.")

  if (!(all(arggs$theta2 > 0) & all(arggs$theta2 < 1)))
    stop("The argument 'theta2' must be between 0 and 1.")

  if(!(arggs$cs > 0))
    stop("Argument 'cs' must be positive.")

  if (!(length(arggs$noise) == 1 | length(arggs$noise) == nrow(TM)))
    stop("The argument 'noise' must be of length 1 or R.")

  if (!(all(arggs$noise >= 0) & all(arggs$noise < 1)))
    stop("The argument 'noise' must be between 0 and 1.")

  return(list("TM" = TM, "prop1" = prop1, "sizes" = sizes))
}

### rcmult
rcmult <- function(n, p, th){
  if (th == 0){
    votes <- stats::rmultinom(1L, n, p)
  } else {
    votes <- stats::rmultinom(1L, n, rDirichelet(th, p))
  }
  return(votes)
}

### rDirichelet
rDirichelet <- function(theta, p) {
  # Generate a random vector from a Dirichlet distribution
  # with parameters theta (within (0,1)) and p = [p(1),...,p(J)].
  # The p's must be positive and sum to 1.

  # Reparametrize
  alpha <- (1 - theta) / theta
  a <- alpha * p

  # Simulation
  R <- rgamma(length(a), a, 1)

  # Normalize to obtain a Dirichlet random vector
  R <- R / sum(R)
  return(R)
}

### SimElection2
# Function to simulate second election results and transfer matrices
SimElection2 <- function(P, th2, cs, X, noise) {
  # P: Transition probability matrix
  # th: Overdispersion parameter (0, 1)
  # cs: Cluster size
  # N: Results of the first election

  # Number of options and polling stations
  R <- nrow(P)
  C <- ncol(P)
  K <- nrow(X)

  # Generate results for the second election
  Y <- matrix(NA, nrow = K, ncol = C) # Matrix to store second election results
  TM.u <- array(NA, dim = c(R, C, K)) # Array to store cross-distributions of votes
  TM.u2 <- matrix(NA, nrow = K, ncol = R*C) # Matrix to store cross-distributions of votes


  if (max(noise) == 0){
    # Loop over each polling station
    for (kk in 1L:K) {
      # Simulate election results for the polling station
      result <- clusag(n = X[kk, ], P = P, th = th2, cs = cs)
      Y[kk, ] <- result$y # Store column totals
      TM.u[, , kk ] <- result$tm.u
      TM.u2[kk, ] <- as.vector(t(result$tm.u))
    }
  } else {
    # Check if noise is a vector, if not, convert to a vector
    if (length(noise) == 1L) {
      noise <- rep(noise, R)
    }
    for (kk in 1L:K){
      n2 <- pmin(X[kk, ], round(X[kk, ]*noise))
      result1 <- clusag(n = X[kk, ] - n2, P = P, th = th2, cs = cs)
      result2 <- sim_constant(n = n2, nc = C)
      Y[kk, ] <- result1$y + colSums(result2) # Store column totals
      TM.u[, , kk ] <- result1$tm.u + result2
      TM.u2[kk, ] <- as.vector(t(result1$tm.u)) + as.vector(t(result2))
    }
  }

  return(list("Y" = Y, "TM.u" = TM.u, "TM.u2" = TM.u2))
}

#### clusag
clusag <- function(n, P, th, cs) {
  # n: Vector with the number of votes for each option in the first election
  # P: Transition probability matrix
  # th: Overdispersion parameter for each option
  # cs: Cluster size

  # Preliminaries
  R <- nrow(P)
  C <- ncol(P)
  Y <- matrix(NA, nrow = R, ncol = C) # transfer matrix

  # Check if th is a vector, if not, convert to a vector
  if (length(th) == 1L) {
    th <- rep(th, R)
  }

  # Number of clusters and clustered votes
  for (i1 in 1L:R) {
    n.c <- ceiling(n[i1] / cs) # Number of clusters: rounds to the largest integer
    if (n.c > 1) {
      m <- stats::rmultinom(1L, n[i1], rep(1, n.c)/n.c) # Generate clustered votes
    } else {
      m <- n[i1]
    }
    n.c <- length(m)
    Yc <- matrix(NA, nrow = n.c, ncol = C) # Matrix to store votes for each cluster

    for (ic in 1L:n.c) {
      Yc[ic, ] <- rcmult(m[ic], P[i1, ], th[i1]) # Generate votes for each cluster
    }
    Y[i1, ] <- colSums(Yc) # Sum votes across clusters for each option
  }
  y <- colSums(Y) # Total votes across all options
  return(list("y" = y, "tm.u" = Y))
}

### sim_constant
# n2: a vector of length R with the number of votes by row
# C: number of columns
sim_constant <- function(n, nc){
  R <- length(n)
  out <- matrix(0, nrow = R, ncol = nc)
  for (rr in 1L:R){
     out[rr, ] <- stats::rmultinom(1L, n[rr], rep(1/nc, R))
  }
  return(out)
}


### add_noise
# init: an output of SimElection2
# noise: a number between 0 and 1 indicating the additional noise to be included in the model
add_noise <- function(init, noise, prop1){
  K <- nrow(init$TM.u2)
  RC <- ncol(init$TM.u2)
  R <- nrow(init$TM.u[, , 1L])
  C <- ncol(init$TM.u[, , 1L])
  for (kk in 1L:K){
    moves <- max(1, floor(sum(init$Y[kk, ])*noise)) # total moves
    moves.r <- as.vector(stats::rmultinom(1L, moves, prop1)) # moves per row
    matrix0 <- init$TM.u[, , kk]
    quitar <- round(matrix0/rowSums(matrix0)*moves.r)
    quitar[is.nan(quitar)] <- 0
    matrix1 <- matrix0 - quitar
    matrix1[matrix1 < 0] <- 0
    agregar <- matrix(0, nrow = R, ncol = C)
    for (rr in 1L:R){
      agregar[rr, ] <- stats::rmultinom(1L, moves.r[rr], rep(1/C, C))
    }
    init$TM.u[, , kk] <- matrix1 + agregar
    init$Y[kk, ] <- colSums(init$TM.u[, , kk])
    init$TM.u2[kk, ] <- as.vector(t(init$TM.u[, , kk]))
  }
  return(init)
}

