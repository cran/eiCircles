#' Simulate RxC Square Tables with Ecological Fallacy Effects Based on Overdispersed-Multinomial Models
#'
#' @description  Generates a set of RxC square (RxR) tables at random, representing the joint distribution of voters in two elections, according to the model proposed by Forcina et al. (2012) as an extension of Brown and Payne (1986), under the assumption that transition probabilities are non-homogeneous across local units. For each unit, a unique transition table is constructed to simulate voter behavior within that unit. Each table is created using a mixture model that considers four latent types of voters: one group following the underlying global transition probabilities of the BPF model, another composed mainly of loyal voters, a third characterized by strategic voting, and a final group whose probability of loyalty to the party they supported in the first election depends on that party's strength in the unit during the first election
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
#'           probabilities for the Overdispersed-Multinomial Model. If the matrix is not row-standardized,
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
#' @param prop.dev Either a two-component vector with positive values between 0 and 1,
#'                 indicating the minimum and maximum proportion of voters (to be simulated)
#'                 that deviate from the base Overdispersed-Multinomial Model in each unit
#'                 or a vector of length `n.units` specifying the proportion of voters deviating
#'                 from the basic model in each unit. If `prop.dev` is a two-component vector,
#'                 the proportion of deviating voters in each unit is randomly assigned using
#'                 a uniform distribution with the specified minimum and maximum values.
#'                 Default, `c(0.4, 0.6)`.
#'
#' @param prop.loyal A KxR matrix where each cell `(k, r)` represents the proportion of voters from party
#'                   `r` in unit `k` who are strongly loyal. These voters are highly likely to vote for
#'                   the same party with near certainty (see the parameter `par.loyal`).
#'                   In contrast, the remaining `prop.dev` percent of the voters from the party
#'                   follow the transition probabilities specified in `TM`. The sum of the matrices `prop.loyal`,
#'                   `prop.strategic`, and `prop.contextual` must equal one for each cell.
#'                   If this condition is not met, the function internally standardizes
#'                   the provided matrices. Default, `matrix(0.34, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM))`.
#'
#' @param prop.strategic A KxR matrix where each cell `(k, r)` represents the proportion of voters
#'                       from party `r` in unit `k` who are strategic voters. These voters are a
#'                       `par.strategic` percent more likely to support parties that improve their
#'                       results in the second election compared to their performance in their
#'                       first election (see the parameter `par.strategic`). In contrast, the remaining
#'                       `prop.dev` percent of the voters from the party  follow the transition
#'                       probabilities specified in `TM`. The sum of the matrices `prop.loyal`,
#'                       `prop.strategic`, and `prop.contextual` must equal one for each cell.
#'                       If this condition is not met, the function internally
#'                       standardizes the provided matrices.
#'                       Default, `matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM))`.
#'
#' @param prop.context A KxR matrix where each cell `(k, r)` represents the proportion of voters
#'                     from party `r` in unit `k` who are influenced by the relative strength
#'                     in their neighborhood of the party they voted for in the first election.
#'                     These voters are a `par.context` multiplied by the party's strength
#'                     in the unit percent more likely to support the same party in the second
#'                     election (see the parameter `par.context`). In contrast, the remaining
#'                     `prop.dev` percent of the voters from the party  follow the transition
#'                     probabilities specified in `TM`. The sum of the matrices `prop.loyal`,
#'                     `prop.strategic`, and `prop.contextual` must equal one for each cell.
#'                     If this condition is not met, the function internally
#'                     standardizes the provided matrices.
#'                     Default, `matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM))`.
#'
#' @param par.loyal A number between 0.9 and 1 indicating the minimum probability with which loyal
#'                  voters will support the same party in the second election as they did in the
#'                  first. For each unit, the probability is randomly chosen between `par.loyal`
#'                  and 1. Default, `0.95`.
#'
#' @param par.strategic A positive number indicating the proportion of increase that
#'                      the initial transfer probabilities in `TM` should be increased
#'                      for those parties improving their support in the second election
#'                      compared to their performance in their first election. Default, `0.5`.
#'
#' @param par.context A positive number indicating the factor by which the proportion of
#'                    support for a party in each unit should be multiplied to increase
#'                    the initial transfer probabilities in `TM` corresponding to that party.
#'                    Default, `0.5`.
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
#' example <- simula_BPF_with_deviations(n.units = 100, TM = TMg, prop1 = c(0.3, 0.3, 0.4),
#'                                       polling.sizes = c(750, 850))


simula_BPF_with_deviations <- function(n.units,
                                       TM,
                                       prop1,
                                       polling.sizes,
                                       theta1 = 0.1,
                                       theta2 = 0.1,
                                       cs = 50,
                                       prop.dev = c(0.4, 0.6),
                                       prop.loyal = matrix(0.34, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM)),
                                       prop.strategic = matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM)),
                                       prop.context = matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM)),
                                       par.loyal = 0.95,
                                       par.strategic = 0.5,
                                       par.context = 0.5,
                                       simplify = FALSE,
                                       ...
){

  inputs <- c(as.list(environment()), list(...))
  arggs <- tests_inputs_simula_BF_with_deviations(inputs)
  TM <- arggs$TM
  prop1 <- arggs$prop1
  sizes <- arggs$sizes
  # K <- n.units
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
  prop.dev <- arggs$prop.dev
  sum.prop <- prop.loyal + prop.strategic + prop.context
  prop.loyal <- prop.loyal/sum.prop
  prop.strategic <- prop.strategic/sum.prop
  prop.context <- prop.context/sum.prop
  prop.loyal[is.nan(prop.loyal)] <- 1/3
  prop.strategic[is.nan(prop.strategic)] <- 1/3
  prop.context[is.nan(prop.context)] <- 1/3
  par.loyal <- array(stats::runif(length(prop.loyal), par.loyal, 1),
                      dim = dim(prop.loyal))
  strat <- which(colSums(TM*prop1) > prop1) # parties benefiting from strategic voting

  TM.units <- array(NA, dim = c(nrow(TM), ncol(TM), K)) # Array to store cross-distributions of votes
  TM.units2 <- matrix(NA, nrow = K, ncol = prod(dim(TM)))
#  votes1 <- matrix(NA, nrow = K, ncol = nrow(TM))
  votes2 <- matrix(NA, nrow = K, ncol = ncol(TM))
  for (kk in 1L:K){
   # votes1[kk, ] <- rcmult(sizes[kk], prop1, theta1)
    TMk.loyal <- TMk_loyal(par.loyal[kk, ], TM) * prop.loyal[kk, ]
    TMk.strategic <- TMk_strategic(par.strategic, TM, strat) * prop.strategic[kk, ]
    TMk.context <- TMk_context(par.context, TM, votes1[kk, ]) * prop.context[kk, ]
    TMk <- TMk.loyal + TMk.strategic + TMk.context
    TMk <- (1 - prop.dev[kk]) * TM + prop.dev[kk] * TMk
    uk <- SimElection2_unit(P = TMk, th2 = theta2, cs = cs, x = votes1[kk, ])
    TM.units[, , kk] <- uk$TM.u
    TM.units2[kk, ] <- uk$TM.u2
    votes2[kk, ] <- uk$Y
  }

  names1 <- paste0("R", 1L:nrow(TM))
  names2 <- paste0("C", 1L:ncol(TM))
  names3 <- paste0("unit", 1L:K)

  # Improving outputs presentation
  rownames(votes1) <- rownames(votes2) <- rownames(TM.units2) <- names3
  colnames(votes1) <- names1
  colnames(votes2) <- names2
  dimnames(TM.units) <- list(names1, names2, names3)

  TM.global <- apply(TM.units, c(1L, 2L), sum)

  if(simplify){
    TM.units <- TM.units2
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
tests_inputs_simula_BF_with_deviations <- function(arggs){

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

  if (ncol(arggs$TM) != nrow(arggs$TM))
    stop("The 'TM' matrix must be a square matrix.")

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

  if(arggs$par.loyal > 1 | arggs$par.loyal < 0.9)
    stop("Argument 'par.loyal' must be between 0.9 and 1.")

  if(arggs$par.strategic <= 0)
    stop("Argument 'par.strategic' must be positive.")

  if(arggs$par.context <= 0)
    stop("Argument 'par.context' must be positive.")

  if(nrow(arggs$prop.loyal) != arggs$n.units |
     ncol(arggs$prop.loyal) != nrow(arggs$TM))
    stop("Argument 'prop.loyal' does not have the proper structure.")

  if(nrow(arggs$prop.strategic) != arggs$n.units |
     ncol(arggs$prop.strategic) != nrow(arggs$TM))
    stop("Argument 'prop.loyal' does not have the proper structure.")

  if(nrow(arggs$prop.context) != arggs$n.units |
     ncol(arggs$prop.context) != nrow(arggs$TM))
    stop("Argument 'prop.context' does not have the proper structure.")

  if (any(is.na(arggs$prop.loyal)))
    stop("NA's are allowed in 'prop.loyal'.")

  if (any(is.na(arggs$prop.strategic)))
    stop("NA's are not allowed in 'prop.strategic'.")

  if (any(is.na(arggs$prop.context)))
    stop("Negative values are not allowed in 'prop.context'.")

  if (any(arggs$prop.loyal < 0))
    stop("Negative values are not allowed in 'prop.loyal'.")

  if (any(arggs$prop.strategic < 0))
    stop("Negative values are not allowed in 'prop.strategic'.")

  if (any(arggs$prop.context < 0))
    stop("Negative values are allowed in 'prop.context'.")

  if (!(length(arggs$prop.dev) == 2 | length(arggs$prop.dev) == arggs$n.units))
    stop("The argument 'prop.dev' must be a vector of length 2 or 'n.units'.")
  if(length(arggs$prop.dev) == arggs$n.units){
    prop.dev <- arggs$prop.dev
  } else {
    prop.dev <- stats::runif(arggs$n.units, arggs$prop.dev[1L], arggs$prop.dev[2L])
  }


  return(list("TM" = TM, "prop1" = prop1, "sizes" = sizes,
              "prop.dev" = prop.dev))
}


### TMk_loyal
### Function to construct a transition table for loyal voters
TMk_loyal <- function(p.loyal, TM){
  p.loyal <- runif(1, p.loyal, 1)
  out <- TM
  diag(out) <- 0
  out <- out*(1 - pmax(p.loyal, diag(TM)))/rowSums(out)
  diag(out) <- pmax(p.loyal, diag(TM))
  return(out)
}

### TTMk_strategic
### Function to construct a transition table for strategic voters
TMk_strategic <- function(par.strategic, TM, strat){
  out <- TM
  if(length(strat) > 0){
    out[, strat] <- 0
    out.s <- TM[, strat, drop = FALSE] * (1 + par.strategic)
    out.s[out.s > 1] <- 1
    out[, -strat] <- out[, -strat, drop = FALSE] *
      (1 - rowSums(out.s))/rowSums(out[, -strat, drop = FALSE])
    out[, strat] <- out.s
  }
  return(out)
}

### TTMk_strategic
### Function to construct a transition table for voters influenced by context
TMk_context <- function(par.context, TM, votes1){
  out <- TM
  diag(out) <- 0
  p.cont <- pmin(1, diag(TM) * (1 + par.context * votes1/sum(votes1)))
  out <- out * (1 - p.cont)/rowSums(out)
  diag(out) <- p.cont
  return(out)
}

### SimElection2
# Function to simulate second election results and transfer matrices
SimElection2_unit <- function(P, th2, cs, x) {
  # P: Transition probability matrix
  # th: Overdispersion parameter (0, 1)
  # cs: Cluster size
  # x: Results of the first election

  # Simulate election results for the polling station
  result <- clusag(n = x, P = P, th = th2, cs = cs)
  TM.u <- result$tm.u
  TM.u2 <- as.vector(t(result$tm.u))
  Y <- colSums(TM.u)

  return(list("Y" = Y, "TM.u" = TM.u, "TM.u2" = TM.u2))
}

