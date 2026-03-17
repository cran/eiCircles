#' Simulate RxC Tables with Mixed Electoral Behaviours Using Overdispersed Multinomial Models
#'
#' @description  Generates a set of RxC electoral contingency tables under a mixture of voting behaviours, including ecological fallacy effects, within the Overdispersed Multinomial model framework proposed by Forcina et al. (2012), an extension of Brown and Payne (1986). The simulated tables represent the joint distribution of voters in two elections across a set of voting units. Each table is generated using a mixture model that incorporates seven latent voter types, where, consistent with the tradition of mixture models, the number of voters of each type in every unit is assumed to follow a multinomial distribution. The seven electoral behaviours considered (ordinary, faithful, trendy, local retrospective strategic, global retrospective strategic, (global) strategic, and economic voters) are specified in the function's arguments and in **Details**.
#'
#' @author Jose M. Pavia, \email{pavia@@uv.es}
#' @author Antonio Forcina, \email{forcinarosara@@gmail.com}
#'
#' @references Brown, P. and Payne, C. (1986). Aggregate data, ecological regression and voting transitions. *Journal of the American Statistical Association*, 81, 453--460. \doi{10.1080/01621459.1986.10478290}
#' @references Forcina, A., Gnaldi, M. and Bracalente, B. (2012). A revised Brown and Payne model of voting behaviour applied to the 2009 elections in Italy. *Statistical Methods & Applications*, 21, 109--119. \doi{10.1007/s10260-011-0184-x}
#'
#' @param n.units Either a positive integer, `K`, indicating the number of polling units to be simulated, or
#'                a `KxR` data.frame (or matrix) giving the number of votes obtained in election 1
#'                for each of the `R` options in each of the `K` units. If `n.units` is a matrix or data.frame
#'                of counts (votes) the values of arguments `prop1` and `theta1` are ignored.
#'
#' @param TP A `RxC` row-standardized matrix of global transition probabilities for the
#'           Overdispersed Multinomial model (ordinary voters). If not row-standardized,
#'           rows are internally normalized. The row-standardized matrix is represented as
#'           \eqn{\mathbf{P} = [p_{rc}] = [\mathbf{p}^{T}_{r}]}.
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
#' @param theta2 Either a single number between 0 and 1 or a vector of length `nrow(TP)`
#'               containing numbers between 0 and 1. The values in `theta2` serve as
#'               overdispersion parameters and are used alongside the row-probability
#'               vectors in `TP` within the underlying Dirichlet distributions.
#'               These distributions are employed to generate probability vectors
#'               for each combination of unit, cluster, and row, which are then used
#'               to simulate vote transfers from the first to the second election.
#'               If `theta2` is a vector, each row is assigned a distinct overdispersion
#'               parameter based on its corresponding value. Default, `0.1`.
#'
#' @param cs A positive number indicating the average number of cluster size. Default, `50`.
#'
#' @param tau An `Rx7` row-standardized matrix with, by rows (`r`), the vectors of probabilities
#'            of the multinomial distributions used to simulate, in each polling unit,
#'            the number of voters by behaviour type among those who chose option `r` in the first
#'            election. Each cell `(r, t)` defines the probability that a voter who chose option
#'            `r` in the first election behaves as type `t` in the second election.
#'            Probabilities corresponding to electoral behaviours are, by columns, in the order:
#'            ordinary, faithful, trendy, local retrospective strategic, global retrospective
#'            strategic, (global) strategic, and economic voters. If not row-standardized, rows are
#'            internally normalized. The row-standardized matrix is represented as
#'            \eqn{\mathbf{\Theta} = [\tau_{rt}] = [\mathbf{\tau}^{T}_{r}]}. By default,
#'            \eqn{\mathbf{\tau}_{r} =[1,0,0,0,0,0,0]}, i.e., all voters are assumed to behave as
#'            ordinary voters.
#'
#' @param TP.f A `RxC` row-standardized matrix of transition probabilities for faithful
#'             (strongly party-identified) voters, who will vote for the same party again with (almost)
#'             probability 1. The matrix is represented as \eqn{\mathbf{F} = [\mathbf{f}^{T}_{r}]}.
#'             By default, it is (initially) the rectangular identity matrix
#'             of size `RxC` (i.e., \eqn{\mathbf{F} = \mathbf{I}_{RxC}}, where \eqn{(\mathbf{I}_{RxC})_{rc} = 1}
#'             if \eqn{r=c \leq \min(R,C)} and 0 otherwise), assuming the same order for the intersecting
#'             options in the first and second elections. If the same voting options are available
#'             in both elections, this will by default be the identity matrix. If an entire row of
#'             `TP.f` consists of zeroes, it is replaced by the corresponding row in `TP`,
#'             so that faithful voters from the first election linked to the row are still
#'             transferred to the second.
#'
#' @param TP.t A non-negative vector of length `C` that sums to 1, representing the transition probabilities
#'             for trendy voters. If not standardized, the vector is internally normalized.
#'             By default, `TP.t` = \eqn{\mathbf{t}} is the vector of expected results in the second
#'             election implied by the matrix `TP`, assuming that all voters behave as ordinary
#'             voters. Formally, trendy voters behave according to the `RxC` matrix of transition
#'             probabilities: \eqn{\mathbf{1}_{R}\mathbf{t}}.
#'
#' @param LRSV.par A `4xR` matrix of parameters governing the behaviour of local retrospective strategic voters,
#'                 who base their decisions on their party's past results in their polling station.
#'                 The first two rows correspond to minor-party voters and the last two to major-party
#'                 voters. Within each block, the first row corresponds
#'                 to proportion thresholds and the second row to beta parameters. Proportion thresholds
#'                 are non-negative numbers not greater than one, while beta parameters are non-negative
#'                 for minor parties and non-positive for major parties.
#'                 A party cannot be treated simultaneously as both a minor and a major party.
#'                 The sign of beta determines the party's size. See **Details** to understand how the
#'                 parameters are combined to define transition probabilities for these voters. By default,
#'                 all beta parameters are set to zero, which is equivalent to assuming that local
#'                 retrospective strategic voters behave as ordinary voters.
#'
#' @param GRSV.par A `4xR` matrix of parameters governing the behaviour of global retrospective strategic voters,
#'                 who base their decisions on their party's global past results.
#'                 The first two rows correspond to minor-party voters and the last two to major-party
#'                 voters. Within each block, the first row corresponds to proportion thresholds and
#'                 the second row to beta parameters. Proportion thresholds are non-negative numbers
#'                 not greater than one, while beta parameters are non-negative for minor parties
#'                 and non-positive for major parties. A party cannot be treated simultaneously as
#'                 both a minor and a major party. The sign of beta determines the party's size.
#'                 See **Details** to understand how the parameters are combined to define transition
#'                 probabilities for these voters. By default, all beta parameters are set to zero, which is equivalent
#'                 to assuming that global retrospective strategic voters behave as ordinary voters.
#'
#' @param GSV.par A `4xmin(R,C)` matrix of parameters governing the behaviour of global strategic voters,
#'                who base their behaviour on expected results in the second election.
#'                It is assumed that the order of parties in the first and second elections coincides
#'                for the first min(R, C) parties (voting options) in both elections.
#'                The first two rows correspond to minor-party voters and the last two to major-party
#'                voters. Within each block, the first row corresponds to proportion thresholds and
#'                the second row to beta parameters. Proportion thresholds are non-negative numbers
#'                not greater than one, while beta parameters are non-negative for minor parties
#'                and non-positive for major parties. A party cannot be treated simultaneously as both a
#'                minor and a major party. The sign of beta determines the party's size. See **Details**
#'                to understand how the parameters are combined to define transition probabilities for these voters.
#'                By default, all beta parameters are set to zero, which is equivalent to assuming
#'                that strategic voters behave as ordinary voters.
#'
#' @param eco.par A list with three vectors governing the behaviour of economic voters.
#'                These voters prioritise economic performance, rewarding or punishing
#'                parties in the governing coalition based on the perceived local change
#'                in the economic situation. The first component is a vector of length K,
#'                whose elements capture the (perceived) variation in the economy across
#'                voting units, with positive values indicating improvement.
#'                The second component is a vector of length
#'                `R` with the non-negative beta parameters that map the scale of economic
#'                performance to the logits of transition probabilities for each party.
#'                The third component is a vector of length `C`, with entries equal to one
#'                for parties in the second election that were part of the governing coalition
#'                between the first and second elections, and zero otherwise. See **Details** to
#'                understand how the parameters are combined to define transition probabilities for these voters.
#'                By default, all beta parameters are set to zero, which is equivalent to
#'                assuming that economic voters behave as ordinary voters.
#'
#' @param simplify A TRUE/FALSE argument indicating whether the simulated RxCxK array of
#'                 counts by polling unit should be rearranged as a matrix of order Kx(RC).
#'                 Default, FALSE.
#'
#' @param ... Other arguments to be passed to the function. Not currently used.
#'
#'
#' @details Description of how parameters for strategic and economic voters are combined.
#' \itemize{
#'   \item{`local retrospective strategic voters`:} These are voters who consider retrospective outcomes
#'        and make tactical decisions to maximize their preferred outcomes, not necessarily their first choice.
#'        Their decisions are assumed to depend on the local strength of the party they supported in the previous election.
#'        (i) If their party was a minor one, they will support it again when it appears sufficiently strong,
#'        or vote for a different option to avoid wasting their vote;
#'        (ii) If their party was a major one, they will support it again when it seems to require their support
#'        in order to remain strong enough; otherwise, they may choose differently.
#'        Formally, let \eqn{\mathbf{f}_{r}} denote the *r*th row of the matrix \eqn{\mathbf{F}} for faithful voters,
#'        and let \eqn{\mathbf{\lambda}_{r}} denote the vector of logits \eqn{\log(\mathbf{p}_{r}/p_{rC})}
#'        based on the matrix of transition probabilities for ordinary voters.
#'        The vector of retrospective-strategy-local-modified logits for voting unit *s* is defined as
#'        \eqn{\mathbf{\lambda}_{sr}^{LRS} = \mathbf{\lambda}_{r} + \beta_{r}(\pi_{sr} - a_{r})},
#'        where \eqn{a_{r}} is the threshold for party *r*, \eqn{\pi_{sr}} is the proportion of votes
#'        gained by party *r* in voting unit *s* in the first election, and \eqn{\beta_{r}} is the corresponding mapping parameter,
#'        non-negative for minor parties and non-positive for major parties.
#'        In words, \eqn{\mathbf{\lambda}_{r}} is the vector of logits for ordinary voters (representing
#'        basic preferences), \eqn{\pi_{sr}} represents the local strength of party *r* in unit *s*,
#'        \eqn{a_{r}} is the threshold parameter that determines the switching point in voter behaviour,
#'        and \eqn{\beta_{r}} adjusts the degree of strategic consideration.
#'        Under this specification, because of the interaction with the difference \eqn{(\pi_{sr} - a_{r})},
#'        a value of \eqn{\beta_{r} > 0} makes voters more likely to support their party
#'        if it appears sufficiently strong and less likely otherwise,
#'        whereas a value of \eqn{\beta_{r} < 0} makes voters less likely to support their party
#'        if it appears sufficiently strong and more likely otherwise.
#'
#'   \item{`global retrospective strategic voters`:} These voters behave similarly to `local` `retrospective` `strategic
#'        voters`, but consider global rather than local results. They take retrospective outcomes into account
#'        and make tactical decisions to maximize their preferred outcomes, not necessarily their first choice.
#'        Their decisions are assumed to depend on the overall strength of the party they supported in the previous election.
#'        (i) If their party was a minor one, they will support it again when it appears sufficiently strong,
#'        or vote for a different option to avoid wasting their vote;
#'        (ii) If their party was a major one, they will support it again when it seems to require their support
#'        in order to remain strong enough; otherwise, they may choose differently.
#'        Formally, let \eqn{\mathbf{f}_{r}} denote the *r*th row of the matrix \eqn{\mathbf{F}} for faithful voters,
#'        and let \eqn{\mathbf{\lambda}_{r}} denote the vector of logits \eqn{\log(\mathbf{p}_{r}/p_{rC})}
#'        based on the matrix of transition probabilities for ordinary voters.
#'        The vector of retrospective-strategy-global-modified logits is defined as
#'        \eqn{\mathbf{\lambda}_{r}^{GRS} = \mathbf{\lambda}_{r} + \beta_{r}(\pi_{r} - b_{r})},
#'        where \eqn{b_{r}} is the threshold for party *r*, \eqn{\pi_{r}} is the total proportion of votes
#'        gained by party *r* in the first election, and \eqn{\beta_{r}} is the corresponding mapping parameter,
#'        non-negative for minor parties and non-positive for major parties.
#'        In words, \eqn{\mathbf{\lambda}_{r}} is the vector of logits for ordinary voters (representing
#'        basic preferences), \eqn{\pi_{r}} represents the global strength of party *r* in the first election,
#'        \eqn{b_{r}} is the threshold parameter that determines the switching point in voter behaviour,
#'        and \eqn{\beta_{r}} adjusts the degree of strategic consideration.
#'        Under this specification, because of the interaction with the difference \eqn{(\pi_{r} - b_{r})},
#'        a value of \eqn{\beta_{r} > 0} makes voters more likely to support their party
#'        if it appears sufficiently strong and less likely otherwise,
#'        whereas a value of \eqn{\beta_{r} < 0} makes voters less likely to support their party
#'        if it appears sufficiently strong and more likely otherwise.
#'
#'   \item{`global strategic voters`:} These voters behave similarly to `global` `retrospective`
#'        `strategic` `voters`, but base their decisions on expected results in the second election.
#'        They consider expected outcomes and make tactical decisions to maximize their preferred
#'        outcomes, not necessarily their first choice. Their decisions are assumed to depend on the
#'        expected overall strength in the second election of the party they supported in the first election,
#'        knowledge that in practice may be obtained from surveys.
#'        (i) If their party was a minor one, they will support it again when it appears sufficiently strong,
#'        or vote for a different option to avoid wasting their vote;
#'        (ii) If their party was a major one, they will support it again when it seems to require their support
#'        to remain strong enough; otherwise, they may choose differently.
#'        Formally, let \eqn{\mathbf{f}_{r}} denote the *r*th row of the matrix \eqn{\mathbf{F}} for faithful voters,
#'        and let \eqn{\mathbf{\lambda}_{r}} denote the vector of logits \eqn{\log(\mathbf{p}_{r}/p_{rC})}
#'        based on the matrix of transition probabilities for ordinary voters.
#'        Assuming the same order of parties in the first and second elections for those parties affected
#'        by strategic voters, the vector of strategy-global-modified logits is defined as
#'        \eqn{\mathbf{\lambda}_{r}^{GS} = \mathbf{\lambda}_{r} + \beta_{r}\left(\sum_{j}\pi_{j} p_{jr} - c_{r}\right)},
#'        where \eqn{c_{r}} is the threshold for party *r*, \eqn{\pi_{j}} is the total proportion of votes
#'        gained by party *j* in the first election, \eqn{p_{jr}} is the transition probability from party *j*
#'        to party *r* for ordinary voters, and \eqn{\beta_{r}} is the corresponding transforming parameter,
#'        non-negative for minor parties and non-positive for major parties.
#'        In words, \eqn{\mathbf{\lambda}_{r}} is the vector of logits for ordinary voters (representing
#'        basic preferences), \eqn{\sum_{j}\pi_{j} p_{jr}} represents the expected global strength of party *r*
#'        in the second election, \eqn{c_{r}} is the threshold parameter that determines the switching point in voter behaviour,
#'        and \eqn{\beta_{r}} adjusts the degree of strategic consideration.
#'        Under this specification, because of the interaction with the difference \eqn{\sum_{j}(\pi_{j} p_{jr}) - c_{r}},
#'        a value of \eqn{\beta_{r} > 0} makes voters more likely to support their party
#'        if it appears sufficiently strong and less likely otherwise,
#'        whereas a value of \eqn{\beta_{r} < 0} makes voters less likely to support their party
#'        if it appears sufficiently strong and more likely otherwise.
#'
#'   \item{`economic voters`:} These voters prioritise economic performance, rewarding or punishing
#'        parties in the governing coalition based on the perceived change in the local economic situation.
#'        Formally, let \eqn{\mathbf{v}_{s}} denote the perceived measure of economic variation in unit *s*,
#'        \eqn{\mathbf{\lambda}_{r}} the vector of logits \eqn{\log(\mathbf{p}_{r}/p_{rC})} based on the
#'        matrix of transition probabilities for ordinary voters, and \eqn{\mathbf{g}} a vector with entries
#'        equal to 1 for parties in the governing coalition and 0 otherwise. The vector of economically
#'        modified logits for voting unit *s* is then defined as
#'        \eqn{\mathbf{\lambda}_{sr}^{E} = \mathbf{\lambda}_{r} + \beta_{r} \mathbf{v}_{s} \mathbf{g}},
#'        with \eqn{\beta_{r} > 0} being the mapping parameter. Under this specification, these voters are more likely to support
#'        government parties if the local economy improves.
#' }
#'
#'
#' @return
#' A list with the following components
#'
#'  \item{votes1}{ A `KxR` matrix with the (simulated) results in each polling unit for the first election.}
#'  \item{votes2}{ A `KxC` matrix with the simulated results in each polling unit for the second election.}
#'  \item{TM.global}{ An `RxC` matrix with the simulated global transfer matrix of counts.}
#'  \item{TM.units}{ An `RxCxK` array with the simulated transfer matrices of votes by polling unit.
#'                   If `simplify = TRUE`, the simulated transfer matrices of votes are returned
#'                   in a `Kx(RC)` matrix.}
#'  \item{TM.by.behaviour}{ A list with seven components, each of which is itself a list containing the
#'                         four simulated elements (`votes1`, `votes2`, `TM.global` and ` TM.units` as
#'                         RxCxK arrays) corresponding to the subgroups of voters by behaviour, in the
#'                         following order: ordinary, faithful, trendy, local retrospective strategic,
#'                         global retrospective strategic, (global) strategic, and economic voters.}
#'  \item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
#'
#' @export
#'
#' @importFrom stats rmultinom rgamma
#'
#' @note Compared with `simula_BPF_with_deviations`, this function (i) is not restricted to square matrices; (ii) considers up to seven voter types; and (iii) because it mixes distributions, it draws from a distribution with larger variance, even when the latent types and their parameters are the same.
#'
#' @seealso \code{\link{simula_BPF}} \code{\link{simula_BPF_with_deviations}}
#'
#' @examples
#' TMg <- matrix(c(0.6, 0.1, 0.3, 0.1, 0.7, 0.2, 0.1, 0.1, 0.8),
#'              byrow = TRUE, nrow = 3)
#' example <- simula_mixture(n.units = 100, TP = TMg, prop1 = c(0.3, 0.3, 0.4),
#'                                       polling.sizes = c(750, 850))


simula_mixture <- function(n.units,
                           TP,
                           prop1,
                           polling.sizes,
                           theta1 = 0.1,
                           theta2 = 0.1,
                           cs = 50,
                           tau,
                           TP.f,
                           TP.t,
                           LRSV.par,
                           GRSV.par,
                           GSV.par,
                           eco.par,
                           simplify = FALSE,
                           ...
){


  inputs <- c(as.list(environment()), list(...))
  arggs0 <- inputs
  if(missing(tau)){
    n.behave <- 7L
    arggs0$tau <- cbind(rep(1, nrow(TP)),
                        matrix(0, nrow = nrow(TP), ncol = n.behave - 1L)  )
  }
  arggs <- tests_inputs_simula_mixture(arggs0)

  TP <- arggs$TP
  prop1 <- arggs$prop1
  sizes <- arggs$sizes
  tau <- arggs$tau

  # K <- n.units
  if(is.null(dim(n.units))){
    K <- n.units
    votes1 <- matrix(NA, nrow = K, ncol = nrow(TP))
    for (kk in 1L:K){
      votes1[kk, ] <- rcmult(sizes[kk], prop1, theta1)
    }
  } else {
    votes1 <- as.matrix(n.units)
    K <- nrow(votes1)
  }

  # Voters by behaviour type
  votes1.bh <- array(0, dim = c(dim(votes1), ncol(tau)))
  for (kk in 1:K){
    for (rr in 1:ncol(votes1)){
      votes1.bh[kk, rr, ] <- stats::rmultinom(1, votes1[kk,rr], tau[rr, ])
    }
  }

  # Ordinary voters
  votes1.o <- votes1.bh[, , 1L]
  out.o <- votes_ord(voters = votes1.o, TP = TP,
                      theta = theta2, cs = cs)[-5L]

  # Faithful voters
  votes1.f <- votes1.bh[, , 2L]
  if (missing(TP.f)){
      TP.f <- matrix(0, nrow = nrow(TP), ncol = ncol(TP))
      diag(TP.f) <- 1
  }
  if (nrow(TP.f) != nrow(TP) | ncol(TP.f) != ncol(TP))
    stop("The dimension of 'TP.f' should be equal to the dimension of 'TP'.")

  if(any(is.na(rowSums(TP.f))))
    stop("At least of row in 'TP.f' contains NA's.")

  if (!(all(TP.f >= 0)))
    stop("Negative values are allowed in argument 'TP.f'.")

  zero.r <- which(rowSums(TP.f) == 0)
  TP.f <- TP.f/rowSums(TP.f)
  TP.f[zero.r, ] <- TP[zero.r, ]

  out.f <- votes_faitful(voters = votes1.f, TP.f = TP.f,
                         theta = theta2, cs = cs, n.c = ncol(TP))

  # Trendy voters
  votes1.t <- votes1.bh[, , 3L]
  if (missing(TP.t))
    TP.t <- t(colSums(votes1)/sum(votes1)) %*% TP
  out.t <- votes_trendy(voters = votes1.t, TP = TP.t,
                         theta = theta2, cs = cs, n.c = ncol(TP))

  # Local retrospective strategic voters
  votes1.lr <- votes1.bh[, , 4L]
  pii <- votes1/rowSums(votes1)
  if (missing(LRSV.par))
    LRSV.par <- rbind(rep(1, nrow(TP)),
                      rep(0, nrow(TP)),
                      rep(1, nrow(TP)),
                      rep(0, nrow(TP)))
  out.lr <- votes_local(voters = votes1.lr, param = LRSV.par, pii = pii, TP = TP,
                        TP.f = TP.f, theta = theta2, cs = cs, n.c = ncol(TP))

  # Global retrospective strategic voters
  votes1.gr <- votes1.bh[, , 5L]
  piig <- colSums(votes1)/sum(votes1)
  if (missing(GRSV.par))
    GRSV.par <- rbind(rep(1, nrow(TP)),
                      rep(0, nrow(TP)),
                      rep(1, nrow(TP)),
                      rep(0, nrow(TP)))
  out.gr <- votes_global_r(voters = votes1.gr, param = GRSV.par, pii = piig, TP = TP,
                           TP.f = TP.f, theta = theta2, cs = cs, n.c = ncol(TP))

  # Global strategic voters
  votes1.s <- votes1.bh[, , 6L]
  if (missing(GSV.par))
    GSV.par <- rbind(rep(1, min(dim(TP))),
                    rep(0, min(dim(TP))),
                      rep(1, min(dim(TP))),
                      rep(0, min(dim(TP))))
  out.s <- votes_strat(voters = votes1.s, param = GSV.par, pii = piig, TP = TP,
                           TP.f = TP.f, theta = theta2, cs = cs, n.c = ncol(TP))

  # Economic voters
  votes1.e <- votes1.bh[, , 7L]
  if (missing(eco.par))
    eco.par <- list(rep(0, nrow(votes1)),
                    rep(0, nrow(TP)),
                    rep(0, ncol(TP)))
  out.e <- votes_eco(voters = votes1.e, param = eco.par, TP = TP,
                       theta = theta2, cs = cs, n.c = ncol(TP))

  # Outputs: Consolidation
  TM.units <- out.o$TM.units + out.f$TM.units + out.t$TM.units +
    out.lr$TM.units + out.gr$TM.units + out.s$TM.units +
    out.e$TM.units
  votes2 <- apply(TM.units, c(3, 2), sum)
  TM.units2 <- t(apply(TM.units, 3, function(x) as.vector(t(x))))
  TM.behaviour <- list("ordinary" = out.o, "faithful" = out.f, "trendy" = out.t,
                       "local" = out.lr, "global" = out.gr, "strategic" = out.s,
                       "economic" = out.e)

  names1 <- paste0("R", 1L:nrow(TP))
  names2 <- paste0("C", 1L:ncol(TP))
  names3 <- paste0("unit", 1L:K)

  # Improving outputs presentation
  rownames(votes1) <- rownames(votes2) <- rownames(TM.units2) <- names3
  colnames(votes1) <- names1
  colnames(votes2) <- names2
  dimnames(TM.units) <- list(names1, names2, names3)

  TM.global <- apply(TM.units, c(1L, 2L), sum)

  if(simplify){
    TM.units <- TM.units2
    names4 <- paste0(rep(names1, each = nrow(TP)), rep(names2, ncol(TP)))
    colnames(TM.units) <- names4
  }

  output <- list("votes1" = votes1, "votes2" = votes2,
                 "TM.global" = TM.global, "TM.units" = TM.units,
                 "TM.by.behaviour" = TM.behaviour, "inputs" = inputs)

  class(output) <- "simula_BPF"
  return(output)
}


######################################
# Auxiliary functions

### tests_inputs_simula_mixture
tests_inputs_simula_mixture <- function(arggs){

  n.behave <- 7L
  if(is.null(dim(arggs$n.units))){
    if(!(arggs$n.units > 0 & (arggs$n.units - round(arggs$n.units) == 0)))
      stop("Argument 'n.units' must be either a positive integer or a matrix of non-negative integers.")
  } else {
    matriz <- as.matrix(arggs$n.units)
    arggs$prop1 <- rep(1, ncol(matriz))
    arggs$theta1 <- 0.1
    arggs$polling.sizes <- rowSums(matriz)
    arggs$n.units <- nrow(arggs$n.units)
    if (!all(matriz >= 0 & matriz == floor(matriz)))
      stop("Argument 'n.units' must be either a positive integer or a matrix of non-negative integers.")
    if (!all(rowSums(matriz) > 0) )
      stop("At least a row in argument 'n.units' has zeroes all its entries.")
  }

  if (!(all(arggs$TP >= 0)))
    stop("Negative values are allowed in argument 'TP'.")
  TP <- arggs$TP/rowSums(arggs$TP)

  if(any(is.na(rowSums(TP))))
    stop("At least of row in 'TP' contains NA's or all its elements are zero.")

  if (length(arggs$prop1) != nrow(TP))
    stop("The argument 'prop1' must have the same length as the number of rows in 'TP'.")

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

  if (!(length(arggs$theta2) == 1 | length(arggs$theta2) == nrow(TP)))
    stop("The argument 'theta2' must be of length 1 or R.")

  if (!(all(arggs$theta2 > 0) & all(arggs$theta2 < 1)))
    stop("The argument 'theta2' must be between 0 and 1.")

  if(!(arggs$cs > 0))
    stop("Argument 'cs' must be positive.")

  if (is.null(arggs$tau))
    stop("Argument 'tau' must be of dim Rx7")

  if (ncol(arggs$tau) != n.behave | nrow(arggs$tau) != nrow(arggs$TP) | any(arggs$tau < 0) | any(rowSums(arggs$tau) == 0))
    stop("Argument 'tau' must be of dim Rx7, with non-negative components and at least one positive component in each row.")
  tau <- arggs$tau/rowSums(arggs$tau)

  if (any(is.na(arggs$tau)))
      stop("NA's are not allowed in 'tau'.")

  return(list("TP" = TP, "prop1" = prop1, "sizes" = sizes,
              "tau" = tau))
}

## Transition of ordinary voters

votes_ord <- function(voters, TP, theta, cs){
  if (sum(voters) == 0){
    out <- list("votes1" = voters,
                "votes2" = matrix(0, nrow = nrow(voters), ncol = ncol(TP)),
                "TM.global" = matrix(0, nrow = ncol(voters), ncol = ncol(TP)),
                "TM.units" = array(0, dim = c(ncol(voters), ncol(TP), nrow(voters)))
               )
  } else {
    votes2 <- matrix(0, nrow = nrow(voters), ncol = nrow(TP))
    TM.units <- array(0, dim = c(ncol(voters), nrow(TP), nrow(voters)))
    sel <- (rowSums(voters) != 0)
    out <- simula_BPF(n.units = voters[sel, ], TM = TP,
                      theta2 = theta, cs = cs)[-5L]
    votes2[sel, ] <- out$votes2
    TM.units[, , sel] <- out$TM.units
    out$votes2 <- votes2
    out$TM.units <- TM.units
    out$votes1 <- voters
  }
  return(out)
}

## Transitions of faithful voters
votes_faitful <- function(voters, TP.f, theta, cs, n.c){
  if (sum(voters) == 0){
    out <- list("votes1" = voters,
                "votes2" = matrix(0, nrow = nrow(voters), ncol = n.c),
                "TM.global" = matrix(0, nrow = ncol(voters), ncol = n.c),
                "TM.units" = array(0, dim = c(ncol(voters), n.c, nrow(voters)))
               )
  } else {
    votes2 <- matrix(0, nrow = nrow(voters), ncol = nrow(TP.f))
    TM.units <- array(0, dim = c(ncol(voters), nrow(TP.f), nrow(voters)))
    sel <- (rowSums(voters) != 0)
    out <- simula_BPF(n.units = voters[sel, ], TM = TP.f,
                      theta2 = theta, cs = cs)[-5L]
    votes2[sel, ] <- out$votes2
    TM.units[, , sel] <- out$TM.units
    out$votes2 <- votes2
    out$TM.units <- TM.units
    out$votes1 <- voters
  }
  return(out)
}

## Transitions of trendy voters
votes_trendy <- function(voters, TP, theta, cs, n.c){
  if (sum(voters) == 0){
    out <- list("votes1" = voters,
                "votes2" = matrix(0, nrow = nrow(voters), ncol = n.c),
                "TM.global" = matrix(0, nrow = ncol(voters), ncol = n.c),
                "TM.units" = array(0, dim = c(ncol(voters), n.c, nrow(voters)))
    )
  } else {
    if (length(TP) != n.c)
      stop("The length of 'TP.t' should be equal to the number of columns of 'TP'.")

    if(any(is.na(TP)))
      stop("At least of component in 'TP.t' is NA.")

    if (!(all(TP >= 0)))
      stop("Negative values are allowed in argument 'TP.t'.")
    TP.t <- TP/sum(TP)
    TP.t <- matrix(rep(TP.t, each = ncol(voters)), nrow = ncol(voters), byrow = FALSE)

    votes2 <- matrix(0, nrow = nrow(voters), ncol = nrow(TP.t))
    TM.units <- array(0, dim = c(ncol(voters), nrow(TP.t), nrow(voters)))
    sel <- (rowSums(voters) != 0)
    out <- simula_BPF(n.units = voters[sel, ], TM = TP.t,
                      theta2 = theta, cs = cs)[-5L]
    votes2[sel, ] <- out$votes2
    TM.units[, , sel] <- out$TM.units
    out$votes2 <- votes2
    out$TM.units <- TM.units
    out$votes1 <- voters
   }
  return(out)
}


## Transitions of local retrospective strategic voters
votes_local <- function(voters, param, pii, TP, TP.f, theta, cs, n.c){
  if (sum(voters) == 0){
    out <- list("votes1" = voters,
                "votes2" = matrix(0, nrow = nrow(voters), ncol = n.c),
                "TM.global" = matrix(0, nrow = ncol(voters), ncol = n.c),
                "TM.units" = array(0, dim = c(ncol(voters), n.c, nrow(voters)))
    )
  } else {
    if (any(dim(param) != c(4L,ncol(voters))))
      stop("The dimension of 'LRSV.par' should be 4x(the number of rows of 'TP').")
    if(any(is.na(param)))
      stop("At least of component in 'LRSV.par' is NA.")

    if (!(all(param[c(1L:3L), ] >= 0)))
      stop("Negative values are allowed in rows 1 to 3 of 'LRSV.par'.")

    if (!(all(param[c(1L, 3L), ] <= 1)))
      stop("Values greater than 1 are allowed in rows 1 and 3 of 'LRSV.par'.")

    if (!(all(param[4L, ] <= 0)))
      stop("Positive values are allowed in row 4 of 'LRSV.par'.")

    if (max((param[4L, ] < 0) + (param[2L, ] > 0)) == 2)
      stop("Error in 'LRSV.par': a party cannot simultaneously be considered both big and small.")

    # TP.f <- TP.f/rowSums(TP.f)
    votes2 <- matrix(0, nrow = nrow(voters), ncol = nrow(TP))
    TM.units <- array(0, dim = c(ncol(voters), nrow(TP), nrow(voters)))
    sel <- which(rowSums(voters) != 0)
    out <- simula_BPF(n.units = voters[sel, ], TM = TP, theta2 = theta, cs = cs)[-5L]
    for (kk in sel){
      newTP <- updating_TP(pii = pii[kk, ], param = param, TP = TP, TP.f = TP.f)
      out.k <- simula_BPF(n.units = t(as.matrix(voters[kk, ])), TM = newTP,
                          theta2 = theta, cs = cs)[-5L]
      out$votes2[kk, ] <- out.k$votes2
      out$TM.units[, , kk] <- out.k$TM.units
    }
    out$TM.global <- apply(out$TM.units, c(1L, 2L), sum)
    votes2[sel, ] <- out$votes2
    TM.units[, , sel] <- out$TM.units
    out$votes2 <- votes2
    out$TM.units <- TM.units
    out$votes1 <- voters
  }
  return(out)
}

# updating_TP
# Function to updating TP considering retrospective strategic voters
updating_TP <- function(pii, param, TP, TP.f){
  newTP <- TP
  for (rr in 1L:ncol(param)){
    tp <- TP[rr, ]
    tp[tp == 0] <- exp(-32)
    tpm <- which.max(tp)
    lambda <- log(tp/tp[tpm])
    lambda <- lambda + param[2L, rr]*(pii[rr] - param[1L, rr])*TP.f[rr, ]
    lambda <- lambda + param[4L, rr]*(pii[rr] - param[3L, rr])*TP.f[rr, ]
    newTP[rr, ] <- exp(lambda)/sum(exp(lambda))
  }
  return(newTP)
}

## Transitions of global retrospective strategic voters
votes_global_r <- function(voters, param, pii, TP, TP.f, theta, cs, n.c){
  if (sum(voters) == 0){
    out <- list("votes1" = voters,
                "votes2" = matrix(0, nrow = nrow(voters), ncol = n.c),
                "TM.global" = matrix(0, nrow = ncol(voters), ncol = n.c),
                "TM.units" = array(0, dim = c(ncol(voters), n.c, nrow(voters)))
    )
  } else {
    if (any(dim(param) != c(4L,ncol(voters))))
      stop("The dimension of 'GRSV.par' should be 4x(the number of rows of 'TP').")
    if(any(is.na(param)))
      stop("At least of component in 'GRSV.par' is NA.")

    if (!(all(param[c(1L:3L), ] >= 0)))
      stop("Negative values are allowed in rows 1 to 3 of 'GRSV.par'.")

    if (!(all(param[c(1L, 3L), ] <= 1)))
      stop("Values greater than 1 are allowed in rows 1 and 3 of 'GRSV.par'.")

    if (!(all(param[4L, ] <= 0)))
      stop("Positive values are allowed in row 4 of 'GRSV.par'.")

    if (max((param[4L, ] < 0) + (param[2L, ] > 0)) == 2)
      stop("Error in 'GRSV.par': a party cannot simultaneously be considered both big and small.")

    #TP.f <- TP.f/rowSums(TP.f)
    newTP <- updating_TP(pii = pii, param = param, TP = TP, TP.f = TP.f)

    votes2 <- matrix(0, nrow = nrow(voters), ncol = nrow(TP))
    TM.units <- array(0, dim = c(ncol(voters), nrow(TP), nrow(voters)))
    sel <- (rowSums(voters) != 0)
    out <- simula_BPF(n.units = voters[sel, ], TM = newTP, theta2 = theta, cs = cs)[-5L]
    votes2[sel, ] <- out$votes2
    TM.units[, , sel] <- out$TM.units
    out$votes2 <- votes2
    out$TM.units <- TM.units
    out$votes1 <- voters
  }
  return(out)
}

## Transitions of global strategic voters
votes_strat <- function(voters, param, pii, TP, TP.f, theta, cs, n.c){
  if (sum(voters) == 0){
    out <- list("votes1" = voters,
                "votes2" = matrix(0, nrow = nrow(voters), ncol = n.c),
                "TM.global" = matrix(0, nrow = ncol(voters), ncol = n.c),
                "TM.units" = array(0, dim = c(ncol(voters), n.c, nrow(voters)))
    )
  } else {
    if (any(dim(param) != c(4L, min(dim(TP)))))
      stop("The dimension of 'GSV.par' is innappropiate, it should be 4x(min(dim(TP))).")

    if(any(is.na(param)))
      stop("At least of component in 'GSV.par' is NA.")

    if (!(all(param[c(1L:3L), ] >= 0)))
      stop("Negative values are allowed in rows 1 to 3 of 'GSV.par'.")

    if (!(all(param[c(1L, 3L), ] <= 1)))
      stop("Values greater than 1 are allowed in rows 1 and 3 of 'GSV.par'.")

    if (!(all(param[4L, ] <= 0)))
      stop("Positive values are allowed in row 4 of 'GSV.par'.")

    if (max((param[4L, ] < 0) + (param[2L, ] > 0)) == 2)
      stop("Error in 'GSV.par': a party cannot simultaneously be considered both big and small.")

    #TP.f <- TP.f/rowSums(TP.f)
    pii <- pii %*% TP
    newTP <- updating_TP(pii = pii, param = param, TP = TP, TP.f = TP.f)

    votes2 <- matrix(0, nrow = nrow(voters), ncol = nrow(TP))
    TM.units <- array(0, dim = c(ncol(voters), nrow(TP), nrow(voters)))
    sel <- (rowSums(voters) != 0)
    out <- simula_BPF(n.units = voters[sel, ], TM = newTP,
                      theta2 = theta, cs = cs)[-5L]
    votes2[sel, ] <- out$votes2
    TM.units[, , sel] <- out$TM.units
    out$votes2 <- votes2
    out$TM.units <- TM.units
    out$votes1 <- voters
  }
  return(out)
}


## Transitions of economic voters

votes_eco <- function(voters, param, TP, theta, cs, n.c){
  if (sum(voters) == 0){
    out <- list("votes1" = voters,
                "votes2" = matrix(0, nrow = nrow(voters), ncol = n.c),
                "TM.global" = matrix(0, nrow = ncol(voters), ncol = n.c),
                "TM.units" = array(0, dim = c(ncol(voters), n.c, nrow(voters)))
    )
  } else {
    if (length(param) != 3L)
      stop("The length of 'eco.par' is innappropiate, it should a list of length 3.")

    # Component 1
    if(any(is.na(param[[1]])))
      stop("At least of component in the first vector of 'eco.par' is NA.")

    if(length(param[[1]]) != nrow(voters))
      stop("The length of the first vector of 'eco.par' must be K.")

    # Component 2
    if(any(is.na(param[[2]])))
      stop("At least of component in the second vector of 'eco.par' is NA.")

    if (!(all(param[[2]] >= 0)))
      stop("Negative values are allowed in the second vector of 'eco.par'.")

    if(length(param[[2]]) != ncol(voters))
      stop("The length of the second vector of 'eco.par' must be equal to the number of rows in 'TP'.")

    # Component 3
    if(any(is.na(param[[3]])))
      stop("At least of component in the third vector of 'eco.par' is NA.")

    if (!all(param[[3]] %in% c(0, 1)))
      stop("Only entries 0 or 1 allowed in the third vector of 'eco.par'.")

    if(length(param[[3]]) != n.c)
      stop("The length of the third vector of 'eco.par' must be equal to the number of columns in 'TP'.")

    votes2 <- matrix(0, nrow = nrow(voters), ncol = nrow(TP))
    TM.units <- array(0, dim = c(ncol(voters), nrow(TP), nrow(voters)))
    sel <- which(rowSums(voters) != 0)
    out <- simula_BPF(n.units = voters[sel, ], TM = TP, theta2 = theta, cs = cs)[-5L]
    for (rr in sel){
      newTP <- updating_TP_eco(param = param, TP = TP, rr = rr)
      out.r <- simula_BPF(n.units = t(as.matrix(voters[rr, ])), TM = newTP,
                          theta2 = theta, cs = cs)[-5L]
      out$votes2[rr, ] <- out.r$votes2
      out$TM.units[, , rr] <- out.r$TM.units
    }
    out$TM.global <- apply(out$TM.units, c(1L, 2L), sum)
    votes2[sel, ] <- out$votes2
    TM.units[, , sel] <- out$TM.units
    out$votes2 <- votes2
    out$TM.units <- TM.units
    out$votes1 <- voters
  }
  return(out)

}

# updating_TP
# Function to updating TP considering economic voters
updating_TP_eco <- function(param, TP, rr){
  newTP <- TP
  for (rr in 1L:nrow(TP)){
    tp <- TP[rr, ]
    tp[tp == 0] <- exp(-32)
    tpm <- which.max(tp)
    lambda <- log(tp/tp[tpm])
    lambda <- lambda + param[[2]][rr]*param[[1]][rr]*param[[3]]
     newTP[rr, ] <- exp(lambda)/sum(exp(lambda))
  }
  return(newTP)
}
