% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simula_BPF_with_deviations.R
\name{simula_BPF_with_deviations}
\alias{simula_BPF_with_deviations}
\title{Simulate RxC Square Tables with Ecological Fallacy Effects Based on Overdispersed-Multinomial Models}
\usage{
simula_BPF_with_deviations(
  n.units,
  TM,
  prop1,
  polling.sizes,
  theta1 = 0.1,
  theta2 = 0.1,
  cs = 50,
  prop.dev = c(0.4, 0.6),
  prop.loyal = matrix(0.34, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)),
    ncol = nrow(TM)),
  prop.strategic = matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units,
    nrow(n.units)), ncol = nrow(TM)),
  prop.context = matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units,
    nrow(n.units)), ncol = nrow(TM)),
  par.loyal = 0.95,
  par.strategic = 0.5,
  par.context = 0.5,
  simplify = FALSE,
  ...
)
}
\arguments{
\item{n.units}{Either a positive integer number, \code{K}, indicating the number of polling units to be simulated, or
a \code{KxR} data.frame of a matrix with the number of votes gained in election 1
for each of the \code{R} options in each of the \code{K} units. If \code{n.units} is a matrix (data.frame) of
counts (votes) the values of arguments \code{prop1} and \code{theta1} are ommitted.}

\item{TM}{A row-standardized RxC matrix with the underlying global transition
probabilities for the Overdispersed-Multinomial Model. If the matrix is not row-standardized,
it is internally row-standardized by the function.}

\item{prop1}{A vector of length R with the initial assumed probabilities of voting (to be simulated)
for each of the R competing options in the first election. If the provided vector
is not a set of probabilities (i.e., a vector of positive numbers adding to 1),
it is internally standardized by the function.}

\item{polling.sizes}{Either a vector of two components with two positive integer
numbers indicating the minimum and maximum number of voters
for each unit or a vector of length \code{n.units} of positive integer
numbers informing about the number of voters in each unit. When
\code{polling.sizes} is a vector of length two, a number of voters is
randomly assigned for each unit using a uniform distribution with
parameters the minimum and maximum values included in \code{polling.sizes}.}

\item{theta1}{A number between 0 and 1 used as the overdispersion parameter.
This parameter is employed by the underlying Dirichlet distribution,
in conjunction with \code{prop1}, to randomly generate vectors of probabilities
for each unit. These vectors are then used to simulate the results
of the first election. The smaller the value of this parameter,
the closer the unit-level marginal distributions for the first election
are to \code{prop1}. Default, \code{0.1}.}

\item{theta2}{Either a single number between 0 and 1 or a vector of length \code{nrow(TM)}
containing numbers between 0 and 1. The values in \code{theta2} serve as
overdispersion parameters and are used alongside the row-probability
vectors in \code{TM} within the underlying Dirichlet distributions.
These distributions are employed to generate probability vectors
for each combination of unit, cluster, and row, which are then used
to simulate vote transfers from the first to the second election.
If \code{theta2} is a vector, each row is assigned a distinct overdispersion
parameter based on its corresponding value. Default, \code{0.1}.}

\item{cs}{A positive number indicating the average number of cluster size. Default, \code{50}.}

\item{prop.dev}{Either a two-component vector with positive values between 0 and 1,
indicating the minimum and maximum proportion of voters (to be simulated)
that deviate from the base Overdispersed-Multinomial Model in each unit
or a vector of length \code{n.units} specifying the proportion of voters deviating
from the basic model in each unit. If \code{prop.dev} is a two-component vector,
the proportion of deviating voters in each unit is randomly assigned using
a uniform distribution with the specified minimum and maximum values.
Default, \code{c(0.4, 0.6)}.}

\item{prop.loyal}{A KxR matrix where each cell \verb{(k, r)} represents the proportion of voters from party
\code{r} in unit \code{k} who are strongly loyal. These voters are highly likely to vote for
the same party with near certainty (see the parameter \code{par.loyal}).
In contrast, the remaining \code{prop.dev} percent of the voters from the party
follow the transition probabilities specified in \code{TM}. The sum of the matrices \code{prop.loyal},
\code{prop.strategic}, and \code{prop.contextual} must equal one for each cell.
If this condition is not met, the function internally standardizes
the provided matrices. Default, \code{matrix(0.34, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM))}.}

\item{prop.strategic}{A KxR matrix where each cell \verb{(k, r)} represents the proportion of voters
from party \code{r} in unit \code{k} who are strategic voters. These voters are a
\code{par.strategic} percent more likely to support parties that improve their
results in the second election compared to their performance in their
first election (see the parameter \code{par.strategic}). In contrast, the remaining
\code{prop.dev} percent of the voters from the party  follow the transition
probabilities specified in \code{TM}. The sum of the matrices \code{prop.loyal},
\code{prop.strategic}, and \code{prop.contextual} must equal one for each cell.
If this condition is not met, the function internally
standardizes the provided matrices.
Default, \code{matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM))}.}

\item{prop.context}{A KxR matrix where each cell \verb{(k, r)} represents the proportion of voters
from party \code{r} in unit \code{k} who are influenced by the relative strength
in their neighborhood of the party they voted for in the first election.
These voters are a \code{par.context} multiplied by the party's strength
in the unit percent more likely to support the same party in the second
election (see the parameter \code{par.context}). In contrast, the remaining
\code{prop.dev} percent of the voters from the party  follow the transition
probabilities specified in \code{TM}. The sum of the matrices \code{prop.loyal},
\code{prop.strategic}, and \code{prop.contextual} must equal one for each cell.
If this condition is not met, the function internally
standardizes the provided matrices.
Default, \code{matrix(0.33, nrow = ifelse(is.null(dim(n.units)), n.units, nrow(n.units)), ncol = nrow(TM))}.}

\item{par.loyal}{A number between 0.9 and 1 indicating the minimum probability with which loyal
voters will support the same party in the second election as they did in the
first. For each unit, the probability is randomly chosen between \code{par.loyal}
and 1. Default, \code{0.95}.}

\item{par.strategic}{A positive number indicating the proportion of increase that
the initial transfer probabilities in \code{TM} should be increased
for those parties improving their support in the second election
compared to their performance in their first election. Default, \code{0.5}.}

\item{par.context}{A positive number indicating the factor by which the proportion of
support for a party in each unit should be multiplied to increase
the initial transfer probabilities in \code{TM} corresponding to that party.
Default, \code{0.5}.}

\item{simplify}{A TRUE/FALSE argument indicating whether the simulated RxCxK array of
counts by polling unit should be rearranged as a matrix of order Kx(RC).
Default, FALSE.}

\item{...}{Other arguments to be passed to the function. Not currently used.}
}
\value{
A list with the following components

\item{votes1}{ A matrix of order KxR with the results simulated in each polling unit for the first election.}
\item{votes2}{ A matrix of order KxC with the results simulated in each polling unit for the second election..}
\item{TM.global}{ A matrix of order RxC with the actual simulated global transfer matrix of counts.}
\item{TM.units}{ An array of order RxCxK with the simulated transfer matrices of votes by polling unit. If
\code{simplify = TRUE} the simulated transfer matrices of votes are returned organized in a Kx(RC) matrix.}
\item{inputs}{ A list containing all the objects with the values used as arguments by the function.}
}
\description{
Generates a set of RxC square (RxR) tables at random, representing the joint distribution of voters in two elections, according to the model proposed by Forcina et al. (2012) as an extension of Brown and Payne (1986), under the assumption that transition probabilities are non-homogeneous across local units. For each unit, a unique transition table is constructed to simulate voter behavior within that unit. Each table is created using a mixture model that considers four latent types of voters: one group following the underlying global transition probabilities of the BPF model, another composed mainly of loyal voters, a third characterized by strategic voting, and a final group whose probability of loyalty to the party they supported in the first election depends on that party's strength in the unit during the first election
}
\examples{
TMg <- matrix(c(0.6, 0.1, 0.3, 0.1, 0.7, 0.2, 0.1, 0.1, 0.8),
             byrow = TRUE, nrow = 3)
example <- simula_BPF_with_deviations(n.units = 100, TM = TMg, prop1 = c(0.3, 0.3, 0.4),
                                      polling.sizes = c(750, 850))
}
\references{
Brown, P. and Payne, C. (1986). Aggregate data, ecological regression and voting transitions. \emph{Journal of the American Statistical Association}, 81, 453–460. \doi{10.1080/01621459.1986.10478290}

Forcina, A., Gnaldi, M. and Bracalente, B. (2012). A revised Brown and Payne model of voting behaviour applied to the 2009 elections in Italy. \emph{Statistical Methods & Applications}, 21, 109–119. \doi{10.1007/s10260-011-0184-x}
}
\seealso{
Other simulators for ecological inference overdispersed-multinomial models: 
\code{\link{simula_BPF}()}
}
\author{
Antonio Forcina, \email{forcinarosara@gmail.com}

Jose M. Pavia, \email{pavia@uv.es}
}
\concept{simulators for ecological inference overdispersed-multinomial models}
