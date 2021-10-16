\name{forward_interface}
\alias{forward_interface}
\title{Forward algorithm}
\description{Efficient implementation of forward algorithm in C++ code,
  for N data and S states.}
\usage{forward_interface(
  log_emission_mat, initial_prob_vec, transition_mat)}
\arguments{
  \item{log_emission_mat}{N x S numeric matrix of log probabilities of
    observing each data point in each state.}
  \item{initial_prob_vec}{S numeric vector of probabilities of observing
    each state at the beginning of the sequence.}
  \item{transition_mat}{S x S numeric matrix; transition_mat[i,j] is the
    probability of going from state i to state j.}
}

\value{list with two elements
\item{log_alpha}{N x S numeric matrix of forward log probabilities.}
\item{log_lik}{numeric scalar log likelihood of data given model parameters.}
}

\author{Toby Dylan Hocking}

\examples{

##simulated data.
seg.mean.vec <- c(2, 0, -1, 0)
data.mean.vec <- rep(seg.mean.vec, each=10)
set.seed(1)
N.data <- length(data.mean.vec)
y.vec <- rnorm(N.data, data.mean.vec)
##model.
n.states <- 3
A.mat <- matrix(1/n.states, n.states, n.states)
state.mean.vec <- c(-1, 0, 1)*0.1
sd.param <- 1
pi.vec <- rep(1/n.states, n.states)
log.emission.mat <- dnorm(
  y.vec,
  matrix(state.mean.vec, N.data, n.states, byrow=TRUE),
  sd.param,
  log=TRUE)
plotHMM::forward_interface(log.emission.mat, pi.vec, A.mat)

}
