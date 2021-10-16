\name{backward_interface}
\alias{backward_interface}
\title{Backward algorithm}
\description{Efficient implementation of backward algorithm in C++ code,
  for N data and S states.}
\usage{backward_interface(
  log_emission_mat, transition_mat)}
\arguments{
  \item{log_emission_mat}{N x S numeric matrix of log probabilities of
    observing each data point in each state.}
  \item{transition_mat}{S x S numeric matrix; transition_mat[i,j] is the
    probability of going from state i to state j.}
}

\value{N x S numeric matrix of backward log probabilities.}

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
log.emission.mat <- dnorm(
  y.vec,
  matrix(state.mean.vec, N.data, n.states, byrow=TRUE),
  sd.param,
  log=TRUE)
plotHMM::backward_interface(log.emission.mat, A.mat)

}
