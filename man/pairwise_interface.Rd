\name{pairwise_interface}
\alias{pairwise_interface}
\title{Pairwise algorithm}
\description{Efficient implementation of pairwise algorithm in C++ code,
  for N data and S states.}
\usage{pairwise_interface(
  log_emission_mat, transition_mat, log_alpha_mat, log_beta_mat)}
\arguments{
  \item{log_emission_mat, log_alpha_mat, log_beta_mat}{
    N x S numeric matrices of log probabilities.}
  \item{transition_mat}{S x S numeric matrix; transition_mat[i,j] is the
    probability of going from state i to state j.}
}

\value{
  S x S x N numeric array of log probabilities.
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
log.emission.mat <- dnorm(
  y.vec,
  matrix(state.mean.vec, N.data, n.states, byrow=TRUE),
  sd.param,
  log=TRUE)
pi.vec <- rep(1/n.states, n.states)
f.list <- plotHMM::forward_interface(log.emission.mat, A.mat, pi.vec)
b.mat <- plotHMM::backward_interface(log.emission.mat, A.mat)
log.gamma.mat <- plotHMM::multiply_interface(f.list$log_alpha, b.mat)
prob.mat <- exp(log.gamma.mat)
plotHMM::pairwise_interface(log.emission.mat, A.mat, f.list$log_alpha, b.mat)

}
