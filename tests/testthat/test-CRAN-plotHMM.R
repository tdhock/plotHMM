library(plotHMM)
library(testthat)
##simulated data.
seg.mean.vec <- c(2, 0, -1, 0)
data.mean.vec <- rep(seg.mean.vec, each=10)
set.seed(1)
N.data <- length(data.mean.vec)
y.vec <- rnorm(N.data, data.mean.vec)
## model.
n.states <- 3
A.mat <- matrix(1/n.states, n.states, n.states)
state.mean.vec <- c(-1, 0, 1)*0.1
sd.param <- 1
pi.vec <- rep(1/n.states, n.states)
## R forward algo.
log.emission.mat <- dnorm(
  y.vec,
  matrix(state.mean.vec, N.data, n.states, byrow=TRUE),
  sd.param,
  log=TRUE)
log.alpha.mat <- matrix(NA, N.data, n.states)
log.alpha.mat[1,] <- elnproduct(
  log(pi.vec),
  log.emission.mat[1,])
for(data.t in 2:N.data){
  log.alpha.vec <- rep(-Inf, n.states)
  for(state.i in 1:n.states){
    prod.vec <- elnproduct(
      log.alpha.mat[data.t-1,state.i],
      log(A.mat[state.i,]))
    log.alpha.vec <- elnsum(log.alpha.vec, prod.vec)
  }
  emission.prob.vec <- log.emission.mat[data.t,]
  log.alpha.mat[data.t,] <- elnproduct(log.alpha.vec, emission.prob.vec)
}
test_that("C++ forward agrees with R", {
  out.list <- plotHMM::forward_interface(log.emission.mat, A.mat, pi.vec)
  expect_equal(out.list$log_alpha, log.alpha.mat)
})

##beta/backward.
log.beta.mat <- matrix(NA, N.data, n.states)
log.beta.mat[N.data, ] <- 0
for(data.t in seq(N.data-1, 1)){
  for(state.i in 1:n.states){
    log.beta <- -Inf
    for(state.j in 1:n.states){
      le <- log.emission.mat[data.t+1,state.j]
      lb <- log.beta.mat[data.t+1,state.j]
      la <- log(A.mat[state.i,state.j])
      inside.product <- elnproduct(le, lb)
      outside.product <- elnproduct(la, inside.product)
      log.beta <- elnsum(log.beta, outside.product)
    }
    log.beta.mat[data.t,state.i] <- log.beta
  }
}
test_that("C++ backward agrees with R", {
  out.mat <- plotHMM::backward_interface(log.emission.mat, A.mat)
  expect_equal(out.mat, log.beta.mat)
})
