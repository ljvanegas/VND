# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

BW_step <- function(data, n_levels, probabilities, P_trans_init) {
    .Call(`_VND_BW_step`, data, n_levels, probabilities, P_trans_init)
}

Log_Likelihood <- function(data, init, probabilities, P_trans) {
    .Call(`_VND_Log_Likelihood`, data, init, probabilities, P_trans)
}

#' Viterbi path
#' 
#' Produces a the estimated viterbi path for given data and parameters.
#' @name Viterbi_simple
#' @param y numeric vector of data
#' @param K integer containing the number of states
#' @param pi vector containing the initial probability of the system
#' @param P_trans KxK matrix of tranisitiong probabilities
#' @param mu vector containing the level for each state
#' @param var vector containing the variance for each state
#' @return vector of length equal to the length of the data containing the viterbi path
#' @export
NULL

Viterbi_simple <- function(y, K, pi, P_trans, mu, var) {
    .Call(`_VND_Viterbi_simple`, y, K, pi, P_trans, mu, var)
}

#' Transition matrix according to the Chung-Kennedy model
#' 
#' Produces a transition matrix acording to the Chung-Kennedy (CK) model for a given vector of parameters.
#' @name ck_model
#' @param parameters vector with parameters of the ck model. 
#' @param l integer containing the number of signals (channels).
#' @return the transition matrix for a ck-model with the given parameters \eqn{\lambda}, \eqn{\eta}, \eqn{k}
#' @examples
#' # For l=2
#' q<-ck_model(c(0.9,0.8,0.5), 2)
#' q
#' @export
NULL

ck_model <- function(p, l) {
    .Call(`_VND_ck_model`, p, l)
}

mydnorm <- function(x, means, sds) {
    .Call(`_VND_mydnorm`, x, means, sds)
}

density_pois <- function(x, lambda) {
    .Call(`_VND_density_pois`, x, lambda)
}

mydnorm2 <- function(x, means, sds) {
    .Call(`_VND_mydnorm2`, x, means, sds)
}

mydnorm3 <- function(x, means, sds) {
    .Call(`_VND_mydnorm3`, x, means, sds)
}

mydnorm4 <- function(x, means, sds) {
    .Call(`_VND_mydnorm4`, x, means, sds)
}

mydnorm5 <- function(x, means, sds, K, r, beta) {
    .Call(`_VND_mydnorm5`, x, means, sds, K, r, beta)
}

mydnorm6 <- function(x, means, sds, K) {
    .Call(`_VND_mydnorm6`, x, means, sds, K)
}

mydnorm7 <- function(x, means, sds, beta) {
    .Call(`_VND_mydnorm7`, x, means, sds, beta)
}

#' Simulation of Markov chains
#' 
#' Simulates \code{n_data} data points from a Markov chain with given transition matrix and intial probability.
#' @name simulate_MC
#' @param n_data integer containing the number of data points to be simulated.
#' @param P_trans transition matrix to be simulated from. It should be a square matrix.
#' @param intitial vector of length equal to the number of rows of \code{P_trans} containing the intital probability distribution.
#' @return vector of length \code{n_data} containing the simulated markov chain with state space {0,1,...,\code{nrow(P_trans)}}.
#' @examples
#' q <- vnd_model(c(0.9,0.8,0.98,0.89))
#' mc <- simulate_MC(100, q, c(1/3,1/3,1/3))
#' plot(mc, type = "s")
#' @export
NULL

simulate_MC <- function(n_data, P_trans, initial) {
    .Call(`_VND_simulate_MC`, n_data, P_trans, initial)
}

#' Transition matrix according to the uncoupled model
#' 
#' Produces a transition matrix acording to the uncoupled (independent) model for a given vector of parameters.
#' @name uc_model
#' @param parameters vector with parameters of the uc model. 
#' @param l integer containing the number of signals (channels).
#' @return the transition matrix for a uc-model with the given parameters \eqn{\lambda}, \eqn{\eta}
#' @examples
#' # For l=2
#' q <- uc_model(c(0.9,0.8), 2)
#' q
#' @export
NULL

uc_model <- function(p, l) {
    .Call(`_VND_uc_model`, p, l)
}

#' Transition matrix according to the VND model
#' 
#' Produces a transition matrix acording to the VND model for a given vector of parameters.
#' @name vnd_model
#' @param parameters vector with parameters of the VND model \eqn{(\lambda_0,\eta_1,\lamda_1,\ldtos,\eta_\ell)}. 
#' @return the transition matrix for a VND model with the given parameters
#' @examples
#' # For l=2
#' q<-vnd_model(c(0.9,0.8,0.98,0.89))
#' q
#' @export
NULL

vnd_model <- function(parameters) {
    .Call(`_VND_vnd_model`, parameters)
}

