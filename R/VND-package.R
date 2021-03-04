#' @docType package
#' @details The VND model is used for a superpoistion of Markov chains that are permutation invariant and conditional independent. 
#' This package allows to estimate parameters and simulate such a model in a HMM setting.
#' @section Simulation:
#' Simulation of data from superposition models, including the VND model, is done with the following functions
#' \describe{
#' \item{\code{vnd_model}}{returns vnd transition matrix from given parameters.}
#' \item{\code{ck_model}}{returns Chung-Kennedy transition matrix from given parameters.}
#' \item{\code{uc_model}}{returns uncoupled transition matrix from given parameters.}
#' \item{\code{simulate_MC}}{simulates data from a HMM with a given transition matrix.}
#' }
#' @section Estimation:
#' Estimation of VND parameters from a VND-HMM is done with the following functions
#' \describe{
#' \item{\code{HMM_analysis}}{returns the estimated model fit, Viterbi path, and log-likelihood for a VND model with gaussian emisions.}
#' \item{\code{HMM_estimation}}{returns the estimated model fit for a VND model with gaussian emisions.}
#' \item{\code{HMM_normal}}{returns the estimated model fit from a custom model with gaussian emisions.}
#' \item{\code{HMM_custom}}{returns the estimated model fit from a custom model with custom emisions.}
#' \item{\code{calculate_log_likelihood}}{returns the estimated the log likelihood for a given model fit.}
#' }
#'
"_PACKAGE"

## usethis namespace: start
#' @useDynLib VND, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats constrOptim
## usethis namespace: end

NULL