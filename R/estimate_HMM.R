# Copyright (C) 2019-2020 Benjamin Eltzner
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
# or see <http://www.gnu.org/licenses/>.

# data:                   vector
# matrix_pars_estimation: function(matrix, function(vector) -> matrix, integer)
#                         -> list('pars': vector, 'loss': float)
# matrix_from_pars:       function(vector) -> matrix
# matrix_pars_init:       vector
# emit_estimation:        function(vector, vector, vector, integer) -> vector
# emit_generation:        function(vector, vector, integer) -> matrix
# emit_init:              vector
# max_it:                 integer
# threshold:              float

#' HMM model fit for Gaussian data with custom transition matrix, custom emission distribution, and associated losses.
#' 
#' Gives a model fit for multiple channels with a custom model of the transition matrix, custom emission distribution, and custom parameter optimization steps.
#' @name HMM_custom
#' @param data numeric vector containing given data.
#' @param matrix_pars_estimation function optimizing the parameters of a tranistion matrix model according to a custom loss function. The function should return 
#' a list of a vector \code{pars} with the optimized parameters, and a double \code{loss} with the minimized loss value.
#' @param matrix_from_pars function describing the transition matrix model given a vector of parameters. The function should return a squared transition matrix.
#' @param matrix_pars_init vector of parameters to initialize the \code{matrix_from_pars} function.
#' @param emit_estimation function optimizing the parameters of an emission model according to a custom loss function. The function should return a vector of
#' the optimized parameters.
#' @param emit_generation function calculating the probabilities of one data point to come from each hidden state, for a given distribution.
#' The function should return a vector of length equal to the sum state space with the probabilities of that data point to be in the corresponding state.
#' @param emit_init vector containing the initial parameters of the emission distribution. The two first components are the base level and step level size. The next \code{n_channels}+1 entries are the initial standard deviations for each level.
#' @param max_it integer containing the maximal number of iterations of the Baum-Welch-type algorithm. Set to defaul of 100.
#' @param threshold double containing the threshold of convergence for the Baum-Welch-type algorithm  
#' @return A list containing a vector \code{matrix_pars} of the estimated parameters for the given \code{model},
#' the estimated transition matrix \code{matrix}, a vector \code{mu} of the estimated levels, a vector \code{sigma} of the estimated standard deviations for each level,
#' a vector \code{pi} of the estimated initial distribution, and a double \code{loglik} of the pseudo-log-likelihood of this model.
#' @export
HMM_custom <- function(data, # vector
                       matrix_pars_estimation, # function which returns list('pars': vector, 'loss': number)
                       matrix_from_pars, #
                       matrix_pars_init,
                       emit_estimation,
                       emit_generation,
                       emit_init,
                       max_it=100,
                       threshold=1e-6) {
  it <- 0;
  matrix_pars <- matrix_pars_init;
  n_matrix_pars <- length(matrix_pars);
  trans_matrix <- matrix_from_pars(matrix_pars);
  emit <- emit_init;
  current_error = 100;
  n_levels <- dim(trans_matrix);
  if ((length(n_levels) < 2) || (n_levels[1] != n_levels[2])) {
    print("ERROR: Output of transition matrix generation function is not square!");
    return(NULL);
  }
  n_levels <- n_levels[1];
  n_emit_pars <- length(emit_init);

  while ((it < max_it) && (current_error > threshold)){
    it <- it+1;

    probabilities <- emit_generation(data, emit, n_levels);

    # Forward-Backward optimization
    tmp <- BW_step(data, n_levels, probabilities, trans_matrix);

    # Parameter Calculation
    emit_new <- emit_estimation(tmp$gamma, tmp$y_gamma, tmp$y_gamma2, n_emit_pars);
    #return(list(emit, emit_new))
    matrix_fit <- matrix_pars_estimation(tmp$xi, matrix_from_pars, n_matrix_pars, matrix_pars);
    matrix_pars <- matrix_fit$pars;
    matrix_new <- matrix_from_pars(matrix_pars);

    # Calculate parameter change
    current_error = max(max(abs(matrix_new - trans_matrix)), max(abs(emit_new - emit)));
    emit <- emit_new;
    trans_matrix <- matrix_new;
  }

  return(list("matrix_pars" = matrix_pars, "matrix" = trans_matrix,
              "emit" = emit, "pi" = tmp$pi, 'matrix_loss' = matrix_fit$loss,
              "weights" = tmp$gamma));
}

#' HMM model fit for Gaussian data with custom transition matrix
#' 
#' Gives a model fit for multiple channels with a custom model of the transition matrix for a given data for a HMM with Gaussian emissions and equidistant levels.
#' @name HMM_normal
#' @param data numeric vector containing given data.
#' @param matrix_from_pars function describing the transition matrix model. The function should return a squared transition matrix.
#' @param matrix_pars_init vector of parameters to initialize the \code{matrix_from_pars} function.
#' @param emit_init vector containing the initial parameters of the emission distribution. The two first components are the base level and step level size. The next \code{n_channels}+1 entries are the initial standard deviations for each level.
#' @param max_it integer containing the maximal number of iterations of the Baum-Welch-type algorithm. Set to defaul of 100.
#' @param threshold double containing the threshold of convergence for the Baum-Welch-type algorithm  
#' @return A list containing a vector \code{matrix_pars} of the estimated parameters for the given \code{model},
#' the estimated transition matrix \code{matrix}, a vector \code{mu} of the estimated levels, a vector \code{sigma} of the estimated standard deviations for each level,
#' a vector \code{pi} of the estimated initial distribution, and a double \code{loglik} of the pseudo-log-likelihood of this model.
#' @export

HMM_normal <- function(data,
                       matrix_from_pars,
                       matrix_pars_init,
                       emit_init,
                       max_it=100,
                       threshold=1e-6) {

  # Determine initial emission parameters
  trans_matrix <- matrix_from_pars(matrix_pars_init);
  n_levels <- dim(trans_matrix)[1];
  step <- (emit_init[2] - emit_init[1]) / (n_levels - 1);
  emit_init <- c(emit_init[1], step, rep(10/step, n_levels));

  # Run generic function
  tmp <- HMM_custom(data, matrix_likelihood, matrix_from_pars, matrix_pars_init,
                    normal_params, normal_probabilities, emit_init, max_it, threshold);

  # Get parameters in standard form
  emit <- tmp$emit;
  mu <- emit[1] + 0:(n_levels-1) * emit[2];
  sigma <- 1 / abs(emit[-(1:2)]);
  log_likelihood <- tmp$matrix_loss + sum(tmp$weights*log(sigma))

  return(list("matrix_pars" = tmp$matrix_pars, "matrix" = tmp$matrix,
              "mu" = mu, "sigma" = sigma, "pi" = tmp$pi,
              "loglik" = log_likelihood));
}

#' HMM model fit for Gaussian data
#' 
#' Gives a model fit for multiple channels with a given data and model for a HMM with Gaussian emissions and equidistant levels.
#' @name HMM_estimate
#' @param data numeric vector containing given data.
#' @param model string of model name to fit to the data. Currently supported "VND", Chung-Kennedy "CK", and uncoupled "independent"
#' @param n_channels integer of number of channels
#' @param emit_init vector containing the initial parameters of the emission distribution. The two first components are the base level and step level size. The next \code{n_channels}+1 entries are the initial standard deviations for each level.
#' @param max_it integer containing the maximal number of iterations of the Baum-Welch-type algorithm. Set to defaul of 100.
#' @param threshold double containing the threshold of convergence for the Baum-Welch-type algorithm  
#' @return A list containing a vector \code{matrix_pars} of the estimated parameters for the given \code{model},
#' the estimated transition matrix \code{matrix}, a vector \code{mu} of the estimated levels, a vector \code{sigma} of the estimated standard deviations for each level,
#' a vector \code{pi} of the estimated initial distribution, and a double \code{loglik} of the pseudo-log-likelihood of this model.
#' @export

HMM_estimate <- function(data,
                         model,
                         n_channels,
                         emit_init,
                         max_it=100,
                         threshold=1e-6) {

  if (model == 'independent') {
    matrix_from_pars <- make_uc_function(n_channels);
    matrix_pars_init <- rep(0.9, 2);
  } else if (model == 'CK') {
    matrix_from_pars <- make_ck_function(n_channels);
    matrix_pars_init <- c(0.9, 0.9,0.1);
  } else if (model == 'VND') {
    matrix_from_pars <- vnd_model;
    matrix_pars_init <- rep(0.9, 2*n_channels);
  } else {
    print(cat("Unknown model type: ", model,
              ". Supported model types: 'independent', 'CK', 'VND'. Defaulting to independent!", sep=""));
    matrix_from_pars <- make_uc_function(n_channels);
    matrix_pars_init <- rep(0.9, 2);
  }

  # Run generic function
  tmp <- HMM_normal(data, matrix_from_pars, matrix_pars_init, emit_init,
                    max_it=100, threshold=1e-6);

  return(tmp);
}

#' Calculate log-likelihood for model fit
#' 
#' Gives the log-likelihood for a given model fit for multiple channels with a given data for gaussian emisions
#' 
#' @param data numeric vector containing given data.
#' @param model_fit a list containing the model_fit. This model fit should contain the estimated mean \code{mu}, estimated variance \code{sigma}, and estimated transition matrix \code{matrix}.
#' @return a double containing the log likelihood of the model for the data.
#' @export

calculate_log_likelihood <- function(data, model_fit) {
  mu <- model_fit$mu;
  sigma <- model_fit$sigma;
  trans_matrix <-model_fit$matrix;
  probabilities <- exp(-0.5*t(t(outer(data,mu,"-")^2)/(sigma^2))) / sigma;
  probabilities <- t(probabilities/rowSums(probabilities));
  return(Log_Likelihood(data, model_fit$pi, probabilities, trans_matrix)$lik);
}

#' HMM model analysis
#' 
#' Gives a complete HMM analysis for multiple channels with a given data and model for a HMM with Gaussian emissions and equidistant levels.
#' @name HMM_analysis
#' @param data numeric vector containing given data
#' @param model string of model name to fit to the data. Currently supported "VND", Chung-Kennedy "CK", and uncoupled "independent"
#' @param n_channels integer of number of channels
#' @param emit_init vector containing the initial parameters of the emission distribution. The two first components are the base level and step level size. The next \code{n_channels}+1 entries are the initial standard deviations for each level.
#' @param max_it integer containing the maximal number of iterations of the Baum-Welch-type algorithm. Set to default to 100.
#' @param threshold double containing the threshold of convergence for the Baum-Welch-type algorithm  
#' @return A list containing the estimated model fit \code{model_fit}, Viterbi path \code{viterbi_path}, and log-likelihood \code{log_likelihood}.\cr
#' \code{model_fit} is a list containing a vector \code{matrix_pars} of the estimated parameters for the given \code{model},
#' the estimated transition matrix \code{matrix}, a vector \code{mu} of the estimated levels, a vector \code{sigma} of the estimated standard deviations for each level,
#' a vector \code{pi} of the estimated initial distribution, and a double \code{loglik} of the pseudo-log-likelihood of this model.\cr
#' \code{viterbi_path} is a vector of length equal to \code{data} with the fitted underlying signal. \code{log_likelihood} is a double of the log-likelihood of this model.
#' @export

HMM_analysis <- function(data,
                         model,
                         n_channels,
                         emit_init,
                         max_it=100,
                         threshold=1e-6) {
  model_fit <- HMM_estimate(data, model, n_channels, emit_init, max_it=100,
                            threshold=1e-6);
  viterbi_path <- Viterbi_simple(data, n_channels+1, model_fit$pi,
                                 model_fit$matrix, model_fit$mu,
                                 model_fit$sigma^2);
  log_likelihood <- calculate_log_likelihood(data, model_fit);
  return(list('model_fit' = model_fit, 'viterbi_path' = viterbi_path,
              'log_likelihood' = log_likelihood));
}

################################################################################
##########################    Auxiliary functions    ###########################
################################################################################

make_uc_function <- function (n_channels) {
  uc_function <- function(parameters) {
    return(uc_model(parameters, n_channels));
  }
  return(uc_function);
}

make_ck_function <- function (n_channels) {
  ck_function <- function(parameters) {
    return(ck_model(parameters, n_channels));
  }
  return(ck_function);
}

matrix_3state <- function (p) {
  a <- c(p[1]*p[2],     p[1]*(1-p[2]), 1-p[1],
         p[3]*(1-p[4]), p[3]*p[4],     1-p[3],
         1-p[5],        p[5]*(1-p[6]), p[5]*p[6]);
  return(t(matrix(a, nrow=3)));
}

#mlikel
matrix_likelihood <- function(BW_matrix, model, n_par, init){
  sol <- constrOptim(init, log_cost, P_trans = BW_matrix, model= model, NULL,
                     ui = rbind(diag(n_par),-diag(n_par)),
                     ci = c(rep(0,n_par), rep(-1,n_par)));
  return(list('pars' = sol$par, 'loss' = sol$value));
}

#log cost
log_cost <- function(x, P_trans, model){
  target <- as.vector(P_trans);
  tmp <- model(x);
  tmp[tmp<=0] = 10^(-10);
  model_vec <- as.vector(log(tmp));
  return(-sum(target*model_vec));
}

normal_probabilities  <- function(data, emit, n_levels) {
  mu <- emit[1] + 0:(n_levels-1) * emit[2];
  sigma <- 1 / abs(emit[-(1:2)]);
  expon <- -0.5*t(outer(data,mu,"-")^2)/(sigma^2);
  probabilities <- exp(expon - matrix(rep(apply(expon,2,max), n_levels), nrow=n_levels, byrow=TRUE)) / sigma;
  return(t(t(probabilities)/colSums(probabilities)));
}

normal_params <- function(g, y, y2, n_par){
  n <- length(g);
  init <- y/g;
  step <- (init[n]-init[1])/(n-1);
  init <- c(init[1], step, rep(10/step, n));
  Mcost <- function(x, arg1, arg2, arg3){
    mu = x[1] + 0:(n-1) * x[2];
    tmp = sum((arg3-2*mu*arg2+arg1*mu^2)*x[-(1:2)]^2 - arg1*log(x[-(1:2)]^2));
    return(tmp);
  }
  sol <- constrOptim(init, Mcost, arg1 = g, arg2 = y, arg3 = y2, NULL,
                     ui = diag(n+2)[2:(n+2),], ci = c(rep(0,n+1)));
  return(sol$par);
}
