#' Nigel Bastian Cendra, S2610920
#' Add your own function definitions on this file.

#' neg_log_lik
#
#' @description Evaluate the negated log-likelihood for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model

neg_log_lik <- function(beta, data, model){
  
  mu <- beta[1] + beta[2]*data[["CAD_Weight"]]
  
  # distinguish between the two models to find the particular standard deviation for the betas
  if(model == "A") {
    sigma <- sqrt(exp(beta[3] + beta[4]*data[["CAD_Weight"]]))
  }else{
    sigma <- sqrt(exp(beta[3])+exp(beta[4]) * (data[["CAD_Weight"]]^2))
  }
  - sum(dnorm(data[["Actual_Weight"]],
              mean = mu,
              sd=sigma,
              log = TRUE))
  
}

#' filament_estimate
#
#' @description Estimate filament models with different variance structure
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @return An estimation object suitable for use with [filament1_predict()]

filament1_estimate <- function(data, model) {
  model <- match.arg(model, c("A", "B"))
  if (model == "A") {
    beta_start <- c(-0.1, 1.07, -2, 0.05)
  } else {
    beta_start <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(beta_start,
               neg_log_lik,
               data = data,
               model = model,
               hessian = TRUE,
               method = "Nelder-Mead",
               control = list(maxit = 5000)
  )
  fit <- list(
    model = model,
    par = opt$par,
    hessian = opt$hessian
  )
  class(fit) <- c("filament1_estimate", "list")
  fit
}

#' filament1_aux_EV
#' 
#' @description Evaluate the expectation and variance for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` containing the required predictors, including `CAD_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @param Sigma_beta : If not NULL, an estimate of the covariance matrix for
#                 the uncertainty of estimated betas
#' @return A list with four elements:
#     E : E(y|beta,x)
#     V : Var(y|beta,x)
#     VE : Var(E(y|beta,x)|x) or NULL
#     EV : E(Var(y|beta,x)|x) or NULL

filament1_aux_EV <- function(beta, data, model = c("A", "B"),
                             Sigma_beta = NULL) {
  
  model <- match.arg(model)
  if (model == "A") {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      EV <- exp(ZV %*% beta + rowSums(ZV * (ZV %*% Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = exp(ZV %*% beta),
      VE = VE,
      EV = EV
    )
  } else {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + I(CAD_Weight^2), data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      # (pmin: Ignore large Sigma_beta values)
      EV <- ZV %*% exp(beta + pmin(0.5^2, diag(Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = ZV %*% exp(beta),
      VE = VE,
      EV = EV
    )
  }
  out
}


#' Predict filament weights using specified model
#'
#' @description Compute predictive distributions and prediction intervals for a new dataset using either Model A or Model B.
#' @param data A `data.frame` used to estimate the model parameters.
#' @param model A character string specifying the model to use, either "A" or "B".
#' @param newdata A `data.frame` containing new observations for prediction.
#' @param alpha The significance level for prediction interval computation (default is 0.05).
#' @return A `data.frame` with the predicted mean, standard deviation, lower and upper bounds of the prediction interval.

filament1_predict <- function(data, model, newdata, alpha = 0.05){
  fit <- filament1_estimate(data, model)
  theta <- fit$par #beta
  
  #cov matrix = inverse hessian
  sigma_theta <- solve(fit$hessian) #sigma_beta
  
  #use filament1_aux_EV to get expectation and variance of model
  aux <- filament1_aux_EV(beta = theta, data = newdata, 
                          model = model, Sigma_beta = sigma_theta)

  pred.df <- data.frame(
    mean = aux$E,
    sd = sqrt(aux$EV + aux$VE),
    lwr = NA,
    upr = NA
  )
  z_score <- qt(1-alpha/2, df = Inf)
  #upper and lower prediction interval
  pred.df$lwr <- pred.df$mean - z_score * pred.df$sd
  pred.df$upr <- pred.df$mean + z_score * pred.df$sd
  
  return (pred.df) 
}

#' Calculate SE and DS Scores
#'
#' @description This function fits a specified model to the given data, makes predictions
#'              on new data, and then calculates the Squared Error (SE) and Dawid-Sebastiani (DS)
#'              scores for the predictions.
#' @param data A `data.frame` used to fit the model.
#' @param model A character string specifying the model to use, either "A" or "B".
#' @param newdata A `data.frame` with the actual weights and other predictors for calculating
#'                the scores.
#' @return A `data.frame` with predicted mean, standard deviation, lower and upper bounds of
#'         the prediction interval, as well as the calculated SE and DS scores for each observation.

fil_score <- function(data, model, newdata){
  #fit the model and get the prediction result
  fit <- filament1_predict(data = data, model = model, newdata = newdata)
  
  #calculate ES and DS score
  act_weight <- newdata$Actual_Weight #actual weight
  mean <- fit$mean
  sd <- fit$sd
  se <- (act_weight - mean)^2
  ds <- (act_weight - mean)^2 / sd^2 + 2*log(sd)
  fit$se <- se
  fit$ds <- ds
  
  return(fit)
}

#' Leave-One-Out Cross-Validation
#'
#' @description Perform leave-one-out cross-validation using the specified model.
#' @param data A `data.frame` containing the data to be used in cross-validation.
#' @param model A character string specifying the model to use, either "A" or "B".
#' @return A `data.frame` with added columns for mean, standard deviation, squared error, and Dawid-Sebastiani score for each observation.

leave1out <- function(data, model){
  n <- nrow(data)
  mean <- numeric(n)
  sd <- numeric(n)
  
  for (i in 1:n){
    training_data <- data[-i, ]
    validation_data <- data[i, ]
    
    fit <- filament1_predict(data = training_data, model = model, newdata = validation_data)
    
    mean[i] <- fit$mean
    sd[i] <- fit$sd
  }
  
  res_df <- data.frame(
    mean = mean, sd = sd
  )
  
  score_res <- cbind(data, res_df) %>%
    mutate(
      se = (Actual_Weight - mean)^2,
      ds = ((Actual_Weight - mean)^2 / sd^2 + 2*log(sd))
    )
  return (score_res)
}


#part 2 codes
#' Calculate log-likelihood for Archaeological Model
#'
#' @description Computes the log-likelihood for a set of observations based on 
#'              the binomial model for archaeological counts.
#' @param df A `data.frame` containing the parameters 'N' and 'phi' for which 
#'           to calculate log-likelihoods.
#' @param y A vector containing the count of left and right femurs, respectively.
#' @return A modified `data.frame` that includes log-likelihoods for each set of parameters.

arch_loglike <- function(df, y) {
  log_likelihoods <- numeric(nrow(df))

  for (i in 1:nrow(df)) {
    N <- df$N[i]
    phi <- df$phi[i]
    #get y1 and y2 from the vector y
    y1 <- y[1]
    y2 <- y[2]

    #calculate the log-likelihood using lgamma for log(gamma)
    log_likelihoods[i] <- -lgamma(y1 + 1) - lgamma(y2 + 1) -
      lgamma(N - y1 + 1) - lgamma(N - y2 + 1) + 2 * lgamma(N + 1) +
      (y1 + y2) * log(phi) + (2 * N - y1 - y2) * log(1 - phi)
  }

  df$log_likelihood <- log_likelihoods
  return(df)
}

#' Estimate Bayesian Posterior via Monte Carlo
#'
#' @description Estimates the posterior probabilities and expectations of 'N' and 'phi' 
#'              using a Monte Carlo integration approach based on a binomial model 
#'              for archaeological counts.
#' @param y A vector containing the count of left and right femurs, respectively.
#' @param xi The probability parameter for the geometric distribution used to sample 'N'.
#' @param a The 'a' parameter for the beta distribution used to sample 'phi'.
#' @param b The 'b' parameter for the beta distribution used to sample 'phi'.
#' @param K The number of Monte Carlo samples to generate.
#' @return A `data.frame` with estimated probabilities and expectations for 'N' and 'phi'.

estimate <- function(y, xi, a, b, K){
  #parameters for df
  N <- rgeom(K, xi)
  phi <- rbeta(K, a, b)
  #df for arch_loglike input containing N and phi
  df <- data.frame(N = N, phi = phi)
  #compute exponential of log likelihood
  exp_arch <- exp(arch_loglike(df, y = y)$log_likelihood)
  py <- (1/K)* sum(exp_arch)
  EN_y <- (1/(K*py)) * sum(N * exp_arch)
  EP_y <- (1/(K*py)) * sum(phi * exp_arch)
  
  df_res <- data.frame(py = py, EN_y = EN_y, EP_y = EP_y)
  return(df_res)
}





