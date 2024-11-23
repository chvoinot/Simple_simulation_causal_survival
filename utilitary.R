library(survival) # Implementation of Kaplan-Meier

library(survRM2)
library(RISCA)
library(grf) # causal_survival_forest function and also 
# survival_forest and probability_forest


library(MASS) # mvrnorm function for simulation
library(rms) # cph function and predict for cph object

library(dplyr)
library(ggplot2)
library(gridExtra) # to plot multiple graph on a page


# Function to calculate the integral of a decreasing function using 
# the rectangle method
# x corresponds to the x coordinate of the function to integrate
# y corresponds to the y coordinate
integral_rectangles <- function(x, y) {
  # Check if the lengths of x and y are the same
  if (length(x) != length(y)) {
    stop("Lengths of x and y must be the same")
  }
  
  # Calculate the width of each rectangle
  dx <- diff(x)
  
  # Initialize the sum
  integral_sum <- 0
  
  # Iterate through each rectangle and sum up the areas
  for (i in 1:(length(x) - 1)) {
    # Calculate the height of the current rectangle
    height <- min(y[i], y[i + 1])
    
    # Multiply the height by the width and add it to the sum
    integral_sum <- integral_sum + height * dx[i]
  }
  mean <- integral_sum + x[1]
  # Return the final integral sum
  return(mean)
}


# Estimate survival function with covariates for each individual at each time Y.grid
# Type of model can be cox or survival forest (n.fold must be completed in this case)

# This function is used also to compute S_c with status = censor.status
estimate_survival_function <- function(data, X.names, 
                                       Y.grid, 
                                       type_of_model = "cox",
                                       T_obs = "T_obs", 
                                       status = "status",
                                       learner = "T-learner",
                                       n.folds = NULL) {
  
  if (learner == "T-learner"){
    # Subset data for A == 0.
    data0 <- data %>% 
      filter(A == 0)
    
    # Subset data for A == 1.
    data1 <- data %>% 
      filter(A == 1)
    # Cox 
    if (type_of_model == "cox") {
      # Formula for cox model (single learner: A as a covariate,
      # T-learner: Stratified fit on A)
      # Here only T-learner
      outcome <- paste0('Surv(', T_obs, ',', status, ')')
      
      # cph do not support notation I(X^2) but X^2 directly (contrary to coxph) 
      X.names <- gsub("I\\((X[0-9]+\\^2)\\)", "\\1", X.names)
      
      
      # Learn Cox regression on two datasets: A|X.
      f <- as.formula(paste(outcome, paste(c(X.names), collapse = " + "), 
                            sep = " ~ "))
      
      # Cox model fitting stratified on A=1
      fitS0 <- cph(f, data = data0, y = TRUE, x = TRUE, times = Y.grid)
      # Cox model fitting stratified on A=0
      fitS1 <- cph(f, data = data1, y = TRUE, x = TRUE, times = Y.grid)
      
      # Predict survival probabilities for each individual at each Y.grid.
      fit.pred1 <- predictCox(fitS1, newdata = data, times = Y.grid, 
                              type = "survival")
      fit.pred0 <- predictCox(fitS0, newdata = data, times = Y.grid, 
                              type = "survival")
      
      # Survival probabilities for each individual at each Y.grid.
      S_hat1 <- fit.pred1$survival
      S_hat0 <- fit.pred0$survival
    } else if (type_of_model == "survival forest") {# Survival forest
      # Initialization
      n <- nrow(data)
      fit.pred1 <- matrix(NA, nrow = n, ncol = length(Y.grid))
      fit.pred0 <- matrix(NA, nrow = n, ncol = length(Y.grid))
      
      if (n.folds > 1) {
        # Split the dataset into n-folds.
        indices <- split(seq(n), sort(seq(n) %% n.folds))
        
        # For each index in each split
        for (idx in indices) {
          # Fit survival forest on observations removed from idx (training set) and A=1
          # A is not included in covariates (T-learner)
          forest.grf1 <- survival_forest(
            X = as.matrix(data[-idx & data[, "A"] == 1, X.names]),
            Y = data[-idx & data[, "A"] == 1, T_obs],
            D = data[-idx & data[, "A"] == 1, status],
            failure.times = Y.grid
          )
          # Fit survival forest on observations removed from idx (training set) and A=0
          # A is not included in covariates (T-learner)
          forest.grf0 <- survival_forest(
            X = as.matrix(data[-idx & data[, "A"] == 0, X.names]),
            Y = data[-idx & data[, "A"] == 0, T_obs],
            D = data[-idx & data[, "A"] == 0, status],
            failure.times = Y.grid
          )
          # Prediction on idx to avoid overfitting
          fit.pred1[idx, ] <- predict(
            forest.grf1, as.matrix(data[idx, X.names]),
            failure.times = Y.grid)$predictions
          
          fit.pred0[idx, ] <- predict(
            forest.grf0, as.matrix(data[idx, X.names]),
            failure.times = Y.grid)$predictions
        }
      } else {# No cross-fitting
        # Fit survival forest on all observations with A=1
        # A is not included in covariates (T-learner)
        forest.grf1 <- survival_forest(
          X = as.matrix(data[data[, "A"] == 1, X.names]),
          Y = data[data[, "A"] == 1, T_obs],
          D = data[data[, "A"] == 1, status],
          failure.times = Y.grid
        )
        # Fit survival forest on all observations with A=0
        # A is not included in covariates (T-learner)
        forest.grf0 <- survival_forest(
          X = as.matrix(data[data[, "A"] == 0, X.names]),
          Y = data[data[, "A"] == 0, T_obs],
          D = data[data[, "A"] == 0, status],
          failure.times = Y.grid
        )
        # Predict on all observations
        fit.pred1 <- predict(forest.grf1, as.matrix(data[, X.names]), 
                             failure.times = Y.grid)$predictions
        fit.pred0 <- predict(forest.grf0, as.matrix(data[, X.names]), 
                             failure.times = Y.grid)$predictions
      }
      S_hat1 <- fit.pred1
      S_hat0 <- fit.pred0
      
    } 
    
  } else if (learner == "S-learner"){
    # Set A=0 for all data
    data0 <- data
    data0$A <- 0
    
    # Set A=1 for all data
    data1 <- data
    data1$A <- 1
    
    if (type_of_model == "cox") {
      outcome <- paste0('Surv(', T_obs, ',', status, ')')
      
      # cph do not support notation I(X^2) but X^2 directly (contrary to coxph) 
      X.names <- gsub("I\\((X[0-9]+\\^2)\\)", "\\1", X.names)
      
      # Learn Cox regression on one datasets and add A as covariate
      f <- as.formula(paste(outcome, paste(c(X.names,"A"), 
                                           collapse = " + "), 
                            sep = " ~ "))
      
      # Fit the two models on the covariates of time Y.grid.
      fitS <- cph(f, data = data, y = TRUE, x = TRUE, times = Y.grid)
      
      # Predict survival probabilities for each individual at each Y.grid.
      fit.pred1 <- predictCox(fitS, newdata = data1, times = Y.grid, 
                              type = "survival")
      fit.pred0 <- predictCox(fitS, newdata = data0, times = Y.grid, 
                              type = "survival")
      
      # Survival probabilities for each individual at each Y.grid.
      S_hat1 <- fit.pred1$survival
      S_hat0 <- fit.pred0$survival
      
    } else if (type_of_model == "survival forest") {
      # Survival forest.
      # Initialize objects
      n <- nrow(data)
      fit.pred1 <- matrix(NA, nrow = n, ncol = length(Y.grid))
      fit.pred0 <- matrix(NA, nrow = n, ncol = length(Y.grid))
      
      if (n.folds > 1) {
        # Split the dataset into n-folds.
        indices <- split(seq(n), sort(seq(n) %% n.folds))
        
        # For all index in each split.
        for (idx in indices) {
          # Fit survival forest on all observations except idx (add A as covariate)
          forest.grf <- survival_forest(
            X = as.matrix(data[-idx, c(X.names,"A")]),
            Y = data[-idx, "T_obs"],
            D = data[-idx, "status"],
            failure.times = Y.grid
          )
          # Predict on idx 
          fit.pred1[idx, ] <- predict(
            forest.grf, as.matrix(data1[idx, c(X.names,"A")]),
            failure.times = Y.grid)$predictions
          
          fit.pred0[idx, ] <- predict(
            forest.grf, as.matrix(data0[idx, c(X.names,"A")]), 
            failure.times = Y.grid)$predictions
        }
      } else {
        # If no cross-fitting 
        # Fit survival forest on all observation (add A as covariate)
        forest.grf <- survival_forest(
          X = as.matrix(data[, c(X.names.outcome,"A")]),
          Y = data[, "T_obs"],
          D = data[, "status"],
          failure.times = Y.grid
        )
        
        # Predict on all observations
        fit.pred1 <- predict(
          forest.grf, as.matrix(data1[, c(X.names.outcome,"A")]), 
          failure.times = Y.grid)$predictions
        fit.pred0 <- predict(
          forest.grf, as.matrix(data0[, c(X.names.outcome,"A")]), 
          failure.times = Y.grid)$predictions
      } 
      S_hat1 <- fit.pred1
      S_hat0 <- fit.pred0
    }
  }
  # Associate the corresponding Survival curve to the observation
  S_hat <- S_hat1 * data$A + (1 - data$A) * S_hat0
  
  return(list('S_hat' = S_hat, "S_hat1" = S_hat1, "S_hat0" = S_hat0, "T" = Y.grid))
}

# Compute the remaining survival function at all time points
Q_t_hat <- function(data, tau, X.names.outcome = c("X1", "X2", "X3", "X4"),
                    nuisance = "cox", n.folds = NULL) {
  # Truncate observed times at tau
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  Y.grid <- sort(unique(data$T_obs_tau))
  
  # Estimate the conditional survival function
  S_hat_all <- estimate_survival_function(
    data = data,
    X.names = X.names.outcome,
    Y.grid = Y.grid,
    type_of_model = nuisance,
    n.folds = n.folds
  )
  
  S.hat <- S_hat_all$S_hat
  
  # Initialize Q.hat matrix
  Y.diff <- diff(c(0, Y.grid))
  Q.hat <- matrix(NA, nrow(S.hat), ncol(S.hat))
  
  # Calculate dot products for conditional expectations
  dot.products <- sweep(S.hat[, 1:(ncol(S.hat) - 1)], 2, Y.diff[2:ncol(S.hat)], "*")
  Q.hat[, 1] <- rowSums(dot.products)
  
  # Update Q.hat backwards to compute conditional expectations
  for (i in 2:(ncol(Q.hat) - 1)) {
    Q.hat[, i] <- Q.hat[, i - 1] - dot.products[, i - 1]
  }
  
  # Normalize by survival probabilities and add back time points
  Q.hat <- Q.hat / S.hat
  Q.hat[is.infinite(Q.hat)] <- 0
  Q.hat <- sweep(Q.hat, 2, Y.grid, "+")
  Q.hat[, ncol(Q.hat)] <- max(Y.grid)
  
  return(Q.hat)
}

# Find the remaining survival function at a specific time y
Q_Y <- function(data, tau, Q.t.hat) {
  # Truncate observed times at tau
  data$T_obs_tau <- ifelse(data$T_obs >= tau, tau, data$T_obs)
  Y.grid <- sort(unique(data$T_obs_tau))
  
  # Find the corresponding Q_t
  Y.index <- findInterval(data$T_obs_tau, Y.grid)
  Q.Y.hat <- Q.t.hat[cbind(seq_along(Y.index), Y.index)]
  
  return(Q.Y.hat)
}


# Function to estimate propensity score
estimate_propensity_score <- function(data, treatment_covariates, 
                                      type_of_model = "glm", n.folds = NULL) {
  # Generalized Linear Model (GLM)
  if (type_of_model == "glm") {
    outcome <- 'A'
    f <- as.formula(paste(outcome, paste(c(treatment_covariates), 
                                         collapse = " + "), sep = " ~ "))
    fitA <- glm(f, data = data, family = binomial(link = "logit"))
    e_hat <- predict(fitA, newdata = data, type = "response")
  }
  
  # Probability Forest (only for continuous variables, 
  # categorical variables need one-hot encoding)
  if (type_of_model == "probability forest" && !is.null(n.folds)) {
    # Initialization
    n <- nrow(data)
    e_hat <- rep(NA, n)
    A <- data$A
    
    # Cross-fitting to avoid overfitting 
    if (n.folds > 1) { 
      # Split the dataset into n folds
      indices <- split(seq(n), sort(seq(n) %% n.folds))
      
      # Learn and predict for each fold
      for (idx in indices) {
        
        # Learn on all data except idx
        propensity_model <- probability_forest(
          as.matrix(data[-idx, treatment_covariates]),
          as.factor(A[-idx]))
        
        # Predict on idx
        e_hat[idx] <- predict(
          propensity_model, 
          newdata = as.matrix(data[idx, treatment_covariates]))$predictions[, 2]
      }
    } 
    # No cross-fitting
    else if (n.folds == 0 | n.folds == 1) {
      
      propensity_model <- probability_forest(
        as.matrix(data[, treatment_covariates]),
        as.factor(A)) 
      
      e_hat <- predict(
        propensity_model,
        newdata = as.matrix(data[, treatment_covariates]))$predictions[, 2]
    }
  }
  return(e_hat)
}


# Compute the area under the survival curve for each individual using the 
# Trapezoidal rule.
# S.hat: predicted survival function for each individual.
expected_survival <- function(S.hat, Y.grid) {
  # Y.grid: vector of time at which to evaluate the survival estimates 
  # (same as S.hat).
  
  # Calculate the distance between each time point.
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
  
  # Compute the area under each survival curve.
  area <- c(base::cbind(1, S.hat) %*% grid.diff)
  
  return(area)
}

# Tool functions 
# Compute hazard function from survival function
estimate_hazard_function <- function(S_hat, Y.grid) {
  
  Y.grid[Y.grid==0]<-0.001
  
  # Calculate differences between successive elements in Y.grid
  Y.diff <- diff(c(0, Y.grid))
  
  # Get the number of columns in S_hat
  grid.length <- ncol(S_hat)
  
  # Compute -log of survival probabilities (cumulative hazard function), 
  # Add 1 as the first value of survival function to ensure that lambda(0)=0
  log.surv.C <- -log(base::cbind(1, S_hat))
  
  # Calculate differences of -log survival probabilities to have 
  # the instantaneous hazard function
  h_hat <- log.surv.C[, 2:(grid.length + 1)] - log.surv.C[, 1:grid.length]
  
  # Divide each column of h_hat by the corresponding element in Y.diff
  h_hat <- sweep(h_hat, 2, Y.diff, "/")
  
  # Return the estimated hazard function
  return(h_hat)
}

integrate <- function(integrand, Y.grid, times) {
  # Create a filter matrix to indicate which elements are within the time 
  # interval
  filter <- sapply(1:length(Y.grid), function(i) {
    return(as.numeric(i <= findInterval(times, Y.grid)))
  })
  
  # Apply the filter to the integrand
  integrand_filtered <- filter * integrand
  
  # Sum the rows of the filtered integrand to get the integrated values
  integrated_value <- rowSums(integrand_filtered)
  
  # Return the integrated values
  return(integrated_value)
}


theta_rmst_survrm2 <- function(data, tau) {
  ATE_pack <- rmst2(data$T_obs, data$status, arm = data$A, tau = tau)
  RMST <- ATE_pack[[5]][1]
  return(RMST)
}

# Function to estimate RMST using unadjusted Kaplan-Meier
RISCA_unadj <- function(data, 
                        tau) {
  # Fit survival curves using Kaplan-Meier stratified by treatment group
  fit <- survfit(Surv(T_obs, status) ~ A, data = data)
  res <- summary(fit)
  
  # Estimate RMST for treatment group A=1
  RMST_A1 <- rmst(
    times = res$time[as.character(res$strata) == "A=1"],
    surv.rates = res$surv[as.character(res$strata) == "A=1"],
    max.time = tau, 
    type = "s" # for step-function
  )
  
  # Estimate RMST for treatment group A=0
  RMST_A0 <- rmst(
    times = res$time[as.character(res$strata) == "A=0"],
    surv.rates = res$surv[as.character(res$strata) == "A=0"],
    max.time = tau, 
    type = "s" # for step-function
  )
  
  # Estimate ATE as the difference in RMST between groups
  ATE_RISCA_unadj <- RMST_A1 - RMST_A0
  return(ATE_RISCA_unadj)
}



# Function to estimate RMST using IPTW Kaplan-Meier
RISCA_iptw <- function(data, 
                       tau, 
                       X.names.propensity, 
                       nuisance_propensity = "glm", 
                       n.folds = NULL) {
  
  # Estimate propensity scores
  e_hat <- estimate_propensity_score(
    data, 
    treatment_covariates = X.names.propensity,
    type_of_model = nuisance_propensity, 
    n.folds = n.folds
  )
  
  # Compute inverse probability weights
  weighted <- (data$A / e_hat) + ((1 - data$A) / (1 - e_hat))
  
  # Fit weighted survival curves
  IPW_pack <- ipw.survival(
    times = data$T_obs, 
    failures = data$status,
    variable = data$A, 
    weights = weighted
  )
  
  # Calculate RMST for treatment group A=1 using weighted survival curve
  RMST_RISCA_A1 <- rmst(
    times = IPW_pack$table.surv$times[IPW_pack$table.surv$variable == 1],
    surv.rates = IPW_pack$table.surv$survival[IPW_pack$table.surv$variable == 1],
    max.time = tau, 
    type = "s"
  )
  
  # Calculate RMST for treatment group A=0 using weighted survival curve
  RMST_RISCA_A0 <- rmst(
    times = IPW_pack$table.surv$times[IPW_pack$table.surv$variable == 0],
    surv.rates = IPW_pack$table.surv$survival[IPW_pack$table.surv$variable == 0],
    max.time = tau, 
    type = "s"
  )
  
  # Compute ATE as the difference in RMST between groups
  ATE_RISCA_IPW <- RMST_RISCA_A1 - RMST_RISCA_A0
  return(ATE_RISCA_IPW)
}


# Function to estimate RMST using single learner G-formula with Cox model
RISCA_gf <- function(data, 
                     tau, 
                     X.names.outcome) {
  
  # Define the outcome formula for the Cox model
  outcome <- paste(c('Surv(', "T_obs", ',', "status", ')'), collapse = "")
  # Single learner : the treatment arm is a predictor
  formula <- as.formula(paste(outcome, paste(c(X.names.outcome, 'A'), 
                                             collapse = " + "), sep = " ~ "))
  
  # Fit the Cox proportional hazards model
  cox.cdt <- coxph(formula, data = data, x = TRUE)
  summary(cox.cdt)
  
  # Compute the effect of the treatment (ATE) using the G-formula
  gc.ate <- gc.survival(
    object = cox.cdt, 
    data = data, 
    group = "A", 
    times = "T_obs",
    failures = "status", 
    max.time = tau, 
    iterations = 100,
    effect = "ATE",
    n.cluster = 1
  )
  
  # Extract the ATE
  ATE_RISCA_gf <- gc.ate$delta[[1]]
  return(ATE_RISCA_gf)
}


# Function to estimate RMST using Causal Survival Random Forest (CSRF)
CSRF <- function(data, X.names, tau) {
  # Fit a causal survival forest
  cf <- causal_survival_forest(X = as.matrix(data[, X.names]), Y = as.matrix(data$T_obs), W = as.matrix(data$A), D = as.matrix(data$status), horizon = tau)
  
  # Predict using the fitted forest
  cf.predict <- predict(cf)
  
  # Estimate the average treatment effect (ATE)
  ATE_csf <- average_treatment_effect(cf)
  
  # Return the estimated ATE
  return(ATE_csf[[1]])
}
