##' The Hirano and Imbens estimator
##'
##' This function estimates the GPS function and estimates the
##' ADRF. The GPS score is based on different treatment models.  The treatment is
##' linearly related to Xs.
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param treat_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.
##' @param outcome_formula is the formula used for fitting the outcome surface. \code{gps} is one of the independent variables to use in the outcome_formula. ie.
##'
##' \code{Y ~ treat+ I(treat^2) + gps + I(gps^2) + treat * gps}
##'
##' or a variation of this.  Use \code{gps} as the name of the variable representing the gps in \code{outcome_formula}.
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param grid_val contains the treatment values to be evaluated.
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model, \code{"Sqrt"} for square-root transformation
##' to a normal treatment, \code{"Poisson"} for Poisson model,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model, \code{"Binomial"} for binomial model.
##' @param link_function For \code{treat_mod = "Gamma"} (fitted using glm) alternatives are "log" or "inverse".
##' For \code{treat_mod = "Binomial"} (fitted using glm) alternatives are "logit", "probit", "cauchit", "log" and "cloglog".
##' @param ... additional arguments to be passed to the outcome lm() function.
##'
##' @return \code{hi_est} returns an object of class "causaldrf",
##' a list that contains the following components:
##' \item{param}{parameter estimates for a hi fit.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{out_mod}{the result of the outcome model fit.}
##' \item{call}{the matched call.}
##'
##' ##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' @export
##'
hi_est <- function (Y,
                    treat,
                    treat_formula,
                    outcome_formula,
                    data,
                    grid_val,
                    treat_mod,
                    link_function,
                    ...)
{
  outcomeData <- data.frame(Y = data[, as.character(Y)],
                            treat = data[, as.character(treat)])
  df <- cbind(data[, which(names(data) %in% Y)],
              data[, which(names(data) %in% treat)],
              data[, -which(names(data) %in% c(Y, treat))])
  colnames(df)[1:2] <- c(Y, treat)
  balanceData <- df[, -1]
  formula_t <- treat_formula
  if (treat_mod == "NegBinom") {
    samp_dat <- balanceData
    result <- MASS::glm.nb(formula = formula_t,
                           link = log,
                           data = samp_dat)
    cond_mean <- result$fitted.values
    cond_var <- cond_mean + cond_mean ^ 2 / result$theta
    prob_nb_est <- (cond_var - cond_mean) / cond_var
    gps_fun_NB <- function(tt) {
      dnbinom(
        x = tt,
        size = result$theta,
        mu = result$fitted.values,
        log = FALSE
      )
    }
    gps_fun <- gps_fun_NB
  }
  else if (treat_mod == "Poisson") {
    samp_dat <- balanceData
    result <- glm(formula = formula_t,
                  family = "poisson",
                  data = samp_dat)
    cond_mean <- result$fitted.values
    samp_dat$gps_vals <-
      dpois(x = outcomeData$treat, lambda = cond_mean)
    gps_fun_Pois <- function(t) {
      dpois(t, lambda = cond_mean)
    }
    gps_fun <- gps_fun_Pois
  }
  else if (treat_mod == "Gamma") {
    samp_dat <- balanceData
    result <-
      glm(
        formula = formula_t,
        family = Gamma(link = link_function),
        data = samp_dat
      )
    est_treat <- result$fitted
    shape_gamma <- as.numeric(MASS::gamma.shape(result)[1])
    theta_given_X <- result$fitted.values / shape_gamma
    theta_treat_X <- outcomeData$treat / shape_gamma
    gps_fun_Gamma <- function(t) {
      dgamma(t, shape = shape_gamma, scale = theta_given_X)
    }
    gps_fun <- gps_fun_Gamma
  }
  else if (treat_mod == "LogNormal") {
    samp_dat <- balanceData
    samp_dat[, as.character(treat)] <-
      log(samp_dat[, as.character(treat)])
    result <- lm(formula = formula_t, data = samp_dat)
    est_log_treat <- result$fitted
    sigma_est <- summary(result)$sigma
    gps_fun_Log <- function(tt) {
      dnorm(log(tt), mean = est_log_treat, sd = sigma_est)
    }
    gps_fun <- gps_fun_Log
  }
  else if (treat_mod == "Sqrt") {
    samp_dat <- balanceData
    samp_dat[, as.character(treat)] <-
      sqrt(samp_dat[, as.character(treat)])
    result <- lm(formula = formula_t, data = samp_dat)
    est_sqrt_treat <- result$fitted
    sigma_est <- summary(result)$sigma
    gps_fun_sqrt <- function(tt) {
      dnorm(sqrt(tt), mean = est_sqrt_treat, sd = sigma_est)
    }
    gps_fun <- gps_fun_sqrt
  }
  else if (treat_mod == "Normal") {
    samp_dat <- balanceData
    result <- lm(formula = formula_t, data = samp_dat)
    gps_fun_Normal <- function(tt) {
      dnorm(tt,
            mean = result$fitted,
            sd = summary(result)$sigma)
    }
    gps_fun <- gps_fun_Normal
  }
  else if (treat_mod == "Binomial") {
    samp_dat <- balanceData
    if (link_function == "logit")
      lf <- binomial(link = logit)
    if (link_function == "probit")
      lf <- binomial(link = probit)
    if (link_function == "cauchit")
      lf <- binomial(link = cauchit)
    if (link_function == "log")
      lf <- binomial(link = log)
    if (link_function == "cloglog")
      lf <- binomial(link = cloglog)
    result <- glm(formula_t, family = lf, data = samp_dat)
    samp_dat$prob_1[outcomeData$treat == 0] <-
      1 - predict.glm(result,
                      type = "response")[outcomeData$treat == 0]
    samp_dat$prob_1[outcomeData$treat == 1] <- predict.glm(result,
                                                           type = "response")[outcomeData$treat == 1]
    gps_fun_Binomial <-
      function(tt)
        ifelse(tt == 1, samp_dat$prob_1,
               1 - samp_dat$prob_1)
    gps_fun <- gps_fun_Binomial
  }
  else {
    stop("No valid treat_mod specified.  Please try again.")
  }
  outcomeData$gps <- gps_fun(outcomeData$treat)
  colnames(outcomeData) <-
    c(as.character(Y),
      as.character(treat),
      "gps")
  outcome_model <- lm(outcome_formula, data = outcomeData, ...)
  outcome_coef <- outcome_model$coef
  mean_outcome_grid <- numeric(length(grid_val))
  for (i in 1:length(grid_val)) {
    temp_matrix <- cbind(numeric(nrow(outcomeData)), grid_val[i],
                         gps_fun(grid_val[i]))
    colnames(temp_matrix) <- c(as.character(Y),
                               as.character(treat),
                               "gps")
    temp_matrix <- data.frame(temp_matrix)
    m_frame <- model.frame(outcome_formula, temp_matrix)
    covar_temp <- model.matrix(outcome_formula, m_frame)
    covar_grid_mat <- covar_temp
    potential_outcomes_at_t <- t(outcome_coef %*% t(covar_grid_mat))
    mean_outcome_grid[i] <- mean(potential_outcomes_at_t)
  }
  z_object <- list(param = mean_outcome_grid,
                   t_mod = result,
                   out_mod = outcome_model,)
  class(z_object) <- "causaldrf"
  z_object
}
