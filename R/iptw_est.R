##' The iptw method or importance weighting method estimates the ADRF by
##' weighting the data with stabilized or non-stabilized weights.
##'
##' @param Y is the the name of the outcome variable contained in \code{data}.
##' @param treat is the name of the treatment variable contained in
##' \code{data}.
##' @param treat_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.
##' @param numerator_formula an object of class "formula" (or one that can be
##' coerced to that class) that regresses \code{treat} on a linear combination
##' of \code{X}: a symbolic description of the model to be fitted.  i.e.
##' \code{treat ~ 1}.
##' @param data is a dataframe containing \code{Y}, \code{treat}, and
##' \code{X}.
##' @param degree is 1 for linear and 2 for quadratic outcome model.
##' @param treat_mod a description of the error distribution to be used in the
##' model for treatment. Options include: \code{"Normal"} for normal model,
##' \code{"LogNormal"} for lognormal model, \code{"Sqrt"} for square-root transformation
##' to a normal treatment, \code{"Poisson"} for Poisson model,
##' \code{"NegBinom"} for negative binomial model, \code{"Gamma"} for gamma
##' model, \code{"Binomial"} for binomial model, \code{"Ordinal"} for ordinal model,
##' \code{"Multinomial"} for multinomial model.
##' @param link_function specifies the link function between the variables in
##' numerator or denominator and exposure, respectively.
##' For \code{treat_mod = "Gamma"} (fitted using glm) alternatives are "log" or "inverse".
##' For \code{treat_mod = "Binomial"} (fitted using glm) alternatives are "logit", "probit", "cauchit", "log" and "cloglog".
##' For \code{treat_mod = "Multinomial"} this argument is ignored, and
##' multinomial logistic regression models are always used (fitted using multinom).
##' For \code{treat_mod = "Ordinal"} (fitted using polr) alternatives are "logit", "probit", "cauchit", and "cloglog".
##' @param ... additional arguments to be passed to the low level treatment regression fitting functions.
##'
##' @return \code{iptw_est} returns an object of class "causaldrf", a list that contains the following components:
##' \item{param}{parameter estimates for a iptw fit.}
##' \item{t_mod}{the result of the treatment model fit.}
##' \item{num_mod}{the result of the numerator model fit.}
##' \item{weights}{the estimated weights.}
##' \item{weight_data}{the weights.}
##' \item{out_mod}{the outcome model.}
##'
##' @references Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a
##' continuous treatment and outcome: alternative estimators for parametric
##' dose-response models. \emph{Manuscript in preparation}.
##'
##' @export
##'
iptw_est <-
  function (Y,
            treat,
            treat_formula,
            numerator_formula,
            data,
            degree,
            treat_mod,
            link_function,
            ...)
  {
    outcomeData <- data.frame(Y = data[, as.character(Y)],
                              treat = data[, as.character(treat)])
    df <- cbind(data[, which(names(data) %in% Y)],
                data[, which(names(data) %in% treat)],
                data[,-which(names(data) %in% c(Y, treat))])
    colnames(df)[1:2] <- c(Y, treat)
    balanceData <- df[,-1]
    formula_t <- treat_formula
    if (treat_mod == "Poisson") {
      samp <- balanceData
      result2 <- stats::glm(formula_t, family = "poisson",
                            data = samp, ...)
      cond_mean <- result2$fitted.values
      samp$gps_vals <-
        dpois(x = outcomeData$treat, lambda = cond_mean)
      result1 <- glm(numerator_formula,
                     family = "poisson",
                     data = samp,
                     ...)
      cond_mean_num <- result1$fitted.values
      samp$num_weight <-
        dpois(x = outcomeData$treat, lambda = cond_mean_num)
      est_import_wt <- samp$num_weight / samp$gps_vals
    }
    else if (treat_mod == "NegBinom") {
      samp <- balanceData
      result2 <- MASS::glm.nb(formula_t, link = log, data = samp,
                              ...)
      cond_mean <- result2$fitted.values
      cond_var <- cond_mean + cond_mean ^ 2 / result2$theta
      prob_nb_est <- (cond_var - cond_mean) / cond_var
      samp$gps_vals <-
        dnbinom(
          x = outcomeData$treat,
          size = result2$theta,
          mu = result2$fitted.values,
          log = FALSE
        )
      result1 <- MASS::glm.nb(numerator_formula, link = link_function,
                              data = samp, ...)
      cond_mean_num <- result1$fitted.values
      cond_var_num <- cond_mean_num + cond_mean_num ^ 2 / result1$theta
      prob_nb_est_num <- (cond_var_num - cond_mean_num) / cond_var_num
      samp$num_weight <-
        dnbinom(
          x = outcomeData$treat,
          size = result1$theta,
          prob = prob_nb_est_num,
          mu = result1$fitted.values,
          log = FALSE
        )
      est_import_wt <- samp$num_weight / samp$gps_vals
    }
    else if (treat_mod == "Gamma") {
      samp <- balanceData
      result2 <- glm(formula_t,
                     family = Gamma(link = link_function),
                     data = samp,
                     ...)
      shape_gamma <- as.numeric(MASS::gamma.shape(result2)[1])
      theta_given_X <- result2$fitted.values / shape_gamma
      samp$gps_vals <- dgamma(outcomeData$treat, shape = shape_gamma,
                              scale = theta_given_X)
      result1 <- glm(numerator_formula,
                     family = Gamma(link = "log"),
                     data = samp,
                     ...)
      shape_gamma_2 <- as.numeric(MASS::gamma.shape(result1)[1])
      theta_given_X_2 <- result1$fitted.values / shape_gamma
      samp$num_weight <-
        dgamma(outcomeData$treat, shape = shape_gamma_2,
               scale = theta_given_X_2)
      est_import_wt <- samp$num_weight / samp$gps_vals
    }
    else if (treat_mod == "LogNormal") {
      samp <- balanceData
      samp[, as.character(treat)] <- log(samp[, as.character(treat)])
      result2 <- lm(formula_t, data = samp, ...)
      samp$gps_vals <- dnorm(samp[, as.character(treat)],
                             mean = result2$fitted,
                             sd = summary(result2)$sigma)
      result1 <- lm(numerator_formula, data = samp, ...)
      samp$num_weight <- dnorm(samp[, as.character(treat)],
                               mean = result1$fitted,
                               sd = summary(result1)$sigma)
      est_import_wt <- samp$num_weight / samp$gps_vals
    }
    else if (treat_mod == "Sqrt") {
      samp <- balanceData
      samp[, as.character(treat)] <- sqrt(samp[,
                                               as.character(treat)])
      result2 <- lm(formula_t, data = samp, ...)
      samp$gps_vals <- dnorm(samp[, as.character(treat)],
                             mean = result2$fitted,
                             sd = summary(result2)$sigma)
      result1 <- lm(numerator_formula, data = samp, ...)
      samp$num_weight <- dnorm(samp[, as.character(treat)],
                               mean = result1$fitted,
                               sd = summary(result1)$sigma)
      est_import_wt <- samp$num_weight / samp$gps_vals
    }
    else if (treat_mod == "Normal") {
      samp <- balanceData
      result2 <- lm(formula_t, data = samp, ...)
      samp$gps_vals <- dnorm(outcomeData$treat,
                             mean = result2$fitted,
                             sd = summary(result2)$sigma)
      result1 <- lm(numerator_formula, data = samp, ...)
      samp$num_weight <-
        dnorm(outcomeData$treat,
              mean = result1$fitted,
              sd = summary(result1)$sigma)
      est_import_wt <- samp$num_weight / samp$gps_vals
    }
    else if (treat_mod == "Binomial") {
      samp <- balanceData
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
      if (is.null(numerator_formula))
        samp$num_weight <- 1
      else {
        result1 <- glm(numerator_formula,
                       family = lf,
                       data = samp,
                       ...)
        samp$num_weight[outcomeData$treat == 0] <-
          1 - predict.glm(result1,
                          type = "response")[outcomeData$treat == 0]
        samp$num_weight[outcomeData$treat == 1] <-
          predict.glm(result1,
                      type = "response")[outcomeData$treat == 1]
        result1$call$formula <- eval(numerator_formula)
        result1$call$family <- link_function
        result1$call$data <- eval(data)
      }
      result2 <- glm(formula_t, family = lf, data = samp,
                     ...)
      samp$denom_weight[outcomeData$treat == 0] <-
        1 - predict.glm(result2,
                        type = "response")[outcomeData$treat == 0]
      samp$denom_weight[outcomeData$treat == 1] <-
        predict.glm(result2,
                    type = "response")[outcomeData$treat == 1]
      result2$call$formula <- eval(formula_t)
      result2$call$family <- link_function
      result2$call$data <- eval(data)
      est_import_wt <- samp$num_weight / samp$denom_weight
    }
    else if (treat_mod == "Multinomial") {
      samp <- balanceData
      if (is.null(numerator_formula))
        samp$num_weight <- 1
      else {
        result1 <- nnet::multinom(formula = numerator_formula,
                                  data = samp, ...)
        pred1 <- as.data.frame(predict(result1, type = "probs"))
        samp$num_weight <- vector("numeric", nrow(outcomeData))
        for (i in 1:length(unique(outcomeData$treat)))
          samp$num_weight[with(outcomeData,
                               treat == sort(unique(outcomeData$treat))[i])] <-
          pred1[outcomeData$treat ==
                  sort(unique(outcomeData$treat))[i], i]
        result1$call$formula <- eval(numerator_formula)
        result1$call$data <- eval(data)
      }
      result2 <- nnet::multinom(formula = formula_t, data = samp,
                                ...)
      pred2 <- as.data.frame(predict(result2, type = "probs"))
      samp$denom_weight <- vector("numeric", nrow(outcomeData))
      for (i in 1:length(unique(outcomeData$treat)))
        samp$denom_weight[with(outcomeData,
                               treat == sort(unique(outcomeData$treat))[i])] <-
        pred2[outcomeData$treat ==
                sort(unique(outcomeData$treat))[i], i]
      result2$call$formula <- eval(formula_t)
      result2$call$data <- eval(data)
      est_import_wt <- samp$num_weight / samp$denom_weight
    }
    else if (treat_mod == "Ordinal") {
      if (link_function == "logit")
        m <- "logistic"
      if (link_function == "probit")
        m <- "probit"
      if (link_function == "cloglog")
        m <- "cloglog"
      if (link_function == "cauchit")
        m <- "cauchit"
      if (is.null(numerator_formula))
        samp$num_weight <- 1
      else {
        result1 <- MASS::polr(
          formula = eval(numerator_formula),
          data = samp,
          method = m,
          ...
        )
        pred1 <- as.data.frame(predict(result1, type = "probs"))
        samp$num_weight <- vector("numeric", nrow(outcomeData))
        for (i in 1:length(unique(outcomeData$treat)))
          samp$num_weight[with(outcomeData,
                               treat == sort(unique(outcomeData$treat))[i])] <-
          pred1[outcomeData$treat ==
                  sort(unique(outcomeData$treat))[i], i]
        result1$call$formula <- eval(numerator_formula)
        result1$call$data <- eval(data)
        result1$call$method <- m
      }
      result2 <- MASS::polr(
        formula = eval(formula_t),
        data = samp,
        method = m,
        na.action = na.fail,
        ...
      )
      pred2 <- as.data.frame(predict(result2, type = "probs"))
      samp$denom_weight <- vector("numeric", nrow(outcomeData))
      for (i in 1:length(unique(outcomeData$treat)))
        samp$denom_weight[with(outcomeData,
                               treat == sort(unique(outcomeData$treat))[i])] <-
        pred2[outcomeData$treat ==
                sort(unique(outcomeData$treat))[i], i]
      result2$call$formula <- eval(formula_t)
      result2$call$data <- eval(data)
      result2$call$method <- m
      est_import_wt <- samp$num_weight / samp$denom_weight
    }
    else {
      stop("Treatment model specified is not valid.  Please try again.")
    }
    if (degree == 2) {
      weight_dat <-
        data.frame(outcomeData$Y,
                   outcomeData$treat,
                   outcomeData$treat ^ 2,
                   est_import_wt)
      colnames(weight_dat) <-
        c("Y", "treat", "treat.2", "est_import_wt")
      design.weight <-
        survey::svydesign(ids = ~ 1,
                          weights = ~ est_import_wt,
                          data = weight_dat)
      import_weight_mod <- survey::svyglm(Y ~ treat + treat.2,
                                          design = design.weight)
      se_est_iptw_coef <- sqrt(diag(vcov(import_weight_mod)))
      import_coefs <- coef(import_weight_mod)
    }
    else if (degree == 1) {
      weight_dat <-
        data.frame(outcomeData$Y, outcomeData$treat, est_import_wt)
      colnames(weight_dat) <- c("Y", "treat", "est_import_wt")
      design.weight <-
        survey::svydesign(ids = ~ 1,
                          weights = ~ est_import_wt,
                          data = weight_dat)
      import_weight_mod <-
        survey::svyglm(Y ~ treat, design = design.weight)
      se_est_iptw_coef <- sqrt(diag(vcov(import_weight_mod)))
      import_coefs <- coef(import_weight_mod)
    }
    else {
      stop("Error: degree needs to be 1 or 2")
    }
    z_object <- list(
      param = import_coefs,
      t_mod = result2,
      num_mod = result1,
      weights = est_import_wt,
      weight_data = weight_dat,
      out_mod = import_weight_mod
    )
    class(z_object) <- "causaldrf"
    z_object
  }
