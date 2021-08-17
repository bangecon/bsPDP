##' Bootstrapped partial dependence plots for the \code{iptw_est} estimator
##'
##' Evaluates a inverse probability of treatment estimator for \code{n_boot} bootstrap iterations to obtain confidence bands for the estimated PDP.
##'
##' @param formula is the outcome formula.
##' @param variable is the treatment variable.
##' @param data is a data frame to be used for the training.
##' @param newdata is an optional data frame of test data for the PDPs.
##' @param grid sets the values of \code{variable} to evaluate the default is 100 values in range of \code{variable}:
##' grid = seq(min(x), max(x), length.out = 100).
##' @param outcome is the outcome class to be predicted for classification problems.
##' @param n_boot is the number of bootstrap replications.
##' @param p_boot is the proportion of the data to select for each bootstrap replication.
##' @param N is the number of observations to select for calculating the PDPs.
##' @param method indicates the type of model (\code{causaldrf}).
##' @param label is a character-string variable label for \code{variable}.
##' @param clock is a logical indicating whether to time each bootstrap replication.
##' @param test is a logical indicating whether to calculate the pdp for both the \code{data} and \code{newdata}.
##' @param seed is a random seed (default is 8675309).
##' @param treatment_formula is a formula for covariate-balancing the treatment (denominator of the IPTW).
##' @param treatment_mod is the type of model for the treatment (default is "Binomial" for factor outcomes or "Normal" otherwise).
##' @param numerator_formula is a formula for the numerator of the inverse probability weight.
##' @param link indicates the type of link function (default is "logit" for \code{treatment_mod = "Binomial"} or \code{treatment_mod = "Ordinal"}).
##' @param degree is the polynomial order for the outcome equation (equal to 1 or 2; default = 1).
##' @param ... additional arguments for the \code{iptw_est} function.
##' @return \code{bsPDPglm} returns an object with class "bsPDP," a list that includes the following components:
##' \item{variable}{the treatment variable.}
##' \item{pdpData}{the estimated average predictions and standard errors along \code{variable}.}
##' \item{marginData}{the estimated average marginal effects and standard errors along \code{variable}.}
##' \item{trainData}{the original training data.}
##' \item{testData}{the test data.}
##' \item{outcome}{the outcome class (\code{NULL} for outcomes with \code{numeric} class).}
##' \item{trControl}{is not applicable for this method (\code{NULL}).}
##' \item{method}{the method (\code{lm})}
##' @export
bsPDP.iptw_est <-
  function(formula,
           variable = NULL,
           data,
           newdata = NULL,
           grid = NULL,
           outcome = NULL,
           n_boot = 100,
           p_boot = 0.6,
           N = 1000,
           label = NULL,
           clock = FALSE,
           test = FALSE,
           seed = 8675309,
           treatment_formula = NULL,
           treatment_mod = NULL,
           numerator_formula = NULL,
           link = NULL,
           degree = 1,
           ...) {
    set.seed(seed)
    data <- model.frame(formula, data)
    x <- data[, which(names(data) %in% variable)]
    Y <- data[, 1]
    if (is.null(newdata)) {
      newdata <- data
    }
    if (is.null(grid)) {
      grid = seq(min(na.omit(x)), max(na.omit(x)), length.out = 100)
    }
    if (is.null(label)) {
      label = variable
    }
    if (is.null(formula)) {
      outcome_formula <-
        formula(Y ~ x + I(x ^ 2) + gps + I(gps ^ 2) + x * gps)
    }
    if (is.null(treatment_formula)) {
      treatment_formula <- formula(Y ~ 1)
      cat("Treatment formula is NULL; estmiating unweighted treatment effects")
    }
    if (is.null(treatment_mod)) {
      if (is.factor(x)) {
        treatment_mod <- "Binomial"
      } else {
        treatment_mod <- "Normal"
      }
    }
    if (is.null(link)) {
      if (treatment_mod == "Gamma") {
        link = "log"
      }
      if (treatment_mod == "Binomial" |
          treatment_mod == "Ordinal") {
        link = "logit"
      }
      if (treatment_mod == "Normal") {
        link = "NULL"
      }
    }
    if (is.null(outcome) & is.factor(data[, 1])) {
      outcome = levels(data[, 1])[2]
    }
    if (is.factor(data[, 1])) {
      data[, 1] <- as.numeric(ifelse(data[, 1] == outcome, 1, 0))
    }
    model <- NULL
    pdps  <- as.matrix(grid)
    pdps.test  <- as.matrix(grid)
    if (clock == TRUE) {
      cat(paste0(
        "Began bootstrap replications at ",
        noquote(strftime(Sys.time(), "%H:%M:%S")),
        ". \n"
      ))
    }
    for (i in 1:n_boot) {
      df <-
        data[sample(nrow(data), size = ceiling(p_boot * nrow(data))), ]
      model <-
        iptw_est(
          Y = Y,
          treat = variable,
          treat_formula = formula(treatment_formula),
          outcome_formula = formula(formula),
          data = df,
          degree = degree,
          treat_mod = treatment_mod,
          link_function = link,
          ...
        )
      if (test == TRUE) {
        newdf.test <-
          newdata[sample(nrow(newdata), size = N, replace = TRUE), -which(names(newdata) %in% variable)]
        xp.test <- as.data.frame(cbind(rep(1, length(grid)), grid))
        colnames(xp.test)[1:2] <- c("Intercept", variable)
        if (degree == 2) {
          xp.test <- cbind(xp.test, grid ^ 2)
          colnames(xp.test)[3] <- paste0(variable, ".2")
        }
        xp.test$yhat <-
          as.matrix(xp.test[, 1:(degree + 1)]) %*%
          as.vector(model$out_mod$coefficients)
        xp.test <-
          as.data.frame(cbind(xp.test[, variable], xp.test[, 'yhat']))
        colnames(xp.test) <- c(variable, 'yhat')
        pdps.test <-
          cbind(pdps.test, aggregate(xp.test, list(xp.test[, 1]), mean)$yhat)
        colnames(pdps.test)[1 + i] <- paste0('yhat', i)
      }
      newdf <-
        df[sample(nrow(df), size = N, replace = TRUE), -which(names(df) %in% variable)]
      xp <- as.data.frame(cbind(rep(1, length(grid)), grid))
      colnames(xp)[1:2] <- c("Intercept", variable)
      if (degree == 2) {
        xp <- cbind(xp, grid ^ 2)
        colnames(xp)[3] <- paste0(variable, ".2")
      }
      xp$yhat <-
        as.matrix(xp[, 1:(degree + 1)]) %*%
        as.vector(model$out_mod$coefficients)
      xp <-
        as.data.frame(cbind(xp[, variable], xp[, 'yhat']))
      colnames(xp) <- c(variable, 'yhat')
      pdps <-
        cbind(pdps, aggregate(xp, list(xp[, 1]), mean)$yhat)
      colnames(pdps)[1 + i] <- paste0('yhat', i)
      rm(df, newdf, xp, newdf.test, xp.test, model)
      if (clock == TRUE) {
        cat(paste0(
          "Completed replication ",
          i,
          " out of ",
          n_boot,
          " at ",
          noquote(strftime(Sys.time(), "%H:%M:%S")),
          ". \n"
        ))
      }
    }
    if (test == TRUE) {
      for (i in 1:length(grid)) {
        if (i == 1) {
          margins.test <-
            matrix(c(grid[i], rep(NA, n_boot)), nrow = 1, ncol = n_boot + 1)
        }
        if (i == length(grid)) {
          margins.test <- rbind(margins.test, c(grid[i], rep(NA, n_boot)))
        }
        if (i < length(grid) & i > 1) {
          margins.test <-
            rbind(margins.test, c(grid[i],
                                  (pdps.test[i + 1,] - pdps.test[i - 1,]) / (grid[i + 1] - grid[i - 1])))
        }
      }
      pdpHat.test <- matrixStats::rowMeans2(pdps.test[, 2:(1 + n_boot)])
      pdpSd.test <- matrixStats::rowSds(pdps.test[, 2:(1 + n_boot)])
    } else{
      pdpHat.test <- c(rep(NA, length(grid)))
      pdpSd.test <- c(rep(NA, length(grid)))
    }
    for (i in 1:length(grid)) {
      if (i == 1) {
        margins <-
          matrix(c(grid[i], rep(NA, n_boot)), nrow = 1, ncol = n_boot + 1)
      }
      if (i == length(grid)) {
        margins <- rbind(margins, c(grid[i], rep(NA, n_boot)))
      }
      if (i < length(grid) & i > 1) {
        margins <-
          rbind(margins, c(grid[i],
                           (pdps[i + 1,] - pdps[i - 1,]) / (grid[i + 1] - grid[i - 1])))
      }
    }
    pdpHat <- matrixStats::rowMeans2(pdps[, 2:(1 + n_boot)])
    pdpSd <- matrixStats::rowSds(pdps[, 2:(1 + n_boot)])
    pdpData = cbind(grid, pdpHat, pdpSd, pdpHat.test, pdpSd.test)
    colnames(pdpData)[1] <- variable
    out <- list(
      trainData = data,
      variable = variable,
      pdpData = pdpData,
      marginData = marginData,
      testData = newdata,
      outcome = c(names(data)[1], outcome),
      trControl = NULL,
      trainMethod = 'iptw_est'
    )
    class(out) <- c('bsPDP', class(out))
    return(out)
  }
