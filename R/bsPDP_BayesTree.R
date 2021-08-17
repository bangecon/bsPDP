##' Bootstrapped partial dependence plots for the \code{BayesTree} estimator for folks who aren't too concerned about doing it quickly. For a faster version, use \code{bsPDP.dbarts}.
##'
##' Evaluates a Bayesian Additive Regression Tree models for \code{n_boot} bootstrap iterations to obtain confidence bands for the estimated PDP.
##'
##' @param formula is an output formula for the outcome.
##' @param variable is the treatment variable.
##' @param data is a data frame to be used for the training.
##' @param newdata is an optional data frame of test data for the PDPs.
##' @param grid sets the values of \code{variable} to evaluate the default is 100 values in range of \code{variable}:
##' grid = seq(min(x), max(x), length.out = 100).
##' @param n_boot is the number of bootstrap replications.
##' @param p_boot is the proportion of the data to select for each bootstrap replication.
##' @param N is the number of observations to select for calculating the PDPs.
##' @param label is a character-string variable label for \code{variable}.
##' @param clock is a logical indicating whether to time each bootstrap replication.
##' @param test is a logical indicating whether to calculate the pdp for both the \code{data} and \code{newdata}.
##' @param seed is a random seed (default is 8675309).
##' @param ... additional arguments for \code{BayesTree}.
##' @return \code{bsPDPbayesTree} returns an object with class "bsPDP," a list that includes the following components:
##' \item{variable}{the treatment variable.}
##' \item{pdpData}{the estimated average predictions and standard errors along \code{variable}.}
##' \item{marginData}{the estimated average marginal effects and standard errors along \code{variable}.}
##' \item{trainData}{the original training data.}
##' \item{testData}{the test data.}
##' \item{trControl}{is not applicable for this method (\code{NULL}).}
##' \item{outcome}{the outcome class (\code{NULL} for outcomes with \code{numeric} class).}
##'
##' @references Yakusheva, Olga; Bang, James T; Bobay, Kathleen; Hughes, Ronda G; Costa, Linda; and Weiss, Marianne (2021, Forthcoming). \emph{Health Services Research.}
##'
##' @export
bsPDP.BayesTree <-
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
           ...) {
    set.seed(seed)
    data <- model.frame(formula, data)
    x <- data[, which(names(data) %in% variable)]
    unregister <- function() {
      env <- foreach:::.foreachGlobals
      rm(list = ls(name = env), pos = env)
    }
    if (is.null(newdata)) {
      newdata <- data
    }
    if (is.null(grid)) {
      grid = seq(min(na.omit(x)), max(na.omit(x)), length.out = 100)
    }
    if (is.null(outcome) & is.factor(data[, 1])) {
      outcome = levels(data[, 1])[2]
    }
    if (is.null(label)) {
      label = variable
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
        BayesTree::bart(x.train = data[, -1],
                     y.train = data[, 1],
                     keeptrees = TRUE,
                     verbose = FALSE,
                     ...)
      if (test == TRUE) {
        newdf.test <-
          newdata[sample(nrow(newdata), size = N, replace = TRUE),-which(names(newdata) %in% variable)]
        xp.test <- cbind(as.data.frame(grid %x% matrix(1, nrow(newdf.test))),
                    newdf[rep(seq_len(nrow(newdf.test)), length(grid)),])
        colnames(xp.test)[1] <- variable
        xp.test$yhat <- colMeans(predict(model, newdata = xp.test))
        xp.test <-
          na.action(cbind(xp.test[, variable], xp.test[, 'yhat']))
        colnames(xp.test) <- c(variable, 'yhat')
        pdps.test <-
          cbind(pdps.test, aggregate(xp.test, list(xp.test[, variable]), mean)$yhat)
        colnames(pdps.test)[1 + i] <- paste0('yhat', i)
      }
      newdf <-
        df[sample(nrow(df), size = N, replace = TRUE), -which(names(df) %in% variable)]
      xp <-
        cbind(as.data.frame(grid %x% matrix(1, nrow(newdf))),
              newdf[rep(seq_len(nrow(newdf)), length(grid)),])
      colnames(xp)[1] <- variable
      xp$yhat <-
        colMeans(predict(model, newdata = xp))
      xp <-
        na.action(cbind(xp[, variable], xp[, 'yhat']))
      colnames(xp) <- c(variable, 'yhat')
      pdps <-
        cbind(pdps, aggregate(xp, list(xp[, variable]), mean)$yhat)
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
      pdpHat.test <-
        matrixStats::rowMeans2(pdps.test[, 2:(1 + n_boot)])
      pdpSd.test <-
        matrixStats::rowSds(pdps.test[, 2:(1 + n_boot)])
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
      trainMethod = 'BayesTree'
    )
    class(out) <- c('bsPDP', class(out))
    return(out)
  }
