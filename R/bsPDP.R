##' Bootstrapped partial dependence plots
##'
##' Evaluates models for \code{n_boot} bootstrap iterations to obtain confidence bands for the estimated PDP.
##'
##' @details \code{bsPDP} currently supports the following methods:
##'
##' \code{lm} for linear models
##'
##' \code{glm} for generalized linear models
##'
##' \code{randomForest} for random forests
##'
##' \code{caret} for pramater tuning with machine-learning models (including tree ensembles and neural networks)
##'
##' \code{dbarts} and \code{BayesTree} for Bayesian Additive Regression Trees (BART)
##'
##' \code{causaldrf} for inverse probability-of-treatment weights estimators and the Hirano-Imbens (2004) covariate balancing estimators.
##'
##' @param formula (required) an output formula for the outcome.
##' @param variable (required) the treatment variable.
##' @param data (required) a data frame to be used for the training.
##' @param newdata an optional data frame of test data for the PDPs (default is \code{data}).
##' @param grid sets the values of \code{variable} to evaluate (the default is 100 values in range of \code{variable}).
##' @param n_boot the number of bootstrap replications (default is 100).
##' @param p_boot the proportion of the data to select for each bootstrap replication (default is 0.6).
##' @param N the number of observations to select for calculating the PDPs (default is 1000).
##' @param method the model method (default is \code{'randomForest'}.
##' @param label a character-string variable label for \code{variable} (default is \code{variable}).
##' @param clock a logical indicating whether to time each bootstrap replication (default is \code{FALSE}.
##' @param test a logical indicating whether to calculate the pdp for both the \code{data} and \code{newdata} (default is \code{FALSE}.
##' @param plot a logical indicating whether to plot the resulting \code{bsPDP} object (default is \code{TRUE}.
##' @param seed a random seed (default is 8675309).
##' @param ... additional arguments specific to the method some additional arguments may be required for some methods, e.g. \code{caret}, \code{iptw_est}, and \code{hi_est}.
##'
##' @return \code{bsPDP} returns an object with class "bsPDP," a list that includes the following components:
##' \item{variable}{the treatment variable.}
##' \item{pdpData}{the estimated average predictions and standard errors along \code{variable}.}
##' \item{trainData}{the original training data.}
##' \item{testData}{the test data.}
##' \item{outcome}{the outcome class (\code{NULL} for outcomes with \code{numeric} class).}
##' \item{method}{the method}
##'
##' @references Yakusheva, Olga; Bang, James T; Bobay, Kathleen; Hughes, Ronda G; Costa, Linda; and Weiss, Marianne (2021, Forthcoming). \emph{Health Services Research.}
##'
##' Breiman, Leo, 2001. Random forests. \emph{Machine learning.} 45(1), pp.5-32.
##'
##' Chipman, H.; George, E.; and McCulloch, R. (2010) BART: Bayesian Additive Regression Trees. \emph{Annals of Statistics.} 4(1), pp. 266-298.
##'
##' Hirano, Keisuke, Imbens, Guido W (2004).  The propensity score with continuous treatments.  \emph{Applied Bayesian modeling and causal inference from incomplete-data perspectives.} 226164, pp. 73-84.
##'
##' Schafer, J.L., Galagate, D.L. (2015).  Causal inference with a continuous treatment and outcome: alternative estimators for parametric dose-response models. \emph{Doctoral dissertation}.
##'
##' @export
##'
bsPDP <-
  function(formula,
           variable = NULL,
           data,
           newdata = NULL,
           grid = NULL,
           outcome = NULL,
           n_boot = 100,
           p_boot = 0.6,
           N = 1000,
           method = 'randomForest',
           label = NULL,
           clock = FALSE,
           test = FALSE,
           seed = 8675309,
           ...) {
    if (method == 'lm') {
      out <- bsPDP.lm(
        formula = formula,
        variable = variable,
        data = data,
        newdata = newdata,
        grid = grid,
        outcome = outcome,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        label = label,
        clock = clock,
        test = test,
        seed = seed,
        ...
      )
    }
    if (method == 'glm') {
      out <- bsPDP.glm(
        formula = formula,
        variable = variable,
        data = data,
        newdata = newdata,
        grid = grid,
        outcome = outcome,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        label = label,
        clock = clock,
        test = test,
        seed = seed,
        ...
      )
    }
    if (method == 'randomForest') {
      out <- bsPDP.randomForest(
        formula = formula,
        variable = variable,
        data = data,
        newdata = newdata,
        grid = grid,
        outcome = outcome,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        label = label,
        clock = clock,
        test = test,
        seed = seed,
        ...
      )
    }
    if (method == 'caret') {
      out <- bsPDP.caret(
        formula = formula,
        variable = variable,
        data = data,
        newdata = newdata,
        grid = grid,
        outcome = outcome,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        label = label,
        clock = clock,
        test = FALSE,
        seed = 8675309,
        ...
      )
    }
    if (method == 'dbarts') {
      out <- bsPDP.dbarts(
        formula = formula,
        variable = variable,
        data = data,
        newdata = newdata,
        grid = grid,
        outcome = outcome,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        label = label,
        clock = clock,
        test = test,
        seed = seed,
        ...
      )
    }
    if (method == 'BayesTree') {
      out <- bsPDP.BayesTree(
        formula = formula,
        variable = variable,
        data = data,
        newdata = newdata,
        grid = grid,
        outcome = outcome,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        label = label,
        clock = clock,
        test = test,
        seed = seed,
        ...
      )
    }
    if (method == 'iptw_est') {
      out <- bsPDP.iptw_est(
        formula = formula,
        variable = variable,
        data = data,
        newdata = newdata,
        grid = grid,
        outcome = outcome,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        label = label,
        clock = clock,
        test = test,
        seed = seed,
        ...
      )
    }
    if (method == 'hi_est') {
      out <- bsPDP.hi_est(
        formula = formula,
        data,
        newdata = newdata,
        outcome = outcome,
        variable = variable,
        label = label,
        grid = grid,
        clock = clock,
        test = test,
        n_boot = n_boot,
        p_boot = p_boot,
        N = N,
        seed = seed,
        ...
      )
    }
    if (!method %in% c('lm',
                       'glm',
                       'randomForest',
                       'caret',
                       'BayesTree',
                       'dbarts',
                       'iptw_est',
                       'hi_est')) {
      stop("Method not currently supported for bsPDP.")
    }
    class(out) <- 'bsPDP'
    return(out)
  }
