##' Plot bootstrapped partial dependence plots
##'
##' Plots objects of class \code{bsPDP} created by bootstrapping a predictive or causal model to show confidence bands for the average predicted partial dependence plot of a variable.
##'
##' @details \code{bsPDP} and \code{plot.bsPDP} currently support the following methods:
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
##' \code{causaldrf} for inverse probability-of-treatment weights estimators \code{iptw_est} and the Hirano-Imbens (2004) covariate balancing estimators \code{hi_est}.
##'
##' @param bsPDP (required) a list created by a \code{bsPDP} function and contains the following elements:
##' 1. \code{trainData}: a data frame containing the original training sample.
##' 2. \code{variable}: the x-variable from \code{trainData} to predict and plot.
##' 3. \code{pdpData}: data frame containing the following variables:
##'    a. \code{variable}: grid values of the x-variable;
##'    b. \code{pdpHat}: the average of the bootstrapped predictions at each grid value in the training sample;
##'    c. \code{pdpSd}: the standard errors of the bootstrapped predictions at each grid value in the training sample;
##'    d. \code{pdpHat.test}: the average of the bootstrapped predictions at each grid value in the test sample;
##'    e. \code{pdpSd.test}: the standard errors of the bootstrapped predictions at each grid in the test sample;
##' 4. \code{testData}: a data frame containing the test sample.
##' 5. \code{outcome}: the outcome class (for classification).
##' 6. \code{trControl}: the \code{trainControl} settings (for tuned \code{caret} methods).
##' 7. \code{trainMethod}: the \code{method} used in the \code{caret} functionl.
##' @param type is one of \code{"response"} or \code{"margins"} indicating whether to plot the predicted outcome (\code{"response"}, the default) or the approximate average marginal effect along \code{variable}
##' @param test is a logical indicating whether to plot the pdp for the test data instead of the train data (default is \code{FALSE}).
##' @param title is a string for the main title of the plot (the default is "Partial Dependence Plot").
##' @param xlim is a vector combination equal to the minimum and maximum values of the x-axis (the default is the minimum and maximum values of the x-variable).
##' @param ylim is a vector combination equal to the minimum and maximum values of the y-axis (the default is the minimum of \code{pdpHat} minus 3 times the maximum value of \code{pdpSd} and the maximum of pdpHat plus 3 times the maximum value of \code{pdpSd}).
##' @param xtitle is a character string for the x-axis title (the default is the name of the x-variable)
##' @param ytitle is a character string for the y-axis title (the default is "Predicted Outcome")
##' @param iqr is a logical indicating whether to include lines at the IQR bounds of x (the default is \code{FALSE})
##' @param hist is a logical indicating whether to include a histogram of x (the default is \code{FALSE})
##'
##' @return \code{plot.bsPDP} returns a \code{ggplot} that can bee manipulated according to \code{ggplot} syntax.
##'
##' @references Yakusheva, Olga; Bang, James T; Bobay, Kathleen; Hughes, Ronda G; Costa, Linda; and Weiss, Marianne (2021, Forthcoming). \emph{Health Services Research.}
##'
##' @export
##'
plot.bsPDP <-
  function(bsPDP,
           type = "response",
           test = FALSE,
           title = "Partial Dependence Plot",
           xlim = NULL,
           ylim = NULL,
           xtitle = NULL,
           ytitle = "Predicted Outcome",
           iqr = FALSE,
           hist = FALSE) {
    # Set some plot parameters
    trainData <- bsPDP$trainData
    if (type == "response") {
      data <- bsPDP$pdpData
      if (test == TRUE) {
        yhat <- 'pdpHat.test'
        ysd  <- 'pdpSd.test'
      } else {
        yhat <- 'pdpHat'
        ysd  <- 'pdpSd'
      }
    } else {
      data <- bsPDP$marginData
      if (test == TRUE) {
        yhat <- 'marginHat.test'
        ysd  <- 'marginSd.test'
      } else {
        yhat <- 'marginHat'
        ysd  <- 'marginSd'
      }
    }
    pdpData <- bsPDP$pdpData
    if (is.null(xlim)) {
      xlim <-
        quantile(pdpData[, 1], probs = c(0, 1))
    }
    if (is.null(ylim)) {
      ylim <-
        c(min(data[, which(colnames(data) %in% yhat)]) - 3 *
            max(data[, which(colnames(data) %in% ysd)]),
          max(data[, which(colnames(data) %in% yhat)]) + 3 *
            max(data[, which(colnames(data) %in% ysd)]))
    }
    if (hist == TRUE) {
      ylim[1] = min(0, ylim[1])
    }
    if (is.null(xtitle)) {
      xtitle <- names(pdpData)[1]
    }
    out <- ggplot() +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(limits = xlim) +
      scale_y_continuous(limits = ylim) +
      labs(x = xtitle, y = ytitle)
    # Plot the mean pdp from the bootstrap
    out <- out +
      geom_line(
        aes(x = data[, 1],
            y = data[, which(colnames(data) %in% yhat)]),
        size = 0.7,
        alpha = 1,
        color = 'black'
      ) +
      geom_ribbon(
        aes(
          x = data[, 1],
          ymin = data[, which(colnames(data) %in% yhat)] +
            qnorm(0.025) * data[, which(colnames(data) %in% ysd)],
          ymax = data[, which(colnames(data) %in% yhat)] -
            qnorm(0.025) * data[, which(colnames(data) %in% ysd)]
        ),
        fill = 'red',
        alpha = 0.1
      )
    # Add markers for the interquartile range of x
    if (iqr == TRUE) {
      range <-
        quantile(trainData[, which(names(trainData) %in% colnames(data)[1])],
                 probs = c(0.25, 0.75))
      out <- out +
        geom_vline(xintercept = range[1],
                   linetype = 'dashed') +
        geom_vline(xintercept = range[2],
                   linetype = 'dashed')
    }
    # Add histogram for the distribution of x
    if (hist == TRUE) {
      out <-
        out + geom_histogram(
          aes(x = trainData[, which(names(trainData) %in% colnames(data)[1])],
              y = stat(count) / sum(count)),
          color = 'black',
          fill = 'black',
          alpha = 0.1
        )
    }
    return(out)
  }
