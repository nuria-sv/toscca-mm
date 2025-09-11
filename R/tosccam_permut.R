
# permutation testing -----------------------------------------------------
#' Computes permutatied cc fot TOSCCA-MM
#'
#' This function estimates sparse canonical vectors for permuted matrices with multiple measurements.
#' @param A A matrix.
#' @param B A matrix.
#' @param nonzero_a Integer. Threshold parameter for A.
#' @param nonzero_b Integer. Threshold parameter for B.
#' @param K Interger. Numner of components.
#' @param alpha_init Character. Type of canonical vector initialisation c("uniform", "random", "eigen")
#' @param toPlot Logical. Indicates if plots will be produced. Default is False.
#' @param draws Integer. Number of draws in permutation.
#' @param cancor Numeric vector of length K with estimated canonical correlations.
#' @param folds Integer. Indicates number of folds to perform.
#' @param parallel_logic Logical. TRUE to parallelise folds. Default is FALSE.
#' @param nuisanceVar Numeric. Number of nuisance variables.
#' @param bootCCA deprecated.
#' @param testStatType Character. Choice of test-statistic c("CC", "Wilks", "Roy"),
#' @param silent Logical. TRUE to keep silent output messages. Default is FALSE.
#' @param model Character. c("lme", "ar"). Model to fit longitudinal latent space.
#' @param lmeformula Character. LME formula. Default is " ~ -1 + time + (1|id)".
#' @param arformula Numeric vector. Choice of ARIMA. Default is c(1,0,0).
#' @return Permutated canonical correlation for ell K and p-values.
#'
#' @importFrom graphics abline hist legend lines par text
#' @importFrom utils install.packages installed.packages
#' @importFrom stats aggregate arima as.formula coefficients cor density predict runif
#' @importFrom foreach %dopar%
#' @export

toscamm.perm = function (A, B, nonzero_a, nonzero_b, K, alpha_init = c("eigen",
                                                                       "random", "uniform"), folds = 1, toPlot = FALSE, draws = 1000,
                         cancor, bootCCA = NULL, silent = TRUE, parallel_logic = TRUE,
                         nuisanceVar = 0, testStatType = "CC", model = "lme", lmeformula = " ~ 0 + poly(time,3) + (1|id)", arformula = NULL )
{
  histNullCCA <- modelCanCor <- x <- y <- `..count..` <- NULL
  # if (!requireNamespace("EnvStats", quietly = TRUE))
  #   install.packages("EnvStats")
  if (!requireNamespace("parallel", quietly = TRUE))
    install.packages("parallel")
  if (!requireNamespace("doParallel", quietly = TRUE))
      install.packages("doParallel")
  if (!requireNamespace("ggplot2", quietly = TRUE))
        install.packages("ggplot2")
  if (!requireNamespace("lme4", quietly = TRUE))
    install.packages("lme4")
  if (!requireNamespace("toscca", quietly = TRUE))
    devtools::install_github("nuria-sv/toscca")
  # library(foreach);library(doParallel)
  perm = matrix(NA, nrow = draws, ncol = K)
  if (isTRUE(parallel_logic)) {
    myCluster <- parallel::makeCluster(parallel::detectCores() -
                                         1, type = "PSOCK")
    doParallel::registerDoParallel(cl = myCluster)
    if (isFALSE(foreach::getDoParRegistered()))
      stop("DoPar not registered. Check cores")
    perm <- foreach::foreach(d = 1:draws, .combine = "rbind", .export = ls(globalenv())) %dopar%
      {
        res_perm = list()
        ASample = data.frame(A[,1:2], A[sample(1:nrow(A), nrow(A)), -c(1,2)])
        X.temp = ASample
        Y.temp = B
        cc_time = matrix(NA, length(unique(Y.temp$time)), K)
        for (k in 1:K) {
          if(k > 1) {
            # residualise for subsequent components
            X.temp = data.frame(X.temp[,c(1,2)],toscca::residualisation(as.matrix(X.temp[,-c(1,2)]), res_perm[[k-1]]$alpha, type = "basic") )
            Y.temp = data.frame(Y.temp[,c(1,2)],toscca::residualisation(as.matrix(Y.temp[,-c(1,2)]), res_perm[[k-1]]$beta, type = "basic") )

            nz_a_gen = as.numeric(table(res_perm[[k-1]]$alpha != 0)[2])
            nz_b_gen = as.numeric(table(res_perm[[k-1]]$beta != 0)[2])
          }

          res_perm[[k]] <- tosccamm(X.temp, Y.temp, folds = 2,
                                 nonzero_a[k], nonzero_b[k],
                                 model = model, lmeformula = lmeformula)

         #  for (s in unique(X.temp$time)) {
         #    w = unique(intersect(X.temp[X.temp$time==s,]$id, Y.temp[Y.temp$time==s,]$id))
         #    g = as.matrix(X.temp[w, -c(1,2)])%*%res_perm[[k]]$alpha; e = as.matrix(Y.temp[w, -c(1,2)])%*%res_perm[[k]]$beta;
         #    cc_time[s, k] = abs(cor(g, e))
         #  }
         #
         # m= colMeans(cc_time)

        }
       sapply(1:K, function(k) res_perm[[k]]$cancor)

      }
  }
  else {
    for (d in 1:draws) {
      # if (isFALSE(silent))
        # progressBar(draws, d)

      res_perm = list()
      ASample = data.frame(A[,1:2], A[sample(1:nrow(A), nrow(A)), -c(1,2)])
      X.temp = ASample
      Y.temp = B
      cc_time = matrix(NA, length(unique(Y.temp$time)), K)
      for (k in 1:K) {
        if(k > 1) {
          # residualise for subsequent components
          X.temp = data.frame(X.temp[,c(1,2)],toscca::residualisation(as.matrix(X.temp[,-c(1,2)]), res_perm[[k-1]]$alpha, type = "basic") )
          Y.temp = data.frame(Y.temp[,c(1,2)],toscca::residualisation(as.matrix(Y.temp[,-c(1,2)]), res_perm[[k-1]]$beta, type = "basic") )

          nz_a_gen = as.numeric(table(res_perm[[k-1]]$alpha != 0)[2])
          nz_b_gen = as.numeric(table(res_perm[[k-1]]$beta != 0)[2])
        }

        res_perm[[k]] <- tosccamm(X.temp, Y.temp, folds = 2,
                                  nonzero_a[k], nonzero_b[k],
                                  model = model, lmeformula = lmeformula, silent = silent)

        #  for (s in unique(X.temp$time)) {
        #    w = unique(intersect(X.temp[X.temp$time==s,]$id, Y.temp[Y.temp$time==s,]$id))
        #    g = as.matrix(X.temp[w, -c(1,2)])%*%res_perm[[k]]$alpha; e = as.matrix(Y.temp[w, -c(1,2)])%*%res_perm[[k]]$beta;
        #    cc_time[s, k] = abs(cor(g, e))
        #  }
        #
        # m= colMeans(cc_time)

      }
      perm[d, ] = sapply(1:K, function(k) res_perm[[k]]$cancor)
    }
  }
  if (isTRUE(parallel_logic))
    parallel::stopCluster(cl = myCluster)
  margin = 0.05
  # cc_time = matrix(NA, length(unique(XX2$time)), K)
  # for (k in 1:K) {
  #   for (s in unique(A$time)) {
  #   w = unique(intersect(XX2[XX2$time==s,]$id, YY2[YY2$time==s,]$id))
  #   g = as.matrix(XX2[w, -c(1,2)])%*%res_perm[[k]]$alpha; e = as.matrix(YY2[w, -c(1,2)])%*%res_perm[[k]]$beta;
  #   cc_time[s, k] = abs(cor(g, e))
  #   }}
  m = (cancor) #sapply(1:K, function(k)mean(cc_time[,k])) #
  testStatistic = m #sapply(1:K, function(k) res_perm[[k]]$cancor) # CCAtStat(cancor, A, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]
  h = hist(perm[, getWhich(testStatistic, max)], breaks = draws/2,
           plot = FALSE)
  if (is.null(bootCCA)) {
    permDensity = density(perm[, getWhich(testStatistic,
                                          max)])
    xlim = c(min(min(perm) - margin, testStatistic - margin),
             max(max(perm) + margin, testStatistic + margin))
    ylim = c(min((modes(permDensity)$y), min(perm[, getWhich(testStatistic,
                                                             max)])), max((modes(permDensity)$y) + margin, max(h$counts) +
                                                                            margin))
    par(mfrow = c(1, 1))
  #   hist(perm[, getWhich(testStatistic, max)], breaks = draws/2,
  #        xlab = "Canonical Correlation", xlim = xlim, ylim = ylim,
  #        main = paste0("Distribution under de Null - ", testStatType,
  #                      " Statistic"), col = scales::alpha("#01768c",
  #                                                         0.2))
  #   lines(permDensity, col = "black", lwd = 2)
  #   EnvStats::epdfPlot(perm[, getWhich(testStatistic, max)],
  #                      discrete = FALSE, density.arg.list = NULL, plot.it = TRUE,
  #                      add = TRUE, epdf.col = "steelblue", epdf.lwd = 3 *
  #                        par("cex"), epdf.lty = 1, curve.fill = FALSE,
  #                      curve.fill.col = "steelblue", main = NULL, xlab = NULL,
  #                      ylab = NULL)
  #   abline(v = testStatistic, col = "#e86d07", lwd = 2)
  #   legend("topleft", c("Empirical pdf", "density", "model canCor"),
  #          col = c("steelblue", "black", "red"), lty = c(2,
  #                                                        1, 1), cex = 0.8)
  #   text(x = as.character(testStatistic), y = 0.9 * par("usr")[4],
  #        labels = as.character(1:K), cex = 0.9)

    col = mpalette[3:7] # viridis::magma(5)
    # library(ggplot2)
    # library(viridis)
    # library(data.table)
    # library(scales)

    # Assuming perm is a data.table or matrix, and testStatistic, margin, h, draws, testStatType, K are defined

    # Extract the relevant values
    perm_values <- perm[, getWhich(testStatistic, max)]
    testStat_val <- testStatistic

    # Compute density
    perm_density <- density(perm_values)

    # Create a data frame for plotting
    df <- data.frame(
      x = perm_values
    )

    density_df <- data.frame(
      x = perm_density$x,
      y = perm_density$y
    )

    # Compute histogram manually for better control in ggplot2
    hist_data <- hist(perm_values, breaks = draws / 2, plot = FALSE)
    hist_df <- data.frame(
      x = hist_data$mids,
      y = hist_data$counts
    )

    # Define x and y limits
    xlim <- c(min(min(perm) - margin, testStat_val - margin),
              max(max(perm) + margin, testStat_val + margin))

    ylim <- c(min(c(min(perm_values), min(perm_density$y))),
              max(c(max(hist_data$counts), max(perm_density$y))) + margin)

    # Plot

   plt= ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..count..),
                     data = df,
                     bins = draws / 2,
                     fill = mpalette[5], # viridis::viridis(1, option = "magma", alpha = 0.2),
                     color = NA) +
      # geom_line(data = density_df, aes(x = x, y = y), color = col[5], size = 1) +
     ggplot2::geom_line(data = density_df, ggplot2::aes(x = x, y = y), color = col[1], size = 1.2, linetype = "dashed") +
     ggplot2::geom_vline(xintercept = testStat_val, color = col[5], size = 1) +
     ggplot2::labs(x = "Canonical Correlation",
           y = "Density / Count",
           title = paste0("")) +
     ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
     ggplot2::theme_minimal(base_size = 14) +
     # ggplot2::scale_fill_viridis_d(option = "magma") +
     ggplot2::scale_fill_manual(values = mpalette) +
     ggplot2::annotate("text", x = testStat_val, y = 0.9 * ylim[2], label = as.character(1:K), size = 4) +
     ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
     ggplot2::theme(legend.position = "none")
print(plt)

  }
  else {
    li = bootCCA$ccaInterval[1]
    ui = bootCCA$ccaInterval[2]
    xlim = c(min(min(perm) - margin, modelCanCor - margin,
                 li - margin), max(max(perm) + margin, modelCanCor +
                                     margin, ui + margin))
    li = bootCCA$ccaInterval[1]
    ui = bootCCA$ccaInterval[2]
    xlim = c(min(min(perm) - margin, modelCanCor - margin),
             max(max(perm) + margin, modelCanCor + margin))
    permDensity = density(perm[, getWhich(testStatistic,
                                          max)])
    par(mfrow = c(1, 1))
    hist(perm[, getWhich(testStatistic, max)], breaks = draws/2,
         xlab = "Canonical Correlation", xlim = xlim, ylim = sort(modes(permDensity)$y),
         main = paste0("Distribution under de Null - ", testStatType,
                       " Statistic"), col = "#93D9D9")
    lines(permDensity, col = "black", lwd = 2)
    lines(x = c(li, li), y = c(max(histNullCCA$density) *
                                 0.5/2, max(histNullCCA$density) * 0.6/2))
    lines(x = c(ui, ui), y = c(max(histNullCCA$density) *
                                 0.5/2, max(histNullCCA$density) * 0.6/2))
    lines(x = c(li, ui), y = rep((max(histNullCCA$density) *
                                    0.5/2 + max(histNullCCA$density) * 0.6/2)/2, 2))
    abline(v = modelCanCor, col = "purple", lwd = 2)
    abline(v = mean(bootCCA$boostCanCor), col = "red", lwd = 2)
    legend("topleft", c("normal", "density", "model canCor",
                        "bootstrap avg canCor"), col = c("steelblue", "black",
                                                         "purple", "red"), lty = c(2, 2, 1, 1), cex = 0.8)
    text(x = as.character(testStatistic), y = 0.9 * par("usr")[4],
         labels = as.character(1:K), cex = 0.9)
  }
  pValues = sapply(1:K, function(k) round(mean(testStatistic[k] <=
                                                 (perm[, getWhich(testStatistic, max)])), 6))
  print(cat("Empirical p-values:", pValues, sep = "\n"))
  return(list(perm_estimate = perm, p_values = pValues))
  }
