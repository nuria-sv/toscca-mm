# scale and stand measurements--------------------------------------------------
#' Standardises matrices with multiple measurements per individual.
#'
#' This function stardadises matrices with multiple measurements
#' w.r.t. a chosen origin.
#'
#' @param mat A matrix.
#' @param origin Measurement of reference for stardadisation.
#' @param centre Logical. TRUE to centre data. Default is FALSE.
#' @export
scale_rm <- function(mat, origin = NULL, centre = F) {
  ncols = ncol(mat)
  nrows = nrow(mat)
  if(is.null(origin)) { # get min collection date
    origin_rows = as.numeric(unlist(sapply(unique(mat$id), function(x) which(mat$time==min(mat[mat$id==x,]$time)))))
  } else {
    origin_rows = which(mat$time == origin)
  }

  # scale first columnwise
  cent = sapply(3:ncols, function(j) mean(as.matrix(mat[origin_rows,j])))
  dist   = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
  mat[,3:ncols] = sapply(3:ncols, function(j) mat[,j]/dist[j-2])                       # faster than sweep
  # centre whole dataset
  if(centre) {
    gen_cent = mean(as.matrix(mat[origin_rows,3:ncols]))
    mat[,3:ncols]  = sapply(3:ncols, function(j) mat[,j] - gen_cent)
  }
  # cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
  # if(centre) {
  #   mat[,3:ncols] <- sapply(3:ncols, function(j) mat[,j] - cent[j-2])
  #   # mat[,3:ncols] <- sweep(mat[3:ncols], 2L, cent, check.margin=FALSE)
  #   dist = sapply(3:ncols, function(j) sqrt(sum(mat[,j]**2)/max(1, length(mat[,j]-1))))
  #
  # } else {
  #   dist = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
  #
  # }
  # mat_st = sweep(mat[,3:ncols], 2, dist, "/", check.margin = FALSE)
  # mat_st = as.matrix(cbind(mat[,1:2], mat_st), nrow = nrows, ncol = ncols)
  # colnames(mat_st) <- colnames(mat)
  return(data.frame(mat))
}


# modes -------------------------------------------------------------------

#' Calculates mode.
#' @param d density object.
modes = function(d) {
  i <- which(diff(sign(diff(d$y))) < 0 ) + 1
  data.frame(x = d$x[i], y = d$y[i])
}


# cpev calc --------------------------------------------------------------------
#' Calculates cpev
#'
#' This function calculates cpev
#' w.r.t. a chosen origin.
#'
#' @param mat A matrix.
#' @param weights A numberic vector of canonical weights.
cpev.toscca = function (mat, weights) {
  if(is.null(dim(weights))) weights = as.matrix(weights, ncol =1)
  cpev.fun = function(mat, matK) {
    # mat = scale(mat); matK = scale(matK)
    var_mat = t(mat) %*% mat
    traceX = sum(rowSums(var_mat))
    var_matK = t(matK) %*% matK
    traceK = sum(rowSums(var_matK))
    cpev = traceK/traceX

    return(cpev)
  }

  # if(var(weights[weights!=0])!=1) {
  #   weights[weights!=0] = scale(weights[weights!=0])
  # }
  if(!is.null(mat$time)) {
    nmax = sum(table(mat$time))
    dist  = sqrt(colSums(weights^2))
    weights = sweep(weights, 2, dist, "/")
    matk = as.matrix(mat[, -c(1, 2)])%*%(weights)

    cpev_ls = sapply(sort(unique(mat$time)),
                     function(s) cpev.fun(as.matrix(mat[mat$time == s , -c(1, 2)]), matk))/nmax
    cpev = (mean(cpev_ls))

  } else {

    cpev = cpev.fun(mat, mat%*%weights)

  }

  return(cpev)
}

# create cv heatmap ------------------------------------------------------------
#' Plot heatmap of cv w.r.t. the penalty parameter perfotmance.
#'
#' This function plots cca for different thresholds
#'
#' @param mat A matrix.
#' @param palette Character. Name of a palette for the heatmap. Default is "Teal".
#' @param xlab Character. Label for x axis.
#' @param ylab Character. Label for y axis.
#' @param show_axes Logic. Default is True.
#' @param show_labels Logic. Default is True.
#' @param K Numeric. Number of components.
#' @export
myHeatmap <- function(mat, palette = "rocket", xlab = "", ylab = "",
                      show_axes = TRUE, show_labels = TRUE, K = 1) {
  if (!is.matrix(mat)) {
    stop("input 'mat' must be a matrix.")
  }

  # colors <- viridis::viridis(100, option = palette, direction = -1)
  #
  # layout(matrix(c(1,2), ncol = 2), widths = c(4, 1))
  #
  # par(mar = c(5, 4, 4, 2))
  # image(1:nrow(mat), 1:ncol(mat), (mat[nrow(mat):1,ncol(mat):1 ]),
  #       col = colors, axes = FALSE, xlab = xlab, ylab = ylab)
  #
  # if (show_axes) {
  #   axis(2, at = 1:ncol(mat), labels = rev(colnames(mat)), las = 1)
  #   axis(1, at = 1:nrow(mat), labels = rev(rownames(mat)), las = 2)
  # }
  #
  # par(mar = c(5, 2, 4, 4))
  # image(
  #   matrix(seq(0,1, length.out = length(mat)), nrow = 1),
  #   col = colors, axes = FALSE)
  # axis(4, at = seq(0, 1, length.out = length(mat)) , labels = round(seq(min(mat), max(mat), length.out = length(mat)), 2), las = 2, tick = T, pos =1)

  mat_df <- as.data.frame(as.table(mat[nrow(mat):1, ncol(mat):1]))
  colnames(mat_df) <- c("X", "Y", "value")

  # Create ggplot heatmap
  plt <- ggplot2::ggplot(mat_df, ggplot2::aes(X, Y, fill = value)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(option = palette, direction = -1, limits = c(min(mat_df$value)-0.01, max(mat_df$value)+0.01)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = if (show_labels) ggplot2::element_text(angle = 0, hjust = 1) else ggplot2::element_blank(),
          axis.text.y = if (show_labels) ggplot2::element_text(hjust = 1) else ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          legend.position = "right") #+
    # ggtitle("cc for different sparsity levels") +
    # theme(
    #     plot.title = element_text(size = 16, face = "bold", color = "black",
    #                               hjust = 0.5, vjust = 1.5, lineheight = 1.2)
    # )

  return(plt)

}
# create sel stability plt -----------------------------------------------------
#' Plot slection stability for penalty parameter perfotmance.
#'
#' This function plots cv for different thresholds
#'
#' @param mat A matrix of canonical correlation for different sparsity levels.
#' @param X nxp matrix. Observation matrix.
#' @param Y A nxq matrix. Observation matrix.
#' @param palette Character. Name of a palette for the heatmap. Default is "Teal".
#' @param mm Logic. Indicates whether there are multiple measurements or not. Default is True.
#' @export
plt.selstab=
  function (mat, X, Y, palette = "magma", mm = T)
  {
    nonz_x = sort(as.numeric(rownames(mat)))
    nonz_y = sort(as.numeric(colnames(mat)))

    if(mm){
      run_cca <- function(i) {
        cca_res = tosccamm::tosccamm(
          X, Y,
          nonzero_a = nonz_x[i],
          nonzero_b = nonz_y[1],
          silent = TRUE
        )
        list(alpha = cca_res$alpha[, 1], beta = cca_res$beta[, 1], cancor = cca_res$cancor[1])
      }

    } else {
      run_cca <- function(i) {
        cca_res = toscca::toscca(
          X, Y,
          nonz_x[i],
          nonz_y[1], K = 1, "uniform",
          combination = FALSE, silent = T, toPlot = F)
        list(alpha = cca_res$alpha[, 1], beta = cca_res$beta[, 1], cancor = cca_res$cancor[1])
      }

    }
    results <- sapply(seq_along(nonz_x), run_cca, simplify = FALSE)
    # results <- mclapply(seq_along(nonz_x), run_cca, mc.cores = detectCores() - 1)

    alpha <- do.call(cbind, lapply(results, `[[`, "alpha"))
    beta <- do.call(cbind, lapply(results, `[[`, "beta"))
    cancor <- unlist(lapply(results, `[[`, "cancor"))

    # Convert matrix to binary
    numToBi <- function(data) {
      1L*(data != 0)
    }

    count_alpha_df = as.data.frame(numToBi(alpha))
    colnames(count_alpha_df) <- as.character(nonz_x)
    count_alpha_df$count <- rowSums(count_alpha_df)


    palette_colors <- viridis::viridis(8, option = palette)
    col_a <- palette_colors[pmin(count_alpha_df$count, 8)]
    col_a[is.na(col_a)] <- "black"

    alpha_points = alpha[, ncol(alpha)]
    alpha_points[alpha_points!=0] = scale(alpha_points[alpha_points!=0], center = T)

    plot_data <- data.frame(
      index = 1:nrow(alpha),
      coefficients =alpha_points,
      count = count_alpha_df$count
    )

    # plt = ggplot(plot_data, aes(x = index, y = coefficients, color = factor(count))) +
    #   geom_point(size = 3, alpha = 0.99) +  # Larger points with transparency
    #   scale_color_viridis_d(option = palette, name = "selection count", direction = -1) +
    #   labs(
    #     title = "",
    #     x = "feature index (p)",
    #     y = "coefficient"
    #   ) +
    #   theme_minimal(base_size = 14) +
    #   theme(
    #     legend.position = "right",
    #     panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    #     panel.grid.minor = element_blank()
    #   ) + ggtitle("cc w.r.t. sparsity levels") + theme(plot.title = element_text(size = 16, face = "bold", color = "black",
    #                                                                           hjust = 0.5, vjust = 1.5,
    #                                                                           lineheight = 1.2))

    n_colors <- length(unique(plot_data$count))
    custom_colors <- viridis::viridis(n_colors + 1, option = palette)[-1]

    plt = ggplot2::ggplot(plot_data, ggplot2::aes(x = index, y = coefficients, color = factor(count))) +
      ggplot2::geom_point(aes(alpha = factor(count)), size = 3) +
      # scale_color_viridis_d(option = palette, name = "Selection Count", direction = -1) +
      ggplot2::scale_color_manual(values = custom_colors, name = "Selection Count") +
      ggplot2::scale_alpha_manual(values = scales::rescale(as.numeric(levels(factor(plot_data$count))), to = c(0.1, 1)), guide = "none") +
      ggplot2::labs(
        title = "",
        x = "Feature index (p)",
        y = "Coefficient"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position = "right",
        panel.grid.major = ggplot2::element_line(color = "gray80", linetype = "dashed"),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::ggtitle("Selection stability for canonical vectors") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold", color = "black",
                                  hjust = 0.5, vjust = 1.5, lineheight = 1.2)
      )

    return(plt)
  }


# get which function ------------------------------------------------------
#' Get location of required.
#'
#' @param data Numeric matrix.
#' @param fun Function to search data.
#'
#' @export
getWhich =
  function (data, fun)
  {
    fun = match.fun(fun)
    position = (which(data == fun(data)))
    position
  }

#
# # longitudinal standardisation
# scale_rm <- function(mat, origin = 1, centre = F) {
#   ncols = ncol(mat)
#   nrows = nrow(mat)
#
#   origin_rows = which(mat$time == origin)
#
#
#   # scale first columnwise
#   cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
#   dist   = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
#   mat[,3:ncols] = sapply(3:ncols, function(j) mat[,j]/dist[j-2])                       # faster than sweep
#
#   # centre whole dataset
#   if(centre) {
#     gen_cent = mean(as.matrix(mat[origin_rows,3:ncols]))
#     mat[,3:ncols]  = sapply(3:ncols, function(j) mat[,j] - gen_cent)
#   }
#
#
#   # cent = sapply(3:ncols, function(j) mean(mat[origin_rows,j]))
#   # if(centre) {
#   #   mat[,3:ncols] <- sapply(3:ncols, function(j) mat[,j] - cent[j-2])
#   #   # mat[,3:ncols] <- sweep(mat[3:ncols], 2L, cent, check.margin=FALSE)
#   #   dist = sapply(3:ncols, function(j) sqrt(sum(mat[,j]**2)/max(1, length(mat[,j]-1))))
#   #
#   # } else {
#   #   dist = sapply(3:ncols, function(j) sqrt(sum((mat[,j] - cent[j-2])**2)/max(1, length(mat[,j]-1))))
#   #
#   # }
#
#   # mat_st = sweep(mat[,3:ncols], 2, dist, "/", check.margin = FALSE)
#
#   # mat_st = as.matrix(cbind(mat[,1:2], mat_st), nrow = nrows, ncol = ncols)
#   # colnames(mat_st) <- colnames(mat)
#
#   return(data.frame(mat))
# }

# toscca.perm = function(A, B, nonzero_a, nonzero_b, K, alpha_init = c("eigen", "random", "uniform"), folds = 1, toPlot = FALSE, draws = 20, cancor, bootCCA = NULL, silent = TRUE, parallel_logic = TRUE, nuisanceVar = 0, testStatType = "CC") {
#   require(EnvStats)
#
#   perm = matrix(NA, nrow = draws, ncol = K)
#
#   if(isTRUE(parallel_logic)) {
#     # create env
#     myCluster <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
#     doParallel::registerDoParallel(cl = myCluster)
#     if(isFALSE(foreach::getDoParRegistered())) stop("DoPar not registered. Check cores")
#
#     perm <- foreach(d = 1:draws, .combine = "rbind", .export = ls(globalenv())) %dopar% {
#       if(isFALSE(silent)) progressBar(draws, d)
#
#       ASample = A[sample(1:nrow(A), nrow(A)),]
#       t(toscca.tStat(toscca(A = ASample, B = B, K = K, alpha_init = alpha_init, combination = FALSE, nonzero_a=nonzero_a, nonzero_b=nonzero_b, toPlot = FALSE, silent = TRUE, parallel_logic = FALSE)$cancor, ASample, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]) #off-sample cancor
#
#     }
#
#   } else {
#     for (d in 1:draws) {
#       # cat("|", rep(".", d), rep(" ", (draws-d)), "|", (d/draws)*100, "%\r")
#       if(isFALSE(silent)) progressBar(draws, d)
#       ASample = A[sample(1:nrow(A), nrow(A)),]
#       perm[d,] = toscca.tStat(toscca(A = ASample, B = B, K = K, alpha_init = alpha_init, combination = FALSE, nonzero_a=nonzero_a, nonzero_b=nonzero_b, toPlot = FALSE, silent = TRUE)$cancor, ASample, B, C = nuisanceVar, type = testStatType)[["tStatistic"]] #off-sample cancor
#
#
#     }
#   }
#
#
#   if(isTRUE(parallel_logic)) parallel::stopCluster(cl = myCluster)
#
#
#   margin = 0.05
#   testStatistic = toscca.tStat(cancor, A, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]
#   h =  hist(perm[, getWhich(testStatistic, max)], breaks = draws/2, plot = FALSE)
#
#   if(is.null(bootCCA)) {
#     permDensity = density(perm[, getWhich(testStatistic, max)])
#
#     xlim = c(min(min(perm) - margin, testStatistic - margin), max(max(perm) + margin, testStatistic + margin))
#     ylim = c(min((modes(permDensity)$y), min(perm[,getWhich(testStatistic, max)])), max((modes(permDensity)$y) + margin, max(h$counts) + margin))
#
#     par(mfrow=c(1,1))
#     hist(perm[, getWhich(testStatistic, max)], breaks = draws/2, xlab = "Canonical Correlation",
#          xlim = xlim, ylim =  ylim, main = paste0("Distribution under de Null - ", testStatType, " Statistic") , col =  scales::alpha("#01768c", 0.2))
#     lines(permDensity, col="black", lwd = 2)
#     EnvStats::epdfPlot(perm[,getWhich(testStatistic, max)], discrete = FALSE, density.arg.list = NULL, plot.it = TRUE,
#                        add = TRUE, epdf.col = "steelblue", epdf.lwd = 3 * par("cex"), epdf.lty = 1,
#                        curve.fill = FALSE, curve.fill.col = "steelblue", main = NULL, xlab = NULL, ylab = NULL)
#     abline(v=testStatistic, col = "#e86d07", lwd = 2)
#     legend("topleft", c("Empirical pdf", "density", "model canCor"), col = c("steelblue", "black", "red"), lty=c(2, 1, 1), cex=0.8)
#     text(x = as.character(testStatistic), y = 0.9*par('usr')[4], labels = as.character(1:K), cex = 0.9)
#
#
#   } else {
#     li = bootCCA$ccaInterval[1]
#     ui = bootCCA$ccaInterval[2]
#
#     xlim = c(min(min(perm) - margin, modelCanCor - margin, li - margin), max(max(perm) + margin, modelCanCor + margin, ui + margin))
#     li = bootCCA$ccaInterval[1]
#     ui = bootCCA$ccaInterval[2]
#
#     xlim = c(min(min(perm) - margin, modelCanCor - margin), max(max(perm) + margin, modelCanCor + margin))
#     permDensity = density(perm[, getWhich(testStatistic, max)])
#
#     par(mfrow=c(1,1))
#     hist(perm[, getWhich(testStatistic, max)], breaks = draws/2, xlab = "Canonical Correlation", xlim = xlim, ylim =  sort(modes(permDensity)$y), main = paste0("Distribution under de Null - ", testStatType, " Statistic") , col = "#93D9D9")
#     lines(permDensity, col="black", lwd = 2)
#     # lines(seq(xlim[1], xlim[2], 0.01), dnorm(seq(xlim[1], xlim[2], 0.01), mean(perm[, getWhich(modelCanCor, max)]), sd(perm[, getWhich(modelCanCor, max)])), col="steelblue", lwd = 2, lty = 4)
#     lines(x = c(li, li), y = c(max(histNullCCA$density)*0.5/2, max(histNullCCA$density)*0.6/2))
#     lines(x = c(ui, ui), y = c(max(histNullCCA$density)*0.5/2, max(histNullCCA$density)*0.6/2))
#     lines(x = c(li, ui), y = rep((max(histNullCCA$density)*0.5/2 + max(histNullCCA$density)*0.6/2)/2, 2))
#     abline(v = modelCanCor, col = "purple", lwd = 2)
#     abline(v=mean(bootCCA$boostCanCor), col = "red", lwd = 2)
#     legend("topleft", c("normal", "density", "model canCor", "bootstrap avg canCor"), col = c("steelblue", "black", "purple", "red"), lty=c(2, 2, 1, 1), cex=0.8)
#     text(x = as.character(testStatistic), y = 0.9*par('usr')[4], labels = as.character(1:K), cex = 0.9)
#   }
#
#   pValues = sapply(1:K, function (k) round(mean(testStatistic[k]<=(perm[, getWhich(testStatistic, max)])), 6))
#
#   print(cat("Empirical p-values:", pValues, sep = "\n"))
#
#   return(list(perm_estimate = perm, p_values = pValues))
#
# }


# get latent  -------------------------------------------------------------
#' Get latent variables
#'
#' Gets latent variables from data and estimates.
#'
#' @param data List containint both observation amtrices.
#' @param alpha px1 numeric vector. canonical weights for X.
#' @param beta qx1 numeric vector. canonical weights for Y.

toscca.lv <- function(data, alpha, beta) {
  X = data[[1]]
  Y = data[[2]]
  gamma  = matrix(X%*%alpha, nrow = nrow(X))
  zeta   = matrix(Y%*%beta, nrow = nrow(X))
  cancor = cor(gamma, zeta)

  return(cancor)
}


