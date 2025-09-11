# toscca-mm folds
#' Computes TOSCCA-MM
#'
#' This function estimates sparse canonical vectors for matrices with multiple measurements
#' and the trajectories of the latent variables.
#'
#' @param A A matrix.
#' @param B A matrix.
#' @param nonzero_a Integer. Threshold parameter for A.
#' @param nonzero_b Integer. Threshold parameter for B.
#' @param folds Integer. Indicates number of folds to perform.
#' @param parallel_logic Logical. TRUE to parallelise folds. Default is FALSE.
#' @param silent Logical. TRUE to keep silent output messages. Default is FALSE.
#' @param ATest_res NULL. Keep NULL.
#' @param BTest_res NULL. Keep NULL.
#' @param model Character. c("lme", "ar"). Model to fit longitudinal latent space.
#' @param lmeformula Character. LME formula. Default is " ~ -1 + time + (1|id)".
#' @param arformula Numeric vector. Choice of ARIMA. Default is c(1,0,0).
#' @return Canonical vectors for k components.
#' @export
tosccamm = function(A, B, nonzero_a, nonzero_b, folds = 1, parallel_logic = FALSE, silent = FALSE, ATest_res = NULL, BTest_res = NULL, model = "lme", lmeformula = " ~ -1 + time + (1|id)", arformula = c(1,0,0)) {
  # N = min(nrow(A), nrow(B)) # observations
  # p = ncol(A) - 2 # predictor variables (not really since CCA is symmetric)
  # q = ncol(B) - 2# response variables (not really since CCA is symmetric)
  # s = rep(1:folds, length.out=length(unique(A$id)))
  # # s = s[1:N]
  # # s = s[sample(1:length(s), length(s))]
  # if(folds == 1) s[sample(1:N, 0.25*N)] = 2
  # nonzeroGrid = expand.grid(nonzero_a, nonzero_b)
  # h = nrow(nonzeroGrid)
  # canCor_a = matrix(NA, folds, h)
  # canCor_b = matrix(NA, folds, h)
  # canCor = matrix(NA, folds, h)
  #
  #
  # alphaMat <- list()
  # betaMat  <- list()
  #
  # for (f in 1:folds) {
  #   fold_ids <- unique(A$id)[s == f]
  #   ATrain = A[A$id %in% fold_ids, ]
  #   BTrain = B[B$id %in% fold_ids, ]
  #   if(!is.null(ATest_res)) ATest = ATest_res[!(A$id %in% fold_ids), ]
  #   if(is.null(ATest_res)) ATest  = A[!(A$id %in% fold_ids), ]
  #   if(!is.null(BTest_res)) BTest = BTest_res[!(B$id %in% fold_ids), ]
  #   if(is.null(BTest_res)) BTest  = B[!(B$id %in% fold_ids), ]
  #
  #   #
  #   #       # check between selected vector vs. one with higher cancor or follow up
  #   #       if(alpha_init == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
  #   #       if(alpha_init == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
  #   #       if(alpha_init == "uniform") alphaInit = standardVar(replicate(1, runif(p)), normalise = TRUE)
  #
  #
  #   # if(isFALSE(silent)) progressBar(folds, f)
  #
  #   resultKFold = tosccamm.core(alphaInit = runif(ncol(A)-2), A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], model = model, lmeformula = lmeformula, silent = silent)
  #
  #   alphaMat[[f]] <- resultKFold$a
  #   betaMat[[f]]  <- resultKFold$b
  #
  #   # gamma_a = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # gamma_a = data.frame(id = ATest$id, time = ATest$time, x = gamma_a)
  #   # gamma_a = data.frame(time = unique(sort(round(gamma_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_a[,j+2], by = list(round(gamma_a$time,2)), FUN =mean)$x)) #scale_rm(gamma_a, centre = FALSE)
  #   # zeta_a = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # zeta_a = data.frame(id = ATest$id, time = ATest$time, x = zeta_a)
  #   # zeta_a  = data.frame(time = unique(sort(round(zeta_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_a[,j+2], by = list(round(zeta_a$time,2)), FUN =mean)$x)) #scale_rm(zeta_a, centre = TRUE)
  #   #
  #   # gamma_b = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # gamma_b = data.frame(id = BTest$id, time = BTest$time, x = gamma_b)
  #   # gamma_b = data.frame(time = unique(sort(round(gamma_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_b[,j+2], by = list(round(gamma_b$time,2)), FUN =mean)$x)) #scale_rm(gamma_b, centre = FALSE)
  #   # zeta_b = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # zeta_b = data.frame(id = BTest$id, time = BTest$time, x = zeta_b)
  #   # zeta_b  = data.frame(time = unique(sort(round(zeta_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_b[,j+2], by = list(round(zeta_b$time,2)), FUN =mean)$x)) #scale_rm(zeta_b, centre = TRUE)
  #
  #   getlv = function(mat, long_mod, new_data, grid) {
  #     lv = sapply(1:ncol(mat), function(j) predict(long_mod[[j]], newdata = new_data, allow.new.levels = TRUE, re.form = NULL))
  #     lv =  data.frame(id = new_data$id, time = new_data$time, x = lv)
  #     lv =  data.frame(time = unique(sort(round(lv$time,2))), sapply(1:nrow(grid), function(j) aggregate(lv[,j+2], by = list(round(lv$time,2)), FUN =mean)$x))
  #
  #     return(lv)
  #
  #   }
  #
  #   gamma_a = getlv(alphaMat[[f]], resultKFold$me_x, data.frame(time = ATest[,"time"], id = ATest[,"id"]), nonzeroGrid)
  #   zeta_a  = getlv(betaMat[[f]], resultKFold$me_y, data.frame(time = ATest[,"time"], id = ATest[,"id"]), nonzeroGrid)
  #   gamma_b = getlv(alphaMat[[f]], resultKFold$me_x, data.frame(time = BTest[,"time"], id = BTest[,"id"]), nonzeroGrid)
  #   zeta_b  = getlv(betaMat[[f]], resultKFold$me_y, data.frame(time = BTest[,"time"], id = BTest[,"id"]), nonzeroGrid)
  #   # gammaT <- aggregate(. ~ id + time, data = gamma, FUN = mean)
  #   # zettaT <- aggregate(. ~ id + time, data = zeta, FUN = mean)
  #
  #
  #   # cc_time = matrix(NA, ncol = ncol(gamma_a)-1, nrow =length(unique(A$time)))
  #   #   for (m in unique(A$time)) {
  #   #     w = unique(intersect(A[A$time==m,]$id, B[B$time==m,]$id))
  #   #     g = as.matrix(A[w, -c(1,2)])%*%resultKFold$a; e = as.matrix(B[w, -c(1,2)])%*%resultKFold$b;
  #   #     cc_time[m, ] = sapply(1:(ncol(gamma_a) -1), function(x) cor(g[,x], e[,x]))
  #   #   }
  #   cc = 1# colMeans(cc_time)
  #
  #   if(ncol(gamma_a) + ncol(zeta_a) > 2) canCor_a[f,]  = abs(sapply(2:ncol(gamma_a), function(j) cor(gamma_a[,j], zeta_a[,j])*cc[j-1]))
  #   if(ncol(gamma_a) + ncol(zeta_a) == 2) canCor_a[f,] = abs(cor(gamma_a, zeta_a_b)*cc)
  #
  #   if(ncol(gamma_b) + ncol(zeta_b) > 2) canCor_b[f,]  = abs(sapply(2:ncol(gamma_b), function(j) cor(gamma_b[,j], zeta_b[,j])*cc[j-1]))
  #   if(ncol(gamma_b) + ncol(zeta_b) == 2) canCor_b[f,] = abs(cor(gamma_b, zeta_b)*cc)
  #
  #   canCor[f,] = sapply(1:length(canCor_a[f,]), function(x) mean(canCor_a[f,x], canCor_b[f,x]))
  #
  #   if(any(is.na(canCor[f,]))) stop("oneis NA")
  #
  # }
  #
  # canCorKMeans = colSums(abs(canCor))/folds
  # select       = getWhich(abs(canCorKMeans), max)
  # if(length(select)>1) select <- select[1]
  # canCorPrint  = canCorKMeans[select]
  #
  # cat("📊 K-Fold CV Max. Canonical Correlation:  ", canCorPrint, "                      ")
  #
  # if (isFALSE(silent)) cat("\n")
  # # if (isFALSE(silent)) ca(canCorPrint)
  #
  # if (isFALSE(silent)) {
  #   a <- nonzeroGrid[select, 1]
  #   b <- nonzeroGrid[select, 2]
  #
  #   cat("\n",
  #       "   ──────────────────────────\n",
  #       sprintf("      Nonzero A: %8d \n", a),
  #       sprintf("      Nonzero B: %8d \n", b),
  #       "   ──────────────────────────\n",
  #       sep = "")
  # }
  #
  # if(nrow(nonzeroGrid) > 1) {
  #   mat = matrix(canCorKMeans, nrow = length(nonzero_a), byrow = F)
  #   # mat = matrix(canCorKMeans, nrow = length(nonzero_a), ncol = length(nonzero_b))
  #   rownames(mat) = nonzero_a
  #   colnames(mat) = nonzero_b
  #   # myHeatmap(mat)
  #
  # }
  #
  # result     = tosccamm.core(alphaInit =  runif(ncol(A)-2), A = A, B = B, nonzero_a = nonzeroGrid[select, 1], nonzero_b = nonzeroGrid[select, 2], silent = TRUE,  model = model, arformula = arformula, lmeformula = lmeformula)
  #
  # cat("done")
  # if(nrow(nonzeroGrid) > 1) {
  #   resultSCCA = list(cancor = canCorPrint,
  #                     alpha  = result$a,
  #                     beta   = result$b,
  #                     me_x   = result$me_x,
  #                     me_y   = result$me_y,
  #                     cv_logical = nrow(nonzeroGrid) > 1,
  #                     mat_cc = mat,
  #                     # alphaMat       = alphaMat,
  #                     # betaMat        = betaMat,
  #                     # position       = select,
  #                     nonzero_a = nonzeroGrid[select, 1],
  #                     nonzero_b = nonzeroGrid[select, 2],
  #                     canCor_grid = mat)
  # } else {
  #   resultSCCA = list(cancor = canCorPrint,
  #                     alpha  = result$a,
  #                     beta   = result$b,
  #                     me_x   = result$me_x,
  #                     me_y   = result$me_y,
  #                     cv_logical = nrow(nonzeroGrid) > 1,
  #                     # alphaMat       = alphaMat,
  #                     # betaMat        = betaMat,
  #                     # position       = select,
  #                     nonzero_a = nonzeroGrid[select, 1],
  #                     nonzero_b = nonzeroGrid[select, 2])
  # }

  N = min(nrow(A), nrow(B)) # observations
  p = ncol(A) - 2 # predictor variables (not really since CCA is symmetric)
  q = ncol(B) - 2# response variables (not really since CCA is symmetric)
  s = rep(1:folds, length.out=length(unique(A$id)))
  # s = s[1:N]
  # s = s[sample(1:length(s), length(s))]
  if(folds == 1) s[sample(1:N, 0.25*N)] = 2
  nonzeroGrid = expand.grid(nonzero_a, nonzero_b)
  h = nrow(nonzeroGrid)
  canCor_a = matrix(NA, folds, h)
  canCor_b = matrix(NA, folds, h)
  canCor = matrix(NA, folds, h)


  alphaMat <- list()
  betaMat  <- list()

  for (f in 1:folds) {
    fold_ids <- unique(A$id)[s == f]
    ATrain = A[A$id %in% fold_ids, ]
    BTrain = B[B$id %in% fold_ids, ]
    if(!is.null(ATest_res)) ATest = ATest_res[!(A$id %in% fold_ids), ]
    if(is.null(ATest_res)) ATest  = A[!(A$id %in% fold_ids), ]
    if(!is.null(BTest_res)) BTest = BTest_res[!(B$id %in% fold_ids), ]
    if(is.null(BTest_res)) BTest  = B[!(B$id %in% fold_ids), ]

    #
    #       # check between selected vector vs. one with higher cancor or follow up
    #       if(alpha_init == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
    #       if(alpha_init == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
    #       if(alpha_init == "uniform") alphaInit = standardVar(replicate(1, runif(p)), normalise = TRUE)


    # if(isFALSE(silent)) progressBar(folds, f)

    resultKFold = tosccamm.core(alphaInit = runif(ncol(A)-2), A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], model = model, lmeformula = lmeformula, silent = silent)

    alphaMat[[f]] <- resultKFold$a
    betaMat[[f]]  <- resultKFold$b

    gamma_a = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    gamma_a = data.frame(id = ATest$id, time = ATest$time, x = gamma_a)
    gamma_a = data.frame(time = unique(sort(round(gamma_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_a[,j+2], by = list(round(gamma_a$time,2)), FUN =mean)$x)) #scale_rm(gamma_a, centre = FALSE)
    zeta_a = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    zeta_a = data.frame(id = ATest$id, time = ATest$time, x = zeta_a)
    zeta_a  = data.frame(time = unique(sort(round(zeta_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_a[,j+2], by = list(round(zeta_a$time,2)), FUN =mean)$x)) #scale_rm(zeta_a, centre = TRUE)

    gamma_b = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    gamma_b = data.frame(id = BTest$id, time = BTest$time, x = gamma_b)
    gamma_b = data.frame(time = unique(sort(round(gamma_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_b[,j+2], by = list(round(gamma_b$time,2)), FUN =mean)$x)) #scale_rm(gamma_b, centre = FALSE)
    zeta_b = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    zeta_b = data.frame(id = BTest$id, time = BTest$time, x = zeta_b)
    zeta_b  = data.frame(time = unique(sort(round(zeta_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_b[,j+2], by = list(round(zeta_b$time,2)), FUN =mean)$x)) #scale_rm(zeta_b, centre = TRUE)

    # gammaT <- aggregate(. ~ id + time, data = gamma, FUN = mean)
    # zettaT <- aggregate(. ~ id + time, data = zeta, FUN = mean)

    tv = sort(union(unique(ATest$time), unique(BTest$time)))
    l = sapply(1:length(tv), function(v) length(intersect(ATest[ATest$time==v, ]$id, BTest[BTest$time==v, ]$id)))
    v = which(l==max(l))[1]
    iv = intersect(ATest[ATest$time==v, ]$id, BTest[BTest$time==v, ]$id)
    g = data.frame(id = ATest$id, time = ATest$time, as.matrix(ATest[,-c(1,2)])%*%alphaMat[[f]]); e = data.frame(id = BTest$id, time = BTest$time, as.matrix(BTest[,-c(1,2)])%*%betaMat[[f]])
    m = sapply(1:ncol(alphaMat[[f]]), function(x) abs(cor(g[g$id %in% iv & g$time==v, x+2], e[e$id %in% iv & e$time==v, x+2])))

    if(ncol(gamma_a) + ncol(zeta_a) > 2) canCor_a[f,]  = abs(sapply(2:ncol(gamma_a), function(j) cor(gamma_a[,j], zeta_a[,j])))
    if(ncol(gamma_a) + ncol(zeta_a) == 2) canCor_a[f,] = abs(cor(gamma_a, zeta_a))

    if(ncol(gamma_b) + ncol(zeta_b) > 2) canCor_b[f,]  = abs(sapply(2:ncol(gamma_b), function(j) cor(gamma_b[,j], zeta_b[,j])))
    if(ncol(gamma_b) + ncol(zeta_b) == 2) canCor_b[f,] = abs(cor(gamma_b, zeta_b))

    canCor[f,] = sapply(1:length(canCor_a[f,]), function(x) mean(canCor_a[f,x], canCor_b[f,x])*m[x]) #

    if(any(is.na(canCor[f,]))) stop("oneis NA")

  }

  canCorKMeans = colSums(abs(canCor))/folds
  select       = getWhich(abs(canCorKMeans), max)
  canCorPrint  = canCorKMeans[select]

  names(canCorPrint) <- ("k-fold cv max. cancor")
  if(isFALSE(silent)) cat("\n")
  if(isFALSE(silent)) print(canCorPrint)
  if(isFALSE(silent)) {
    cat("\n ........................................ \n",
        paste0("# nonzero A: ", nonzeroGrid[select, 1],   "\n",
               " # nonzero B: ", nonzeroGrid[select, 2],
               "\n ........................................ \n"))
  }

  if(nrow(nonzeroGrid) > 1) {
    mat = matrix(canCorKMeans, nrow = length(nonzero_a), byrow = F)
    # mat = matrix(canCorKMeans, nrow = length(nonzero_a), ncol = length(nonzero_b))
    rownames(mat) = nonzero_a
    colnames(mat) = nonzero_b
    myHeatmap(mat)

  }

  result     = tosccamm.core(alphaInit =  runif(ncol(A)-2), A = A, B = B, nonzero_a = nonzeroGrid[select, 1], nonzero_b = nonzeroGrid[select, 2], silent = TRUE,  model = model, arformula = arformula, lmeformula = lmeformula)


  if(nrow(nonzeroGrid) > 1) {
    resultSCCA = list(cancor = canCorPrint,
                      alpha  = result$a,
                      beta   = result$b,
                      me_x   = result$me_x,
                      me_y   = result$me_y,
                      # alphaMat       = alphaMat,
                      # betaMat        = betaMat,
                      # position       = select,
                      nonzero_a = nonzeroGrid[select, 1],
                      nonzero_b = nonzeroGrid[select, 2],
                      canCor_grid = mat)
  } else {
    resultSCCA = list(cancor = canCorPrint,
                      alpha  = result$a,
                      beta   = result$b,
                      me_x   = result$me_x,
                      me_y   = result$me_y,
                      # alphaMat       = alphaMat,
                      # betaMat        = betaMat,
                      # position       = select,
                      nonzero_a = nonzeroGrid[select, 1],
                      nonzero_b = nonzeroGrid[select, 2])
  }

  return(resultSCCA)
}
