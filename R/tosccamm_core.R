# toscca-mm core
tosccamm.core = function(alphaInit, A, B, nonzero_a, nonzero_b, iter = 20, tol = 10^(-6), silent = FALSE, model = c("arima", "lme"), arformula = c(1,0,0), lmeformula = " ~ -1 + time + (1|id)")
{
  list.of.packages <- c("lme4", "forecast")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, library, character.only = TRUE)


  # checks
  if(ncol(B) <= max(nonzero_b)) {
    message("At least one of the nonzero options for B is not sparse. Changing to meet criteria")
    nonzeroB = nonzero_b[nonzero_b < (ncol(B) - 2)]
  }

  if(ncol(A) <= max(nonzero_a)) {
    message("At least one of the nonzero options for A is not sparse. Changing to meet criteria")
    nonzero_a = nonzero_a[nonzero_a < (ncol(A) - 2)]
  }


  #Create the matrix A
  alpha = sapply(nonzero_a, function(x) c(alphaInit))

  varTol1 = matrix(0, nrow = nrow(A), ncol = length(nonzero_a))
  varTol2 = matrix(0, nrow = nrow(B), ncol = length(nonzero_b))
  i = 0
  e = 10


  # format data
  id_a = A[,"id"]
  id_b = B[,"id"]
  time_a = A[,"time"]
  time_b = B[,"time"]

  A = as.matrix(A[, 3:ncol(A)])
  B = as.matrix(B[, 3:ncol(B)])

  while (e > tol & i <= iter) {
    i = i +1

    # refresh
    if(i > 1) varTol1 = gamma
    if(i > 1) varTol2 = zeta


    dist = sqrt(colSums(alpha^2))
    alpha = sweep(alpha, 2, dist, "/")

    gamma =  A %*% alpha
    # dist  = sqrt(colSums(gamma^2))
    # gamma = sweep(gamma, 2, dist, "/")
    gamma = as.matrix(scale_rm(data.frame(id = id_a, time = time_a, gamma = gamma))[,-c(1,2)])

    if(length(model) != 1 | !(model %in% c("arima", "lme"))) {
      model = "lme"
      lmeformula = " ~ -1 + time + (1|id)"

      print("Model not correctly specify. \n Default is lme with formula ~ -1 + time + (1|id)")
    }

    if(model == "arima") {
      me = list()
      pred_me = matrix(NA, nrow(B), 1)

      if(is.null(arformula)){
        for(n in unique(id_a)){
          me[[n]] = auto.arima(gamma[n == id_a],max.p = 5,max.q = 5,max.P = 5,max.Q = 5,max.d = 3,seasonal = FALSE,ic = 'aicc')
          pred_me[which(id_b ==n)] = as.numeric(forecast(me[[n]], h = length(time_b[which(id_b == n)]))$fitted)

        }
      } else {
        for(n in unique(id_a)){
          me[[n]] = arima(gamma[n == id_a], order = arformula, method = "ML")
          pred_me[which(id_b ==n)] = as.numeric(predict(me[[n]], n.ahead = length(time_b[which(id_b == n)]))$pred)

        }
      }

      # for(n in unique(id_a)){
      #   pred_me[which(id_b ==n)] = as.numeric(predict(me[[n]], n.ahead = length(time_b[which(id_b == n)]))$pred)
      # }
    }

    if(model == "lme") {
      me = sapply(1:ncol(alpha), function(j) lmer(as.formula(paste("gamma", lmeformula)), data = data.frame(gamma = gamma[,j], time = time_a, id = id_a), REML = TRUE))
      pred_me = sapply(1:ncol(alpha), function(j) predict(me[[j]], newdata = data.frame(time = time_b, id = id_b), allow.new.levels = TRUE, re.form = NULL))

    }

    beta = t(B) %*% pred_me #NAS IN PREDME

    rm(me, pred_me)

    beta = apply(rbind(beta,nonzero_b), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })

    dist = sqrt(colSums(beta^2))
    beta = sweep(beta, 2, dist, "/")

    zeta = B %*% beta
    # dist = sqrt(colSums(zeta^2))
    # zeta = sweep(zeta, 2, dist, "/")
    zeta = as.matrix(scale_rm(data.frame(id = id_b, time = time_b, zeta = zeta))[,-c(1,2)])


    if(model == "arima") {
      me = list()
      pred_me = matrix(NA, length(gamma), 1)

      if(is.null(arformula)){
        for(n in unique(id_b)){
          me[[n]] = auto.arima(zeta[n == id_b],max.p = 5,max.q = 5,max.P = 5,max.Q = 5,max.d = 3,seasonal = FALSE,ic = 'aicc')
          pred_me[which(id_a ==n)] = as.numeric(forecast(me[[n]], h = length(time_a[which(id_a == n)]))$fitted)

        }
      } else {
        for(n in unique(id_b)){
          me[[n]] = arima(zeta[n == id_b], order = arformula, method = "ML")
          pred_me[which(id_a ==n)] = as.numeric(predict(me[[n]], n.ahead = length(time_a[which(id_a == n)]))$pred)

        }
      }

      # for(n in unique(id_b)){
      #   pred_me[which(id_a ==n)] = as.numeric(predict(me[[n]], n.ahead = length(time_a[which(id_a == n)]))$pred)
      # }
    }

    if(model == "lme") {
      me = sapply(1:ncol(beta), function(j) lmer(as.formula(paste("zeta", lmeformula)), data = data.frame(zeta = zeta[,j], time = time_b, id = id_b), REML = TRUE))
      pred_me = sapply(1:ncol(beta), function(j) predict(me[[j]], newdata = data.frame(time = time_a, id = id_a), allow.new.levels = TRUE, re.form = NULL))

    }
    alpha = t(A) %*% pred_me

    rm(me, pred_me)

    alpha = apply(rbind(alpha,nonzero_a), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })


    if(length(nonzero_a) == 1) e = mean(abs(gamma - varTol1)) + mean(abs(zeta - varTol2))
    if(length(nonzero_a) > 1) e  = mean(colMeans(abs(gamma - varTol1))) + mean(colMeans(abs(zeta - varTol2)))

    textSCCA = paste0(" Common convergence error: ", round(e, 5), " & Iterations: ", i)
    if(isFALSE(silent) & (e<= tol || i > iter)) cat(textSCCA, "\r")

  }

  if(model == "lme") {
    me_x = sapply(1:ncol(alpha), function(j) lmer(as.formula(paste("gamma", lmeformula)), data = data.frame(gamma = gamma[,j], time = time_a, id = id_a), REML = TRUE))
    # pred_x = sapply(1:ncol(alpha), function(j) predict(me_x[[j]], newdata = data.frame(time = time_b, id = id_b), allow.new.levels = TRUE, re.form = NULL))
    # lmer(as.formula(paste("gamma", lmeformula)), data = data.frame(gamma = gamma, time = time_a, id = id_a), REML = TRUE)
    me_y = sapply(1:ncol(beta), function(j) lmer(as.formula(paste("zeta", lmeformula)), data = data.frame(zeta = zeta[,j], time = time_b, id = id_b), REML = TRUE))
    # pred_y = sapply(1:ncol(beta), function(j) predict(me_y[[j]], newdata = data.frame(time = time_a, id = id_a), allow.new.levels = TRUE, re.form = NULL))
    # me_y = lmer(as.formula(paste("zeta", lmeformula)), data = data.frame(zeta = zeta, time = time_b, id = id_b), REML = TRUE)
  } else {
    me_x = NULL
    me_y = NULL
  }

  return(list(a = alpha, b = beta, conv = e, iter = i, me_x = me_x, me_y = me_y))
}
