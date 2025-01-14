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

# create cv heatmap ------------------------------------------------------------
#' Plot heatmap of cv w.r.t. the penalty parameter perfotmance.
#'
#' This function stardadises matrices with multiple measurements
#' w.r.t. a chosen origin.
#'
#' @param mat A matrix.
#' @param palette Character. Name of a palette for the heatmap. Default is "Teal".
#' @export
myHeatmap=
  function (mat, palette = "Teal", coln = 12, xlab = "", ylab = "",
            axes = FALSE)
  {
    par(fig = c(0, 7/10, 0, 1))
    image(mat, col = hcl.colors(length(mat), palette, rev = TRUE),
          axes = axes, xlab = xlab, ylab = ylab)
    if (isFALSE(axes)) {
      axis(1, seq(0, 1, length = nrow(mat)), rownames(mat))
      axis(2, seq(0, 1, length = ncol(mat)), colnames(mat),
           las = 1)
    }}
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
