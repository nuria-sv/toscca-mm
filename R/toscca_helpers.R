# toscca helpers

# model summary ----------------------------------------------------------------
#' @export
print.toscca <- function(x, ...) {
  K <- sum(sapply(x, is.list))

  cat("\n toscca summary\n")
  cat("---------------------------------\n")
  cat("  Components:", K, "\n")
  cat("---------------------------------\n")

  if (K > 1) {
    alphas <- sapply(x, function(obj) sum(obj$alpha != 0))
    betas  <- sapply(x, function(obj) sum(obj$beta  != 0))

    cat(sprintf("  Supp(w1): %s\n", paste(alphas, collapse = ", ")))
    cat(sprintf("  Supp(w2): %s\n", paste(betas, collapse = ", ")))

  } else if (K == 1) {
    cat("  Supp(w1):", sum(x[[1]]$alpha!=0), "\n")
    cat("  Supp(w2):", sum(x[[1]]$beta!=0), "\n")

  } else {
    cat("  Empty object\n")
  }

  cat("---------------------------------\n")
}


# summary -----------------------------------------------------------------
#' @export
summary.tosccamm <- function(object, ...){
  K = sum(sapply(object, is.list))
  cat(K, "components: \n")

  if(K > 1) {
    for (k in 1:K) {
      cat("  K=",k, ":\n")

      cat("   Summary estimated model for mm in X: \n")
      print(summary(object[[k]]$me_x))
      # cat("\n")
      cat("   Summary estimated model for mm in Y: \n")
      print(summary(object[[k]]$me_y))
      cat("\n")
    }

  } else if (K == 1) {

    cat("  Summary estimated model for mm in X: \n")
    print(summary(object$me_x))
    cat("\n")
    cat("  Summary estimated model for mm in Y: \n")
    print(summary(object$me_y))
    cat("\n")

  } else {
    cat("Invalid object")
  }

  # cat("Supp(w1) =", x$alpha, "\n")
  # cat("Supp(w2) =", x$beta, "\n")
}

# plot -------------------------------------------------------------------------
#' @export
# plot cross validated penalty x
# plot selection stability     x
# cpev                         x
# plot clustring K1 K2
# plot latent paths
plot.toscca <- function(x, data_list = NULL, show=c("coefficients","cca-grid"),
                      Z=NULL, values=NULL, gather.plts = TRUE, mm = NULL, palette = "magma",
                      cent = 3,
                      groupsets=NULL, codataweights=FALSE, ...){
  if(is.null(mm)) stop("mm must be TRE/FALSE. is True for multiple measurements.")
  show_ggplot <- FALSE
  y<-NULL
  # if(requireNamespace("ggplot2")&requireNamespace("ggpubr")) show_ggplot <- TRUE
  # if (!require("RColorBrewer", character.only=T, quietly=T)) {
  #   install.packages("RColorBrewer")
  # }
  # library("RColorBrewer", character.only=T)
  # if (!require("viridis", character.only=T, quietly=T)) {
  #   install.packages("viridis")
  # }
  # library("viridis", character.only=T)
  # if (!require("gridExtra", character.only=T, quietly=T)) {
  #   install.packages("gridExtra")
  # }
  # library("gridExtra", character.only=T)
  # if (!require("ClusterR", character.only=T, quietly=T)) {
  #   install.packages("ClusterR")
  # }
  # library("ClusterR", character.only=T)
  # if (!require("cluster", character.only=T, quietly=T)) {
  #   install.packages("cluster")
  # }
  # library("cluster", character.only=T)
  K = sum(sapply(x, is.list))

  x.temp = x[[1]]

  # plt 1: best threshold and selection stability ------------------------------
  if(x.temp$cv_logical){

    print("Printing penalty cross-validation ans selection stability")

    # cat("  K=",k, ":\n")
     if(gather.plts) {
       mats = lapply(1:K, function(k) x[[k]]$mat_cc)
       plt_ls = lapply(mats, function(mat) myHeatmap(mat, palette = palette))
       layout_title <- grid::textGrob("cc w.r.t. sparsity levels", gp = grid::gpar(fontsize = 18, fontface = "bold", col = "black", family = "Helvetica"))

       gridExtra::grid.arrange(grobs = plt_ls, ncol = K, top = layout_title)

     } else {

      x.k = x[[1]]
      plt = myHeatmap(x.k$mat_cc, show_labels = T, K, palette)
      plt = plt  + ggplot2::ggtitle("cc w.r.t. sparsity levels") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", color = "black",
                                                                                          hjust = 0.5, vjust = 1.5,
                                                                                          lineheight = 1.2))
      print(plt)

       }
    if (!is.null(data_list)) {
      print(plt.selstab(x.k$mat_cc[c(1,2, 5, 7, 8), ], X = data_list[[1]], Y= data_list[[2]], mm = mm, palette = palette))
    } else {
      message("to plot selection stability provide data_list with respective matrices.")
    }


  }

  # plt 2: CPEV ----------------------------------------------------------------
  print("cpev plots")
  alpha = (sapply(1:K, function(k) x[[k]]$alpha))
  cpev_toscca = cumsum(sapply(1:K, function(k) cpev.toscca(data_list[[1]], alpha[,k])))
  cat("cpev is: ", cpev_toscca, "\n")
  if(K>1) {
    auto_cor = stats::cor(alpha)[2:(K+1)]
    adj_cpev_toscca = c(cpev_toscca[1],
                        sapply(2:K, function(k) cpev_toscca[k]*prod(1-abs(auto_cor[k-1:k]))))
    cat("adjusted cpev is: ", adj_cpev_toscca)

    df <- data.frame(
      K_index = 0:K,
      cpev_toscca = c(0,cpev_toscca),
      adj_cpev_toscca = c(0,adj_cpev_toscca)
    )
    custom_colors <- c("cpev" = "#5E177F",  # Bright Yellow
                       "adj cpev" = "#F0703C")  # Magenta-Purple
    linetypes <- c("cpev" = "solid",  # Bright Yellow
                       "adj cpev" = "dashed")
    p = ggplot2::ggplot(df, ggplot2::aes(x = K_index)) +
      ggplot2::geom_line(ggplot2::aes(y = cpev_toscca, color = "cpev", linetype = "cpev"), linewidth = 1) +
      ggplot2::geom_line(ggplot2::aes(y = adj_cpev_toscca, color = "adj cpev", linetype = "adj cpev"), linewidth = 1) +
      ggplot2::scale_color_manual(name = "Legend", values = custom_colors) +  # Manually set colors
      ggplot2::scale_linetype_manual(name = "Legend", values =linetypes) +  # Set different linetypes
      ggplot2::scale_x_continuous(breaks = seq(1, K, by = 1)) +  # Force integer breaks
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "K",
        y = "%"
      )
    p = p+   ggplot2::ggtitle("cpev toscca") + ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", color = "black",
                                                                                     hjust = 0.5, vjust = 1.5,
                                                                                     lineheight = 1.2))
    print(p)
  }

  # plt 3: latent space --------------------------------------------------------
  if(!mm) {
    print("score plots")
    K = dim(x$a)[2]
    n.plts = sum(1:K)
  }
  if(mm) {
    print("latent trajectory plots")

    # plt.cluster = function(lv){
    #   if(!is.data.frame(lv)) {
    #     lv = data.frame(lv)
    #     colnames(lv) <- paste0("K", 1:ncol(lv))
    #   }
    #
    #   p = ggplot2::ggplot(data = lv, mapping = aes(x = K1, y = K2,  fill = Cluster))  +
    #     ggplot2::stat_ellipse(type = "t", geom = "polygon", alpha = 0.2) +
    #     ggplot2::geom_point(data = lv,
    #                mapping = aes(x = K1, y = K2, color = Region, shape = Region)) +
    #     ggplot2::theme_bw() +  ggplot2::ggtitle(TeX("TOSCCA 100")) +
    #     ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 10), axis.title.x = ggplot2::element_blank())
    # }

    # nk = combn(1:K)
    # lv = data_list[[1]]%*%x$alpha
    # kmeans_toscca_lva<- lapply(1:ncol(nk), function(i) kmeans(lv[,c(nk[1,i], nk[2,i])], centers = cent, nstart = 2))
    # clustera = lapply(1:length(kmeans_toscca_lva), function(i) as.factor(kmeans_toscca_lva[[i]]$cluster))
    #
    # # par(mfrow = c(dim(nk)))
    # plt_ls = lapply(1:length(kmeans_toscca_lva), function(i) plt.cluster(kmeans_toscca_lva[[i]]))
    # layout_title <- textGrob("cc w.r.t. sparsity levels", gp = grid::gpar(fontsize = 18, fontface = "bold", col = "black", family = "Helvetica"))
    # grid.arrange(grobs = plt_ls, ncol = ncol(nk), top = layout_title)
    #
    # lv = data_list[[2]]%*%x$beta
    # kmeans_toscca_lvn<- lapply(1:ncol(nk), function(i) kmeans(lv[,c(nk[1,i], nk[2,i])], centers = cent, nstart = 2))
  }
  if(length(show)==2) show <- "coefficients"
  if(show=="coefficients"){
    if(show_ggplot){
      p1 <- ggplot2::ggplot(data.frame(x=x$sigmahat/x$penalties, y=x$beta^2))+
        ggplot2::aes(x=x,y=y)+
        ggplot2::geom_point()+
        ggplot2::labs(x="Prior variance", y="Squared regression coefficients")+
        ggplot2::geom_abline(intercept=0,slope=1)
      p1
      return(p1)
    }
  #   par(mfrow=c(1,1))
  #   #par(mfrow=c(1,3))
  #   #plot(log(x$penalties), x$beta)
  #   #lot(x$sigmahat/x$penalties, x$beta)
  #   plot(x$sigmahat/x$penalties, x$beta^2,
  #        xlab= "Prior variance", ylab="Squared regression coefficients")
  #   abline(0,1)
  # }else if(show=="priorweights"){
  #   if(is.null(Z)&&is.null(groupsets)){
  #     stop("Provide Z or groupsets used to fit ecpc")
  #   }else if(!is.null(Z)&!is.null(groupsets)){
  #     stop("Provide either Z or groupsets, not both")
  #   }else{
  #     if(!is.null(groupsets)) Z <- lapply(groupsets,createZforGroupset)
  #     if(is.null(values)){
  #       print("Prior variance contribution per group and co-data source is plotted.")
  #       print("Note that for splines, the original continuous co-data values should be provided in values")
  #     }else{
  #       print("Prior variance contribution per group/continuous values and co-data source is plotted.")
  #     }
  #     if(is.null(names(Z))) names(Z) <- paste("Co-data set",1:length(Z))
  #
  #     df <- data.frame()
  #     p1 <- list()
  #     nrows <- floor(sqrt(length(Z)))
  #     ncols <- ceiling(length(Z)/nrows)
  #     par(mfrow=c(nrows,ncols))
  #     for(g in 1:length(Z)){
  #       ind_g <- attributes(x$gamma)$codataSource==g
  #       if(is.null(values[[g]])){
  #         if(is.null(colnames(Z[[g]])) || any(colnames(Z[[g]])=="")){
  #           if(codataweights){
  #             if(show_ggplot){
  #               temp <- data.frame(x=factor(1:dim(Z[[g]])[2]),
  #                                  y=x$gamma[ind_g]*x$tauglobal* x$w[g],
  #                                  Codatasource=names(Z)[g])
  #               #df <- rbind(df,temp)
  #               p1[[g]] <- ggplot2::ggplot(temp)+ggplot2::aes(x=x,y=y)+
  #                 ggplot2::geom_point()+
  #                 ggplot2::labs(x="Co-data variable",y="Prior variance weight", title=names(Z)[g])
  #             }else{
  #               plot(factor(1:dim(Z[[g]])[2]), x$gamma[ind_g]*x$tauglobal* x$w[g],
  #                    main=names(Z)[g], xlab="Co-data variable",
  #                    ylab="Prior variance weight")
  #             }
  #           }else{
  #             if(show_ggplot){
  #               temp <- data.frame(x=factor(1:dim(Z[[g]])[2]),
  #                                  y=x$gamma[ind_g]*x$tauglobal,
  #                                  Codatasource=names(Z)[g])
  #               #df <- rbind(df,temp)
  #               p1[[g]] <- ggplot2::ggplot(temp)+ggplot2::aes(x=x,y=y)+
  #                 ggplot2::geom_point()+
  #                 ggplot2::labs(x="Co-data variable",y="Prior variance weight", title=names(Z)[g])
  #             }else{
  #               plot(factor(1:dim(Z[[g]])[2]), x$gamma[ind_g]*x$tauglobal,
  #                    main=names(Z)[g], xlab="Co-data variable",
  #                    ylab="Prior variance weight")
  #             }
  #           }
  #         }else{
  #           if(codataweights){
  #             if(show_ggplot){
  #               temp <- data.frame(x=factor(colnames(Z[[g]])),
  #                                  y=x$gamma[ind_g]*x$tauglobal* x$w[g],
  #                                  Codatasource=names(Z)[g])
  #               #df <- rbind(df,temp)
  #               p1[[g]] <- ggplot2::ggplot(temp)+ggplot2::aes(x=x,y=y)+
  #                 ggplot2::geom_point()+
  #                 ggplot2::labs(x="Co-data variable",y="Prior variance weight", title=names(Z)[g])
  #             }else{
  #               plot(factor(colnames(Z[[g]]), levels=colnames(Z[[g]])),
  #                    x$gamma[ind_g] * x$tauglobal * x$w[g],
  #                    main=names(Z)[g], xlab="Co-data variable",
  #                    ylab="Prior variance weight")
  #             }
  #           }else{
  #             if(show_ggplot){
  #               temp <- data.frame(x=factor(colnames(Z[[g]])),
  #                                  y=x$gamma[ind_g]*x$tauglobal,
  #                                  Codatasource=names(Z)[g])
  #               #df <- rbind(df,temp)
  #               p1[[g]] <- ggplot2::ggplot(temp)+ggplot2::aes(x=x,y=y)+
  #                 ggplot2::geom_point()+
  #                 ggplot2::labs(x="Co-data variable",y="Prior variance weight", title=names(Z)[g])
  #             }else{
  #               plot(factor(colnames(Z[[g]]), levels=colnames(Z[[g]])),
  #                    x$gamma[ind_g] * x$tauglobal,
  #                    main=names(Z)[g], xlab="Co-data variable",
  #                    ylab="Prior variance weight")
  #             }
  #           }
  #         }
  #       }else{
  #         if(codataweights){
  #           if(show_ggplot){
  #             temp <- data.frame(x=values[[g]],
  #                                y= Z[[g]]%*%x$gamma[ind_g] * x$tauglobal * x$w[g],
  #                                Codatasource=names(Z)[g])
  #             #df <- rbind(df,temp)
  #             p1[[g]] <- ggplot2::ggplot(temp)+ggplot2::aes(x=x,y=y)+
  #               ggplot2::geom_line()+
  #               ggplot2::labs(x="Continuous co-data variable",y="Prior variance weight", title=names(Z)[g])
  #           }else{
  #             plot(values[[g]], Z[[g]]%*%x$gamma[ind_g] * x$tauglobal * x$w[g],
  #                  main=names(Z)[g], xlab="Continuous co-data variable",
  #                  ylab="Prior variance weight")
  #           }
  #         }else{
  #           if(show_ggplot){
  #             temp <- data.frame(x=values[[g]],
  #                                y= Z[[g]]%*%x$gamma[ind_g] * x$tauglobal,
  #                                Codatasource=names(Z)[g])
  #             #df <- rbind(df,temp)
  #             p1[[g]] <- ggplot2::ggplot(temp)+ggplot2::aes(x=x,y=y)+
  #               ggplot2::geom_line()+
  #               ggplot2::labs(x="Continuous co-data variable",y="Prior variance weight", title=names(Z)[g])
  #           }else{
  #             plot(values[[g]], Z[[g]]%*%x$gamma[ind_g] * x$tauglobal,
  #                  main=names(Z)[g], xlab="Continuous co-data variable",
  #                  ylab="Prior variance weight")
  #           }
  #         }
  #       }
  #     }
  #
  #     if(show_ggplot){
  #       #df$Codatasource <- factor(df$Codatasource, levels=unique(df$Codatasource),
  #       #                          labels=unique(df$Codatasource))
  #       ylims <- sapply(1:length(Z), function(g){
  #         ind_g <- attributes(x$gamma)$codataSource==g
  #         if(is.null(values[[g]])){
  #           if(codataweights){
  #             return(range(x$gamma[ind_g] * x$tauglobal *x$w[g]))
  #           }else{
  #             return(range(x$gamma[ind_g]* x$tauglobal))
  #           }
  #         }else{
  #           if(codataweights){
  #             return(range(Z[[g]]%*%x$gamma[ind_g] * x$tauglobal *x$w[g]))
  #           }else{
  #             return(range(Z[[g]]%*%x$gamma[ind_g]* x$tauglobal))
  #           }
  #         }
  #       })
  #       ylims <- range(c(ylims))
  #       #ylims <- ylims + 0.01*diff(ylims)*c(-1,1)
  #       for(g in 1:length(Z)){
  #         p1[[g]] <- p1[[g]]+ggplot2::ylim(ylims)
  #       }
  #       p2 <- ggpubr::ggarrange(plotlist=p1, nrow=nrows)
  #       p2
  #       return(p2)
  #     }
  #   }
  # }
  }
}
