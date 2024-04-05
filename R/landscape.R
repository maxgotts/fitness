#' @title landscape
#'
#' @description Creates a fitness landscape from a kernel density estimation
#'
#' @param x Vector of the range of the phenotypic axis of interest
#' @param de Vector of phenotype observations, of density of observations over axis x
#' @param eta Tuning parameter of weight of high-frequency penalty
#' @param epsilon Measure of heritability of phenotype (standard deviation of offpsring by parent)
#' @param sigma Tuning parmeter for intensity of high-frequency penalty
#' @param kde Boolean for if a kernel density estimation should be applied to de or not
#' @param clean Boolean; if true, will force program to re-generate auxilary matrices, even if they are appropriate to use
#' @param save.parameters Boolean; if true, will save the auxiliary matrices for reuse when appropriate
#' @param scale Method for rescaling ("mean","stretch") for comparison of frequency and fitness, or NULL
#' @param gg Boolean; produce ggplot of de and w, or not
#' @param bootstraps If NULL, perform no bootstraps, otherwise perform the number of bootstraps. Currently requires that de be a vector of phenotype observations, as resampling from phenotype density is not supported.
#' @param verbose Boolean; if false, warnings and updates will be suppressed
#'
#' @return If gg is false, then it returns a dataframe with x, de, and w. If gg is true, it returns a list with data the aforementioned dataframe, and gg, a ggplot of de and w
#' @examples
#' data(iris)
#' len <- 100
#' x <- seq(4.5,7,length.out=len)
#' rawde <- iris[iris$Species=="virginica","Petal.Length"]
#' fl <- fitness::landscape(x, rawde, eta=100, epsilon=sd(rawde), kde=TRUE, scale="stretch", verbose=FALSE)
#' plot(x,fl$de,type="n", xlab="Petal length (cm)", ylab="Value of phenotype")
#' lines(x,fl$de,col="red")
#' lines(x,fl$w,col="blue")
#' legend(6.2,1.5,legend=c("Density","Fitness"), col=c("red","blue"), lty=1, cex=0.8)
#' cat("Note: n=50 phenotypic observations is not necessarily sufficient to produce trustworthy fitness landscapes using this method.\n")
#' @export
#' @importFrom ks "kde"
#' @importFrom reshape2 "melt"
#' @importFrom dplyr "%>%", "mutate"
#' @importFrom ggplot2 "ggplot", "geom_line", "scale_color_manual", "facet_wrap", "theme", "strip.background", "labs", "theme_classic"
#'
#'

landscape <- function(x, de, epsilon, eta=1, sigma=2, kde=FALSE, bootstraps=NULL, clean=FALSE, save.parameters=TRUE, scale=NULL, gg=FALSE, verbose=TRUE, seed=NULL) {
  if (!is.null(bootstraps)) {
    if (!kde) {
      cat("Error: must provide de as raw phenotypic observations, not as density plot.\nFuture versions of `fitness` may include resampling from a provided distribution, but currently this feature is unsupported.")
      return(data.frame(x=x, de=de))
    }
    if (!is.null(seed)) set.seed(seed)
    if (!save.parameters) { # Store old parameters
      if (exists(".x.316d38f30e1cf94cd75afef3f05ffa50c659c48a"))
        x.save <- .x.316d38f30e1cf94cd75afef3f05ffa50c659c48a
      if (exists(".epsilon.316d38f30e1cf94cd75afef3f05ffa50c659c48a"))
        epsilon.save <- .epsilon.316d38f30e1cf94cd75afef3f05ffa50c659c48a
      if (exists(".HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a"))
        HTH.save <- .HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a
      if (exists(".HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a"))
        HT.save <- .HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a
      if (exists(".M.316d38f30e1cf94cd75afef3f05ffa50c659c48a"))
        M.save <- .M.316d38f30e1cf94cd75afef3f05ffa50c659c48a
      if (exists(".sigma.316d38f30e1cf94cd75afef3f05ffa50c659c48a"))
        sigma.save <- .sigma.316d38f30e1cf94cd75afef3f05ffa50c659c48a
    }
    fls <- as.data.frame(matrix(NA, ncol=4, nrow=0))
    colnames(fls) <- c("x","de","w","Bootstrap")
    fl.real <- landscape(x, de, eta=eta, epsilon=epsilon, sigma=sigma, kde=TRUE, scale=scale, clean=clean, verbose=verbose, save.parameters=T, gg=F)
    for (boot in 1:bootstraps) {
      resampled.de <- sample(de, size=length(de), replace=TRUE)
      fl.boot <- landscape(x, de=resampled.de, eta=eta, epsilon=epsilon, sigma=sigma, kde=TRUE, bootstraps=NULL, clean=FALSE, save.parameters=TRUE, scale=scale, gg=FALSE, verbose=verbose)
      fl.boot$Bootstrap <- boot
      fls <- rbind(fls, fl.boot)
    }
    if (!save.parameters) { # Replace new parameters with old ones
      if (exists("x.save"))
        .x.316d38f30e1cf94cd75afef3f05ffa50c659c48a <- x.save
      if (exists("epsilon.save"))
        .epsilon.316d38f30e1cf94cd75afef3f05ffa50c659c48a <- epsilon.save
      if (exists("HTH.save"))
        .HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a <- HTH.save
      if (exists("HT.save"))
        .HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a <- HT.save
      if (exists("M.save"))
        .M.316d38f30e1cf94cd75afef3f05ffa50c659c48a <- M.save
      if (exists("sigma.save"))
        .sigma.316d38f30e1cf94cd75afef3f05ffa50c659c48a <- sigma.save
    }
    if (!gg) {
      return(fls)
    } else {
      linewidth <- 1.8/(1+exp(bootstraps/75))+0.1 # logistic
      alpha <- 1.6/(1+exp(bootstraps/75))+0.2 # logisic
      return(list(
        data=rbind(fl.real %>% mutate(Bootstrap="Real"),fls),
        gg=fls %>%
          melt(id=c("x","Bootstrap")) %>%
          mutate(
            variable=ifelse(variable=="de","Distribution",ifelse(variable=="w","Fitness",NA)),
            variable=factor(variable, levels=c("Distribution","Fitness")),
          ) %>%
          ggplot(aes(x=x, y=value, color=variable, fill=as.factor(Bootstrap)))+
          geom_line(linewidth=linewidth, alpha=alpha)+
          scale_color_manual(values=c("#FE654F","#47B9FF"))+
          facet_wrap(~variable, ncol=2)+
          theme_classic()+
          theme(legend.position="none",
                strip.background=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank())+
          labs(x="", y="",color="")+
          geom_line(data=fl.real %>%
                      melt(id="x") %>%
                      mutate(
                        variable=ifelse(variable=="de","Distribution",ifelse(variable=="w","Fitness",NA)),
                        variable=factor(variable, levels=c("Distribution","Fitness")),
                      ),
                    aes(x=x, y=value, color=variable),
                    color="grey4", linetype="dashed",
                    inherit.aes=FALSE)
      ))
    }
  } else {
    # If `x` is specified as NULL, then generate a 100-length vector from
    # min to max of the observed phenotype values
    if (is.null(x))
      x <- seq(min(de), max(de), length.out=100)

    # Perform kernel density estimation, if `kde` option selected
    if (kde) {
      if (verbose) cat("Warning: performing kernel density estimation (KDE) on input vector.\nIf this is a mistake, change `kde=FALSE`.\n")
      de <- predict(ks::kde(de), x=x)
    }

    # Catch error if density data is not length of x-axis
    if (length(x)!=length(de)) {
      cat("Error: mismatched lengths of x and de.\nIf your data has not been converted to density data using KDE, change `kde=TRUE`.\nOutput is list of x and de.\n")
      return(list(x=x,de=de))
    }

    # Determine if the stored x and epsilon are identical to new inputs
    # This will indicate whether HTH and HT need to be recalculated
    x.epsilon.existNotIdentical <- existsNotIdentical(".x.316d38f30e1cf94cd75afef3f05ffa50c659c48a",x) |
      existsNotIdentical(".epsilon.316d38f30e1cf94cd75afef3f05ffa50c659c48a",epsilon)

    # Calculate HTH
    if (!exists(".HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a") |
        existsDimensionsWrong(".HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a", de) |
        x.epsilon.existNotIdentical |
        clean
    ) { # If (1) stored HTH does not exist, (2) has the wrong dimensions,
      # (3) x or epsilon have changed, or (4) the user forces a clean run, calculate
      # HTH from raw parameters
      HTH <- generateHTH(x, epsilon)
      if (save.parameters) { # If storing variables option is selected, store HTH
        # Also store x and epsilon so we can check these later
        assign(".x.316d38f30e1cf94cd75afef3f05ffa50c659c48a", x, envir=.GlobalEnv)
        assign(".epsilon.316d38f30e1cf94cd75afef3f05ffa50c659c48a", epsilon, envir=.GlobalEnv)
        assign(".HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a", HTH, envir=.GlobalEnv)
      }
    } else { # If (1) HTH exists, (2) has the right dimensions, (3) is from the right
      # x and epsilon, and (4) the user has not forced a clean run, use the stored HTH
      HTH <- .HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a
    }

    # Calculate M
    if (!exists(".M.316d38f30e1cf94cd75afef3f05ffa50c659c48a") |
        existsDimensionsWrong(".M.316d38f30e1cf94cd75afef3f05ffa50c659c48a", de) |
        existsNotIdentical(".sigma.316d38f30e1cf94cd75afef3f05ffa50c659c48a",x) |
        clean
    ) { # If (1) stored M does not exist, (2) has the wrong dimensions,
      # or (3) the user forces a clean run, calculate M from raw parameters
      M <- generateM(sigma, length(de))
      if (save.parameters) { # If storing variables option is selected, store M
        assign(".sigma.316d38f30e1cf94cd75afef3f05ffa50c659c48a", sigma, envir=.GlobalEnv)
        assign(".M.316d38f30e1cf94cd75afef3f05ffa50c659c48a", M, envir=.GlobalEnv)
      }
    } else { # If (1) M exists, (2) has the right dimensions, and (3) the user
      # has not forced a clean run, use the stored M
      M <- .M.316d38f30e1cf94cd75afef3f05ffa50c659c48a
    }

    # Calculate HT
    if (!exists(".HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a") |
        existsDimensionsWrong(".HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a", de) |
        x.epsilon.existNotIdentical |
        clean
    ) { # If (1) stored HT does not exist, (2) has the wrong dimensions,
      # (3) x or epsilon have changed, or (4) the user forces a clean run, calculate
      # HT from raw parameters
      HT <- generateHT(x, epsilon)
      if (save.parameters) { # If storing variables option is selected,save HT
        # We do not need to save x and epsilon because save HTH will have been triggered
        # above if this code is running
        assign(".HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a", HT, envir=.GlobalEnv)
      }
    } else { # If (1) HT exists, (2) has the right dimensions, (3) is from the right
      # x and epsilon, and (4) the user has not forced a clean run, use the stored HT
      HT <- .HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a
    }

    # Calculate the diag(de) matrix
    diagde <- diag(de)

    # Calculate w using
    w <- rowSums(solve(diagde%*%HTH%*%diagde+eta*M)%*%diagde%*%HT%*%diagde)
    # w <- rowSums(solve(M)%*%diagde%*%HT%*%diagde)
    # `rowSums` sums up all of the entries in the kth row and makes that the kth
    # entry in a vector, which is about 500ms faster at len=1000 than multiplying
    # by the vector `rep(1,times=len)`

    if (is.null(scale)) {
    } else if (scale=="mean") { # If the `scale` option is set to "mean", both de and w will be scaled so
      # their means will be equal to 1
      if (verbose) cat("Warning: Scaling de and w so their means will be equal to 1.\nIf this is unexpected, change `scale=NULL`.\n")
      de <- de/mean(de)
      w <- w/mean(w)
    } else if (scale=="stretch") { # If the `scale` option is set to "stretch",
      # w will be shifted and scaled to better compare the shape of density and fitness
      if (verbose) cat("Warning: Stretching (scaling and shifting) w so its shape is comparable with de.\nIf this is unexpected, change `scale=NULL`.\n")
      w <- max(de)*(w-min(w))/(max(w)-min(w))
    }

    # Return the x, de, and w values
    if (gg==FALSE) {
      return(data.frame(
        x=x,
        de=de,
        w=w
      ))
    } else {
      return(list(data=data.frame(
        x=x,
        de=de,
        w=w
      ),gg=data.frame(x=x, de=de, w=w) %>%
          melt(id=c("x")) %>%
          mutate(
            variable=ifelse(variable=="de","Distribution",ifelse(variable=="w","Fitness",NA)),
            variable=factor(variable, levels=c("Distribution","Fitness")),
          ) %>%
          ggplot(aes(x=x, y=value, color=variable))+
          geom_line(linewidth=0.8)+
          scale_color_manual(values=c("#FE654F","#47B9FF"))+
          theme(legend.position="top")+
          labs(x="", y="",color="")+
          theme_classic()+
          theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank())
      ))
    }
  }
}






