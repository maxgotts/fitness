#' @title landscape
#'
#' @description Creates a fitness landscape from a kernel density estimation
#'
#' @param dat A data set object
#' @param id Scope (e.g., country codes or individual IDs)
#' @param time Time (e.g., time periods are given by years, months, ...)
#'
#' @return A data frame object that contains a summary of a sample that
#'     can later be converted to a TeX output using \code{overview_print}
#' @examples
#' data(iris)
#' output_table <- overview_tab(dat = toydata, id = ccode, time = year)
#' @export
#' @importFrom ks "kde"
#'
#'

landscape <- landscapes <- function(x, de=NULL, eta=1, mu=(max(x)-min(x))*3/4, sigma=2, clean=FALSE, formula=NULL, save.parameters=TRUE, kde=FALSE, mean.scale=FALSE, verbose=TRUE) {
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

  # Determine if the stored x and mu are identical to new inputs
  # This will indicate whether HTH and HT need to be recalculated
  x.mu.existNotIdentical <- existsNotIdentical(".x.316d38f30e1cf94cd75afef3f05ffa50c659c48a",x) |
    existsNotIdentical(".mu.316d38f30e1cf94cd75afef3f05ffa50c659c48a",mu)

  # Calculate HTH
  if (!exists(".HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a") |
      existsDimensionsWrong(".HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a", de) |
      x.mu.existNotIdentical |
      clean
  ) { # If (1) stored HTH does not exist, (2) has the wrong dimensions,
    # (3) x or mu have changed, or (4) the user forces a clean run, calculate
    # HTH from raw parameters
    HTH <- generateHTH(x, mu)
    if (save.parameters) { # If storing variables option is selected, store HTH
      # Also store x and mu so we can check these later
      assign(".x.316d38f30e1cf94cd75afef3f05ffa50c659c48a", x, envir=.GlobalEnv)
      assign(".mu.316d38f30e1cf94cd75afef3f05ffa50c659c48a", mu, envir=.GlobalEnv)
      assign(".HTH.316d38f30e1cf94cd75afef3f05ffa50c659c48a", HTH, envir=.GlobalEnv)
    }
  } else { # If (1) HTH exists, (2) has the right dimensions, (3) is from the right
    # x and mu, and (4) the user has not forced a clean run, use the stored HTH
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
      x.mu.existNotIdentical |
      clean
  ) { # If (1) stored HT does not exist, (2) has the wrong dimensions,
    # (3) x or mu have changed, or (4) the user forces a clean run, calculate
    # HT from raw parameters
    HT <- generateHT(x, mu)
    if (save.parameters) { # If storing variables option is selected,save HT
      # We do not need to save x and mu because save HTH will have been triggered
      # above if this code is running
      assign(".HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a", HT, envir=.GlobalEnv)
    }
  } else { # If (1) HT exists, (2) has the right dimensions, (3) is from the right
    # x and mu, and (4) the user has not forced a clean run, use the stored HT
    HT <- .HT.316d38f30e1cf94cd75afef3f05ffa50c659c48a
  }

  # Calculate the diag(de) matrix
  diagde <- diag(de)

  # Calculate w using
  w <- rowSums(solve(diagde%*%HTH%*%diagde+eta*M)%*%diagde%*%HT%*%diagde)
  # `rowSums` sums up all of the entries in the kth row and makes that the kth
  # entry in a vector, which is about 500ms faster at len=1000 than multiplying
  # by the vector `rep(1,times=len)`

  # If the `mean.scale` option is set to true, both de and w will be scaled so
  # their means will be equal to 1
  if (mean.scale) {
    if (verbose) cat("Scaling de and w so their means will be equal to 1.\nIf this is unexpected, change `mean.scale=FALSE`.\n")
    de <- de/mean(de)
    w <- w/mean(w)
  }

  # Return the x, de, and w values
  return(list(
    x=x,
    de=de,
    w=w
  ))
}



data(iris)
len <- 100
x <- seq(4.5,7,length.out=len)
rawde <- iris[iris$Species=="virginica","Petal.Length"]
fl <- landscape(x, rawde, eta=1, mu=2, kde=TRUE, mean.scale=TRUE, verbose=FALSE)
plot(x,fl$de,type="n", xlab="Petal length (cm)", ylab="Value of phenotype")
lines(x,fl$de,col="red")
lines(x,fl$w,col="blue")
legend(6.2,1.5,legend=c("Density","Fitness"), col=c("red","blue"), lty=1, cex=0.8)






