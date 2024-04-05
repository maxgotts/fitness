#' @title landscape
#'
#' @description Creates a fitness landscape from a kernel density estimation
#'
#' @param de Vector of the phenotype observatoins
#' @param h2 Narrow-sense heritability estimate h^2
#' @param SE Narrow-sense heritability standard error
#' @param n Narrow-sense heritability sample size
#'
#' @return Value of epsilon, the standard deviation of
#' @examples
#' set.seed(42)
#' de <- rnorm(100,mean=0, sd=0.1)
#' h2 <- seq(0.5,1,length.out=100)
#' epsilon <- c()
#' for (k in 1:100) {
#'   epsilon[k] <- approximateEpsilon(de, h2=h2[k], SE=0.1, n=50)
#' }
#' plot(h2, epsilon, type="l")
#' @export
#'

approximateEpsilon <- function(de, h2, SE, n) {
  return(sqrt(n^2/(h2*(n-1))*SE^2-(1-h2)^2/h2*sd(de)^2))
}




