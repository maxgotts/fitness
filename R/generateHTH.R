generateHTH <- function(x, mu) {
  len <- length(x)
  gauss <- dnorm(x, mean=(min(x)+max(x))/2, sd=mu)
  hauss <- as.vector(matrix(nrow=1,ncol=len))
  for (k in 1:len) hauss[k] <- sum(gauss*gauss[wrap(len-k+1,len)])
  HTH <- matrix(NA, nrow=len, ncol=len)
  for (k in 1:len) {
    HTH[,k] <- hauss[wrap(len-k+1,len)]
  }
  return(HTH)
}
