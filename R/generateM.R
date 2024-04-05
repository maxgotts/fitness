generateM <- function(sigma, len) {
  M <- matrix(0, nrow=len, ncol=len)
  for (r in 1:len) {
    for (c in 1:r) {
      etheta <- exp(sigma/(len-1)-2*pi*1i*(r-c)/len)
      M[r,c] <- M[c,r] <- Re((1-etheta^len)/(1-etheta))
    }
  }
  return(M)
}
