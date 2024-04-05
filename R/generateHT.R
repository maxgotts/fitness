generateHT <- function(x, mu) { # make transpose!
  len <- length(x)
  gauss <- dnorm(x, mean=(min(x)+max(x))/2, sd=mu)
  HT <- matrix(NA, nrow=len, ncol=len)
  for (k in 1:len) { # 1, 10, 9, 8
    index <- c(1,len:2)[k]
    HT[index,] <- gauss[wrap(k,len)]
  }
  return(HT)
}
