require("pracma")

landConv <- landscapeByConvolutionTheorem <- function(x, de, epsilon) {
  gaussian <- dnorm(x, mean=(max(x)+min(x))/2, sd=epsilon)
  wde <- ifft(fft(de)/fft(gaussian))
  w <- wde/de
  return(w)
}
