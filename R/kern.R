#'@importFrom stats dnorm pnorm
kern <- function(kernel = "epanechnikov") {

  ## CHOOSE KERNEL
  #  input:   kernel .. one of:'epan'  ...Epanechnikov (parabolic)
  #                            'gauss' ...Gaussian
  #                            'quart' ...Quartic (biweight)
  #                            'rect'  ...Rectangular (uniform)
  #                            'trian' ...Triangular
  #
  #  outputs: K .. kernel function, type function
  #           W .. definite integral of K on interval (-supp, x)->function of x, type function
  #           details .. supp  .. support interval [-supp,supp], where beta and var are calculated
  #                      beta  .. value of k-th moment of kernel
  #                      var   .. variation of kernel
  #                      k     .. order of kernel
  #                      name  .. name od kernel in function
  #                      coef  .. coefficients
  ##


  if (kernel== "epanechnikov"){
    K <- function(x) {3/4 * (1-x^2)*(-1<=x & x<=1)}
    W <- function(x) {1/4*(-x^3+3*x+2)*(-1<=x & x<=1)+(x>1)}
    details <- list(supp=1, beta=1/5, var=3/5, k=2, name="epan",
                    coef=c(-0.75,0,0.75))
  }
  else if (kernel == "gaussian"){
   # K <- function(x) {1/sqrt(2*pi) * exp(-1/2*x^2)}
   # W <- function(x) {1/2* erf(x/sqrt(2))+1/2} # library: pracma
    K<-function(x) {dnorm(x)}
    W<-function(x) {pnorm(x)}
    details <- list(supp=Inf, beta=1, var=1/(2*sqrt(pi)), k=2, name="gauss",
                    coef=1)
  }
  else if (kernel == "rectangular"){
    K <- function(x) {1/2*(-1<=x & x<=1)}
    W <- function(x) {(1/2 * x+1/2)*(-1<=x & x<=1)+(x>1)}
    details <- list(supp=1, beta=1/3, var=1/2, k=2, name="rect",
                    coef=0.5)
  }
  #else if (kernel == "triangular"){
  #  K <- function(x) {(1-abs(x))*(-1<=x & x<=1)}
  #  W <- function(x) {(x-abs(x)*x/2+1/2)*(-1<=x & x<=1) + (x>1)}
  #  details <- list(supp=1, beta=1/6, var=2/3, k=2, name="trian",
  #                  coef=NA)
  #}
  else if (kernel == "quartic"){
    K <- function(x) {(15/16 * (1-x^2)^2)*(-1<=x & x<=1)}
    W <- function(x) {(1/16 * (3*x^5-10*x^3+15*x)+1/2)*(-1<=x & x<=1)+(x>1)}
    details <- list(supp=1, beta=1/7, var=5/7, k=2, name="qvart",
                    coef=c(0.9375, 0, -1.8750, 0, 0.9375))
  }
  else {stop("kernel have to be one of: 'epanechnikov','gaussian','rectangular','quartic'")}

  return(list(K=K, W=W, details=details))
}
