#'@description Kernel estimate of conditional hazard function for right-censored data with one covariate. Options include two methods for bandwidth selection.
#'@title Kernel estimate of conditional hazard function for right-censored data
#'@param times vector of observed times
#'@param delta vector of censoring indicator. 0 - censored, 1 - uncensored (dead)
#'@param covariate vector of covariate
#'@param h bandwidth vector of length 2, first element is bandwidth for time and second for covariate. If missing, h is found using some bandwidth selection method.
#'@param t vector of time points at which estimate is evaluated
#'@param x vector of covariate points at which estimate is evaluated
#'@param tx data frame of t and x at which estimate is evaluated
#'@param t.length number of grid points of time
#'@param x.length number of grid points of covariate
#'@param tmin,tmax minimum/maximum values for grid of time
#'@param xmin,xmax minimum/maximum values for grid of covariate
#'@param kernel kernel function, possible values are: "epanechnikov" (default), "gaussian", "rectangular", "quartic".
#'@param type Type of kernel estimate. Possible types are:  "exterior", "interior" (default).
#'@param type.w Type of weights. Default are Nadaraya-Watson weights.
#'@param parallel  allows parallel computation. Default is FALSE.
#'@param h.method method for bandwidth selection. Possible methods are: "crossval" (default), "maxlike".
#'@param optim.method method for numerical optimization of the crossvalidation or log-likelihood function. Possible methods are: "ga"(default).
#'@param tol the desired accuracy of optimization algorithm
#'@param run the number of consecutive generations without any improvement in the best fitness value before the GA is stopped.
#'@param ... additional arguments of GA algorithm
#'@return Returns an object of class 'khazardcond' which is a list with fields
#'\item{time.points}{vector of time points at which estimate is evaluated}
#'\item{covariate.points}{vector of covariate points at which estimate is evaluated}
#'\item{hazard}{matrix of hazard function values on grid or data.frame of time and covariate points and appropriate hazard values if hx is defined}
#'\item{h}{bandwidth vector}
#'\item{CVML}{value of crossvalidation or log-likelihood at h}
#'\item{details}{description of used methods}
#'\item{GA.result}{output of ga, object of class ga-class}
#'@seealso \code{\link{plot.khazardcond}}, \code{\link[GA]{ga}}
#'@examples
#'library(survival)
#'fit<-khazardcond(times = lung$time,delta = lung$status-1,covariate = lung$age,h=c(200,20))
#'@details External type of kernel estimator is defined as the ratio of kernel estimator of the conditional subdensity of the uncensored observations to the conditional survival function of the observable time. Internal type of kernel estimator is based on a convolution of the kernel function with a nonparametric estimator of the cumulative conditional hazard function.
#'@references Selingerova, I., Dolezelova, H., Horova, I., Katina, S., and Zelinka, J. (2016). Survival of Patients with Primary Brain Tumors: Comparison of Two Statistical Approaches. PloS one, 11(2), e0148733.
#'@import foreach
#'@import doParallel
#'@import GA
#'@export
khazardcond <-
  function(times,
           delta,
           covariate,
           h = NULL,
           t = NULL,
           x = NULL,
           tx=NULL,
           t.length = 100,
           x.length=100,
           tmin = NULL,
           tmax = NULL,
           xmin = NULL,
           xmax =NULL,
           kernel = "epanechnikov",
           type = "interior",
           type.w = "nw",
           parallel = FALSE,
           #value = "CVML",
           h.method = "crossval",
           optim.method = "ga",
           tol = ifelse(h.method == "crossval", 10 ^ (-6), 1),
           run = 2,
           ...) {

    possible_kernels <-
      c('epanechnikov', 'gaussian', 'rectangular', 'quartic')
    kernel <-
      possible_kernels[pmatch(tolower(kernel), possible_kernels)]
    if (is.na(kernel)) {
      stop("kernel have to be one of: 'epanechnikov','gaussian','rectangular','quartic'")
    }

    possible_type <- c('exterior', 'interior')
    type <- possible_type[pmatch(tolower(type), possible_type)]
    if (is.na(type)) {
      stop("type have to be one of: 'exterior','interior'")
    }

    possible_type.w <- c('nw')
    type.w <- possible_type.w[pmatch(tolower(type.w), possible_type.w)]
    if (is.na(type.w)) {
      stop("type.w have to be one of: 'nw'")
    }

    # possible_value <- c('cvml', 'hazard', 'both')
    # value <- possible_value[pmatch(tolower(value), possible_value)]
    # if (value == 'cvml') {
    #   value = 'CVML'
    # }
    # if (is.na(value)) {
    #   stop("value have to be one of: 'CVML','hazard','both'")
    # }

    possible_h.method <- c('crossval', 'maxlike')
    h.method <-
      possible_h.method[pmatch(tolower(h.method), possible_h.method)]
    if (is.na(h.method)) {
      stop("h.method have to be one of: 'crossval','maxlike'")
    }

    possible_optim.method <- c('ga')
    optim.method <-
      possible_optim.method[pmatch(tolower(optim.method), possible_optim.method)]
    if (is.na(optim.method)) {
      stop("optim.method have to be one of: 'ga'")
    }

    if (missing(times) | missing(delta)|missing(covariate))
      stop("Parameter times, delta or covariate is missing")
    if ((length(times)+ length(delta)+length(covariate))!=3*length(delta)) {
      stop('times, delta, and covariate have to be of same length')
    }

    notNA <- !(is.na(times) | is.na(delta)|is.na(covariate))
    times <- times[notNA]
    delta <- delta[notNA]
    covariate <- covariate[notNA]
    if (!(is.numeric(times) &
          sum(times < 0 |
              times == Inf) == 0)) {
      stop('times have to be nonnegative and finite numeric')
    }
    if (!(is.numeric(delta) &
          sum(delta == 0 |
              delta == 1) == length(delta))) {
      stop('delta have to contains only 0 or 1')
    }
    if (!(is.numeric(covariate) &
          sum(covariate == Inf) == 0)) {
      stop('covariate have to be finite numeric')
    }

    ##
    if (!is.null(tx)) {
      if (dim(tx)[2]!=2|(!is.numeric(tx[,1])|!is.numeric(tx[,2]))){
        stop('tx have to be matrix or data frame with 2 numeric columns.')

      }
      x<-tx[,2]
      t<-tx[,1]

      t <- na.omit(t)
      if (!(is.numeric(t) &
            sum(t < 0 |
                t == Inf) == 0)) {
        stop('First column of tx have to be nonnegative and finite numeric.')
      }
      x <- na.omit(x)
      if (!(is.numeric(x) &
            sum(  x == Inf) == 0)) {
        stop('Second column of tx have to be finite numeric')
      }
      if(! length(t)==length(x)){stop('Both columns have to have same length after omiting Nas.')}
    }
    ##



    if (is.null(t)&is.null(tx)) {
      if (is.null(tmin)) {
        tmin = 0
      }
      else if (!is.numeric(tmin) |
               length(tmin) > 1 |
               min(tmin) < 0 |
               max(tmin) == Inf) {
        stop('tmin have to be nonnegative finite number')
      }
      if (is.null(tmax)) {
        tmax = max(times)
      }
      else if (!is.numeric(tmax) |
               length(tmax) > 1 |
               max(tmin) > min(tmax) |
               max(tmax) == Inf) {
        stop('tmax have to be nonnegative finite number bigger than tmin, if indicated')
      }
      if (!is.numeric(t.length) |
          length(t.length) > 1 |
          min(t.length) < 1 |
          max(t.length) == Inf) {
        stop('t.length have to be possitive integer')
      }
      if (is.numeric(t.length) & length(t.length) == 1) {
        if (abs(t.length - round(t.length)) > 0) {
          stop('t.length have to be possitive integer')
        }
      }
      t <- seq(tmin, tmax, length.out = t.length)
    }
    t <- na.omit(t)
    if (!(is.numeric(t) &
          sum(t < 0 |
              t == Inf) == 0)) {
      stop('t have to be nonnegative and finite numeric')
    }


    if (is.null(x)&is.null(tx)) {
      if (is.null(xmin)) {
        xmin = min(covariate)
      }
      else if (!is.numeric(xmin) |
               length(xmin) > 1 |
               max(xmin) == Inf) {
        stop('xmin have to be finite number')
      }
      if (is.null(xmax)) {
        xmax = max(covariate)
      }
      else if (!is.numeric(xmax) |
               length(xmax) > 1 |
               max(xmin) > min(xmax) |
               max(xmax) == Inf) {
        stop('xmax have to be finite number bigger than xmin, if indicated')
      }
      if (!is.numeric(x.length) |
          length(x.length) > 1 |
          min(x.length) < 1 |
          max(x.length) == Inf) {
        stop('x.length have to be possitive integer')
      }
      if (is.numeric(x.length) & length(x.length) == 1) {
        if (abs(x.length - round(x.length)) > 0) {
          stop('x.length have to be possitive integer')
        }
      }
      x <- seq(xmin, xmax, length.out = x.length)
    }
    x <- na.omit(x)
    if (!(is.numeric(x) &
          sum(  x == Inf) == 0)) {
      stop('x have to be finite numeric')
    }

##
    if (!is.null(tx)) {
      if (dim(tx)[2]!=2|(!is.numeric(tx[,1])|!is.numeric(tx[,2]))){
        stop('tx have to be matrix or data frame with 2 numeric columns')

      }
      x<-tx[,2]
      t<-tx[,1]

      t <- na.omit(t)
      if (!(is.numeric(t) &
            sum(t < 0 |
                t == Inf) == 0)) {
        stop('First column of tx have to be nonnegative and finite numeric')
      }
      x <- na.omit(x)
      if (!(is.numeric(x) &
            sum(  x == Inf) == 0)) {
        stop('Second column of tx have to be finite numeric')
      }
      if(! length(t)==length(x)){stop('Both columns have to have same length after omitting Nas.')}
    }
##



    if (h.method == "crossval") {
      h.FUN = crossval_haz_cond
    }
    if (h.method == "maxlike") {
      h.FUN = maxlike_haz_cond
    }

    h_original <- h
    if (is.null(h)) {
      a <- c((max(times) - min(times)) / length(times),(max(covariate) - min(covariate)) / length(covariate))
      b <- c(max(times),max(covariate))

      optfce <- function(h) {
        optfce <-
          h.FUN(
            h,
            Y = times,
            X= covariate,
            d = delta,
            kernel = kernel,
            type = type,type.w=type.w
          ) * (ifelse(h.method == "crossval", -1, 1))
        return(optfce)
      }


      if (optim.method == "ga") {
        gaControl("eps" = tol)
        GA.result <-
          ga(
            type = "real-valued",
            fitness = function(x)
               optfce(x),
            lower = a,
            upper = b,
            parallel = parallel,
            monitor = FALSE,
            run = run,
            ...
          )
        h <- as.numeric(summary(GA.result)$solution)
        CVML <-
          summary(GA.result)$fitness * (ifelse(h.method == "crossval", -1, 1))
      }


        hazard_est <-
          hazard_cond(
            t = t,
            x=x,
            Y = times,
            X=covariate,
            d = delta,
            kernel = kernel,
            h = h,
            type = type,
            type.w=type.w
          )



    } else if (is.numeric(h) & length(h) == 2 & min(h) > 0 &
               max(h) < Inf) {

      hazard_est <-
        hazard_cond(
          t = t,
          x=x,
          Y = times,
          X=covariate,
          d = delta,
          kernel = kernel,
          h = h,
          type = type,
          type.w=type.w
        )



        CVML <-
          h.FUN(
            h=h,
            Y = times,
            X= covariate,
            d = delta,
            kernel = kernel,
            type = type,type.w=type.w
          )



    } else {
      stop('h have to be positive finite two-dimensional vector')
    }



    if (is.null(h_original)) {
      details <-
        paste0(
          "h is optimal, estimated by statistical method ",
          h.method,
          " and numerical method ",
          optim.method,
          "."
        )
    }
    if (!is.null(h_original)) {
      details <- paste0("h is user choice. There is ",
                        h.method,
                        " computed in CVML.")
    }
    details <-
      paste0(
        details,
        " hazard is ",
        type,
        " kernel estimate of hazard function. Used kernel is ",
        kernel,
        "."
      )

    if (!exists("GA.result")) {GA.result=NULL}


    if (is.null(tx)) {
      hazlist <-
        list(
          time.points = t,
          covariate.points = x,
          hazard = hazard_est,
          h = h,
          CVML = CVML,
          details = details,
          GA.result = GA.result
        )
    } else {
      hazlist <-
        list(
          time.points = t,
          covariate.points = x,
          hazard = data.frame(t=t,x=x,hazard=diag(hazard_est)),
          h = h,
          CVML = CVML,
          details = details,
          GA.result = GA.result
        )

    }


    class(hazlist) <- "khazardcond"
    return(hazlist)
  }
