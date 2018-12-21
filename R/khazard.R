#'@description Kernel estimate of (unconditional) hazard function for right-censored data. Options include two methods for bandwidth selection.
#'@title Kernel estimate of hazard function for right-censored data
#'@param times vector of observed times
#'@param delta vector of censoring indicator. 0 - censored, 1 - uncensored (dead)
#'@param h bandwidth (scalar or vector). If missing, h is found using some bandwidth selection method.
#'@param t vector of time points at which estimate is evaluated
#'@param t.length number of grid points
#'@param tmin,tmax minimum/maximum values for grid
#'@param kernel kernel function, possible values are: "epanechnikov" (default), "gaussian", "rectangular", "quartic".
#'@param type Type of kernel estimate. Possible types are:  "exterior", "interior" (default).
#'@param parallel  allows parallel computation. Default is FALSE.
#'@param value If h parameter is vector, this option controls output values. If "CVML" (default), the crossvalidation or log-likelihood values only are calculated. If "hazard", the hazard functions only are calculated. If "both" the crossvalidation or log-likelihood values and hazard function are calculated.
#'@param h.method method for bandwidth selection. Possible methods are: "crossval" (default), "maxlike".
#'@param optim.method method for numerical optimization of the crossvalidation or log-likelihood function. Possible methods are: "optimize" (default), "ga".
#'@param tol the desired accuracy of optimization algorithm
#'@param run the number of consecutive generations without any improvement in the best fitness value before the GA is stopped.
#'@param ... additional arguments of GA algorithm
#'@return Returns an object of class 'khazard' which is a list with fields
#'\item{time.points}{vector of time points at which estimate is evaluated}
#'\item{hazard}{data frame of time points, hazard function values and bandwidth}
#'\item{h}{bandwidth}
#'\item{CVML}{value of crossvalidation or log-likelihood at h}
#'\item{details}{description of used methods}
#'\item{GA.result}{output of ga, object of class ga-class}
#'@seealso \code{\link{plot.khazard}}, \code{\link[GA]{ga}}, \code{\link{optimize}}
#'@examples
#'library(survival)
#'fit<-khazard(times = lung$time,delta = lung$status-1)
#'@details External type of kernel estimator is defined as the ratio of kernel estimator of the subdensity of the uncensored observations to the survival function of the observable time. Internal type of kernel estimator is based on a convolution of the kernel function with a nonparametric estimator of the cumulative hazard function (Nelson-Aalen estimator).
#'@references Selingerova, I., Dolezelova, H., Horova, I., Katina, S., and Zelinka, J. (2016). Survival of Patients with Primary Brain Tumors: Comparison of Two Statistical Approaches. PloS one, 11(2), e0148733.
#'@import foreach
#'@import doParallel
#'@import GA
#'@importFrom stats na.omit optimize
#'@importFrom parallel detectCores
#'@export
khazard <-
  function(times,
           delta,
           h = NULL,
           t = NULL,
           t.length = 100,
           tmin = NULL,
           tmax = NULL,
           kernel = "epanechnikov",
           type = "interior",
           parallel = FALSE,
           value = "CVML",
           h.method = "crossval",
           optim.method = "optimize",
           tol = ifelse(h.method == "crossval", 10 ^ (-6), 1),
           run = 2,
           ...) {
    ### verification of inputs + omission of NAs

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

    possible_value <- c('cvml', 'hazard', 'both')
    value <- possible_value[pmatch(tolower(value), possible_value)]
    if (value == 'cvml') {
      value = 'CVML'
    }
    if (is.na(value)) {
      stop("value have to be one of: 'CVML','hazard','both'")
    }

    possible_h.method <- c('crossval', 'maxlike')
    h.method <-
      possible_h.method[pmatch(tolower(h.method), possible_h.method)]
    if (is.na(h.method)) {
      stop("h.method have to be one of: 'crossval','maxlike'")
    }

    possible_optim.method <- c('optimize', 'ga')
    optim.method <-
      possible_optim.method[pmatch(tolower(optim.method), possible_optim.method)]
    if (is.na(optim.method)) {
      stop("optim.method have to be one of: 'optimize','ga'")
    }

    if (missing(times) | missing(delta))
      stop("Parameter times or delta is missing")
    if (length(times) != length(delta)) {
      stop('times and delta have to be of same length')
    }

    notNA <- !(is.na(times) | is.na(delta))
    times <- times[notNA]
    delta <- delta[notNA]
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


    if (is.null(t)) {
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


    if (h.method == "crossval") {
      h.FUN = crossval_haz
    }
    if (h.method == "maxlike") {
      h.FUN = maxlike_haz_ln
    }

    h_original <- h
    if (is.null(h)) {
      a <- (max(times) - min(times)) / length(times)
      b <- max(times)
      optfce <- function(h) {
        optfce <-
          h.FUN(
            h,
            Y = times,
            d = delta,
            kernel = kernel,
            type = type
          ) * (ifelse(h.method == "crossval", -1, 1))
        return(optfce)
      }

      if (optim.method == "optimize") {
        op.result <- optimize(optfce, c(a, b), maximum = TRUE, tol = tol)
        h <- as.numeric(op.result$maximum)
        CVML <- op.result$objective * (ifelse(h.method == "crossval", -1, 1))
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
          hazard(
            t = t,
            Y = times,
            d = delta,
            kernel = kernel,
            h = h,
            type = type
          )
        hazard_df <-
          data.frame(t = t,
                     hazard = hazard_est,
                     h = h)


    } else if (is.numeric(h) & length(h) == 1 & min(h) > 0 &
               max(h) < Inf) {

        hazard_est <-
          hazard(
            t = t,
            Y = times,
            d = delta,
            kernel = kernel,
            h = h,
            type = type
          )
        hazard_df <-
          data.frame(t = t,
                     hazard = hazard_est,
                     h = h)


        CVML <-
          h.FUN(
            h = h,
            Y = times,
            d = delta,
            kernel = kernel,
            type = type
          )


    }  else if (is.numeric(h) & length(h) > 1 & min(h) > 0 &
                max(h) < Inf) {
      if (value == 'hazard' | value == "both") {
        if (parallel == TRUE) {
          no_cores <- detectCores() - 1
          registerDoParallel(no_cores)
          hazard_est <-
            foreach(
              h = h,
              .combine = 'c',
              .packages = "foreach",
              .export = c("hazard", "kern")
            ) %dopar% {
              hazard(
                t = t,
                Y = times,
                d = delta,
                kernel = kernel,
                h = h,
                type = type
              )


            }
          stopImplicitCluster()
        }
        else {
          hazard_est <- rep(NA, length(h) * length(t))
          for (i in 1:length(h)) {
            hazard_est[((i - 1) * length(t) + 1):(i * length(t))] <-
              hazard(
                t = t,
                Y = times,
                d = delta,
                kernel = kernel,
                h = h[i],
                type = type
              )
          }
        }
        hazard_df <- data.frame(
          t = rep(t, times = length(h)),
          hazard = hazard_est,
          h = rep(h, each = length(t))
        )
      }
      if (value == 'CVML' |
          value == "both") {
        CVML <-
          h.FUN(
            h = h,
            Y = times,
            d = delta,
            kernel = kernel,
            type = type,
            parallel = parallel
          )
      }


    }  else if (!is.numeric(h) | min(h) <= 0 | max(h) == Inf) {
      stop('h have to be positive finite numeric')
    }


    ####
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
    if (length(h_original)>1 & value=="hazard") {
      details <- paste0("h is user choice.")
    }

    if (value == 'hazard' |
        value == "both" |
        (is.null(h_original) | length(h_original)==1) ) {
      details <-
        paste0(
          details,
          " hazard is ",
          type,
          " kernel estimate of hazard function. Used kernel is ",
          kernel,
          "."
        )
    }
    ####
    if (!exists("GA.result")) {GA.result=NULL}

    if (value == 'both'|length(h)==1) {
      hazlist <-
        list(
          time.points = t,
          hazard = hazard_df,
          h = h,
          CVML = CVML,
          details = details,
          GA.result = GA.result
        )
    }
    if (value == 'CVML'&length(h)>1) {
      hazlist <- list(time.points = NULL,
                      h = h,
                      hazard = NULL,
                      CVML = CVML,
                      details = details,
                      GA.result = GA.result)
    }
    if (value == 'hazard'&length(h)>1) {
      hazlist <- list(time.points = t,
                      h = h,
                      hazard = hazard_df,
                      CVML=NULL,
                      details = details,
                      GA.result = GA.result)
    }



    class(hazlist) <- "khazard"
    return(hazlist)
  }
