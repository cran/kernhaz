#'@title Plot of kernel hazard estimate from an object of class khazard
#'@description Plot of kernel hazard estimate from an object of class khazard
#'
#'@param x Object of class khazard
#'@param h bandwidth for which hazard function estimate will be plot if x$h is vector
#'@param ylim Limits for the y axis.
#'@param type type argument for plot.
#'@param xlab Label for the x axis.
#'@param ylab Label for the y axis.
#'@param main Title of plot.
#'@param ... Additional arguments.
#'@seealso \code{\link{khazard}}
#'@examples library(survival)
#'fit<-khazard(times = lung$time,delta = lung$status-1)
#'plot(fit)
#'
#'fit<-khazard(times = lung$time,delta = lung$status-1,h=c(100,150,200,250), value="both")
#'plot(fit,h=200)
#'@method plot khazard
#'@export
plot.khazard<-function (x,h=NULL, ylim, type, xlab, ylab, main, ...)
{
  if (is.null(x$hazard)){
    stop("values of hazard rate are unknown")
  }


  if (!is.null(h)){
    if (min(h%in%x$h) & length(h)==1){
      xx<-x$hazard$t[x$hazard$h==h]
      y <- x$hazard$hazard[x$hazard$h==h]
    }
    else {
      stop("Parameter h is longer than one or values of hazard rate are unknown for this parameter h.")
    }

  }

  if (is.null(h)& length(x$h)==1){
    xx<-x$hazard$t
    y <- x$hazard$hazard
    h<-signif(x$h,4)

  }
  if (is.null(h)& length(x$h)>1){
    stop("parameter h is necessary")
  }


  if (missing(ylim))
    ylim <- c(0, max(y))
  if (missing(type))
    type <- "l"
  if (missing(xlab))
    xlab <- "Time"
  if (missing(ylab))
    ylab <- "Hazard Rate"
  if (missing(main))
    main <- paste("Kernel estimation of hazard function, h =",signif(h, digits = 3))
  plot(xx, y, type, ylim = ylim, xlab = xlab, ylab = ylab,main=main,
       ...)
  return(invisible())
}
