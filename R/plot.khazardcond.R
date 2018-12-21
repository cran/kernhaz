#'@title Plot of kernel conditional hazard estimate from an object of class khazardcond
#'@description Plot of kernel conditional hazard estimate from an object of class khazardcond
#'
#'@param x Object of class khazardcond
#'@param type type of plot. Possible types are:  "persp" (default), "persp3d", "contour".
#'@param zlim Limits for the z axis.
#'@param xlab Label for the x axis.
#'@param ylab Label for the y axis.
#'@param zlab Label for the z axis.
#'@param main Title of plot.
#'@param ... Additional arguments.
#'@seealso \code{\link{khazardcond}}
#'@examples library(survival)
#'fit<-khazardcond(times = lung$time,delta = lung$status-1,covariate = lung$age,h=c(200,20))
#'plot(fit)
#'@import rgl
#'@method plot khazardcond
#'@export
#'@importFrom graphics persp contour

plot.khazardcond<-function (x,type="persp", zlim, xlab, ylab,zlab,main, ...)
{

  possible_type <- c('persp', 'persp3d','contour')
  type <- possible_type[pmatch(tolower(type), possible_type)]
  if (is.na(type)) {
    stop("type have to be one of: 'persp', 'persp3d','contour'")
  }

  if (is.data.frame(x$hazard)){
    stop("values of hazard rate have to be in matrix form corresponding grid of time and covariate")
  }

  xx<-x$time.points
  y<-x$covariate.points
  z<-t(x$hazard)
  h<-x$h



  if (missing(zlim))
    zlim <- c(0, max(z))

  if (missing(xlab))
    xlab <- "Time"
  if (missing(ylab))
    ylab <- "Covariate"
  if (missing(zlab))
    zlab <- "Hazard Rate"
  if (missing(main))
    main <- paste0("Kernel estimation of hazard function, h = (",signif(h[1], digits = 3),",",signif(h[2], digits = 3),").")
  if (type=='persp') {
    persp(xx, y, z, zlim = zlim, xlab = xlab, ylab = ylab,zlab=zlab,main=main,
          ...)
  }

  if (type=='persp3d') {
    persp3d(xx, y, z, zlim = zlim, xlab = xlab, ylab = ylab,zlab=zlab,main=main,
          ...)
  }

  if (type=='contour') {
    contour(xx, y, z, zlim = zlim, xlab = xlab, ylab = ylab,main=main,
          ...)
  }


  return(invisible())
}
