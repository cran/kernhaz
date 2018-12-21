#'@importFrom parallel detectCores
maxlike_haz_ln <- function(h,Y,d,kernel = "epanechnikov",type="interior",parallel=FALSE) {
  ## MAXLIKE
  # inputs: h       ...smoothing parameters, vector
  #         Y       ...observed times, vector
  #         d       ...censoring indicator, vector, 0..censor, 1..event
  #         kernel  ...one of: 'epanechnikov','gaussian','rectangular','quartic'
  #                    'epanechnikov'... default
  #         type    ...type of estimate,
  #                    possible: 'exterior' ... exterior
  #                              'interior' ... interior, default
  #         parallel...way of calculating
  #                    FALSE ... common way
  #                    TRUE  ... use more cores of computer, default
  # output:  ML     ...vector of values of log likelihood for smoothing parameters in h
  #
  ##






  n<-length(Y)

  nh<-length(h)

  ML<-numeric(nh)
  exp_ML<-numeric(nh)

  ################
  ################   parallel==FALSE
  ################


  if(parallel==FALSE){


  for (i in 1:nh) {

    S<-0

    for (j in 1:n) {
      S<-S+log(hazard(Y[j],Y[-j],d[-j],h[i],kernel,type)^d[j]*exp(-trapezoid(hazard,0,Y[j],500,Y[-j],d[-j],h[i],kernel,type)))
    }
    ML[i]<-S


  }

  }



  ################
  ################   parallel==TRUE
  ################

  if(parallel==TRUE){
    no_cores <- detectCores() - 1
    registerDoParallel(no_cores)


    ML<-foreach(h=h, .combine='c',.packages="foreach",.export=c("hazard","kern","trapezoid")) %dopar%{
      S<-0
      for (j in 1:n) {
        S<-S+log(hazard(Y[j],Y[-j],d[-j],h,kernel,type)^d[j]*exp(-trapezoid(hazard,0,Y[j],500,Y[-j],d[-j],h,kernel,type)))

      }
      S

    }
    stopImplicitCluster()




  }


  ################
  ################
  ################
  ML[ML==-Inf]<--10^300

  return(ML)

}
