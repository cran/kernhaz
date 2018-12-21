#'@importFrom parallel detectCores
crossval_haz <- function(h,Y,d,kernel = "epanechnikov",type="interior",parallel=FALSE) {
  ## CROSSVALIDATION
  # inputs: h       ...smoothing parameters, vector
  #         Y       ...observed times, vector
  #         d       ...censoring indicator, vector, 0..censor, 1..event
  #         kernel  ...one of: 'epanechnikov','gaussian','rectangular','quartic'
  #                    'epan'... default
  #         type    ...type of estimate,
  #                    possible: 'exterior' ... exterior
  #                              'interior' ... interior, default
  #         parallel...way of calculating
  #                    FALSE ... common way
  #                    TRUE  ... use more cores of computer
  # output: vector of values of crossvalidation for smoothing parameters in h
  ##





  a=0


  n<-length(Y)

  nh<-length(h)

  type_s<-ifelse(type=='interior',"empirical","kernsmooth")


  haz_sqr<-function(t,Y,d,h,kernel,type){
    hazard(t,Y,d,h,kernel,type)^2
  }

  ################
  ################   parallel==FALSE
  ################


  if(parallel==FALSE){


  CV<-numeric(nh)



    for (i in 1:nh) {

      S<-0
      for (j in 1:n) {
        S<-S+hazard(Y[j],Y[-j],d[-j],h[i],kernel,type)*d[j]/surv(Y[j],Y,d,h[i],kernel,type_s)
          #L_n(Y[j],Y)

      }
      b=max(Y)+h[i]
      CV[i]<-trapezoid(haz_sqr,a,b,1000,Y,d,h[i],kernel,type)-2/n*S

    }

  }


  ################
  ################   parallel==TRUE
  ################

  if(parallel==TRUE){

    no_cores <- detectCores() - 1
    registerDoParallel(no_cores)



    CV<-foreach(h=h, .combine='c',.packages="foreach",.export=c("hazard","kern","trapezoid","surv")) %dopar%{
      S<-0
      for (j in 1:n) {
        S<-S+hazard(Y[j],Y[-j],d[-j],h,kernel,type)*d[j]/surv(Y[j],Y,d,h,kernel,type_s)

      }
      b=max(Y)+h
      trapezoid(haz_sqr,a,b,1000,Y,d,h,kernel,type)-2/n*S
    }
    stopImplicitCluster()



  }


  ################
  ################
  ################


  return(CV)

}
