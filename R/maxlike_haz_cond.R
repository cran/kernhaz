#'@importFrom parallel detectCores
maxlike_haz_cond <- function(h,Y,X,d,kernel = "epanechnikov",type="interior",type.w="nw",parallel=FALSE) {
  ## Log of likelihood function
  # inputs: h       ...smoothing parameters, matrix with two columns (first column for time, second for covariate)
  #         Y       ...observed times, vector
  #         d       ...censoring indicator, vector, 0..censor, 1..event
  #         kernel  ...one of: 'epan','gauss','quart','rect','trian'
  #                    'epan'... default
  #         type    ...type of estimate,
  #                    possible: 'ex' ... exterior
  #                              'in' ... interior, default
  #         parallel...way of calculating
  #                    FALSE ... common way
  #                    TRUE  ... use more cores of computer
  # output: vector of values of loglikelihood for smoothing parameters in h
  ##

  ##dodelat popisky a kontroly vstupu, pozor na vektor h osetrit



  if(is.vector(h)){h<-matrix(h,ncol=2)}



  #a<-min(Y)
  a=0
  #b<-max(Y)

  n<-length(Y)

  nh<-dim(h)[1]




  ################
  ################   parallel==FALSE
  ################


  if(parallel==FALSE){


    ML<-numeric(nh)



    for (i in 1:nh) {


      S<-0
      for (j in 1:n) {

        S<-S+log(hazard_cond(t=Y[j],x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w)^d[j]*exp(-trapezoid(hazard_cond,a,Y[j],500,x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w)))


      }

      ML[i]<-S

    }

  }


  ################
  ################   parallel==TRUE
  ################

  if(parallel==TRUE){
    i<-1:nh
    no_cores <- detectCores() - 1
    registerDoParallel(no_cores)


    ML<-foreach(i=i, .combine='c',.packages="foreach",.export=c("hazard_cond","kern","trapezoid","wmatrix")) %dopar%{
      S<-0
      for (j in 1:n) {
        S<-S+log(hazard_cond(t=Y[j],x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w)^d[j]*exp(-trapezoid(hazard_cond,a,Y[j],500,x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w)))

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
