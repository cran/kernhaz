#'@importFrom parallel detectCores
crossval_haz_cond <- function(h,Y,X,d,kernel = "epanechnikov",type="interior",type.w="nw",parallel=FALSE) {
  ## CROSSVALIDATION
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
  # output: vector of values of crossvalidation for smoothing parameters in h
  ##



  if(is.vector(h)){h<-matrix(h,ncol=2)}


  #a<-min(Y)
  a=0
  #b<-max(Y)

  n<-length(Y)

  nh<-dim(h)[1]




  I1<-function(t,x,Y,X,d,h,kernel,type,type.w){
    ##t is vector of equidistant time point with step hh
    hh<-t[2]-t[1]
    haz<-hazard_cond(t=t,x=x,Y=Y,X=X,d=d,h=h,kernel=kernel,type=type,type.w=type.w)
    n<-length(t)
    M<-matrix(0,nrow = n,ncol = n)
    M[upper.tri(M, diag = FALSE)]<-1
    diag(M)<-1/2
    M[1,]<-1/2
    M[1,1]<-0


    inthaz<-(haz%*%M)*hh

    I1<-haz^3*exp(-inthaz)
    return(I1)
  }

  ################
  ################   parallel==FALSE
  ################


  if(parallel==FALSE){


    CV<-numeric(nh)



    for (i in 1:nh) {
      b=max(Y)+h[i,1]
      I1pom<-trapezoid(I1,a,b,500,x=X,Y=Y,X=X,d=d,h=h[i,],kernel=kernel,type=type,type.w=type.w)
      S<-0
      for (j in 1:n) {

        S<-S+1/n*I1pom[j]-
          2/n*hazard_cond(t=Y[j],x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w)^2/surv_cond(t=Y[j],x=X[j],Y=Y,X=X,d=d,h=h[i,],kernel =kernel,type.w=type.w)*d[j]*exp(-trapezoid(hazard_cond,a,Y[j],500,x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w))

      }

      CV[i]<-S

    }

  }


  ################
  ################   parallel==TRUE
  ################


  if(parallel==TRUE){
    i<-1:nh

    no_cores <- detectCores() - 1
    registerDoParallel(no_cores)

    CV<-foreach(i=i, .combine='c',.packages="foreach",.export=c("hazard_cond","kern","trapezoid","surv_cond","wmatrix")) %dopar%{
      b=max(Y)+h[i,1]
      I1pom<-trapezoid(I1,a,b,500,x=X,Y=Y,X=X,d=d,h=h[i,],kernel=kernel,type=type,type.w=type.w)
      S<-0
      for (j in 1:n) {
        S<-S+1/n*I1pom[j]-
          2/n*hazard_cond(t=Y[j],x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w)^2/surv_cond(t=Y[j],x=X[j],Y=Y,X=X,d=d,h=h[i,],kernel =kernel,type.w=type.w)*d[j]*exp(-trapezoid(hazard_cond,a,Y[j],500,x=X[j],Y=Y[-j],X=X[-j],d=d[-j],h=h[i,],kernel=kernel,type=type,type.w=type.w))

      }
      S
    }
    stopImplicitCluster()



  }


  ################
  ################
  ################


  return(CV)

}
