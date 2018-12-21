hazard <- function(t,Y,d,h,kernel = "epanechnikov", type='interior') {
  ## KERNEL ESTIMATION OF HAZARD FUNCTION
  # inputs: t     ...time points in which we want calculate hazard, vector
  #         Y     ...observed times, vector
  #         d     ...censoring indicator, vector, 0..censor, 1..event
  #         h     ...smoothing parameter, number
  #         kernel...one of: 'epan','gauss','quart','rect','trian'
  #                  'epan'... default
  #         type  ...type of estimate, 
  #                  possible: 'ex' ... exterior
  #                            'in' ... interior, default
  # output: vector of hazard values in times t 
  ##
  

  
  K<-kern(kernel)$K
  W<-kern(kernel)$W
  n<-length(Y)
  p<-length(t)
  ##############################
  ##############################
  ##############################
  
  if(type=='exterior'){
  Ymat<-(matrix(rep(t,n),ncol=n)-matrix(rep(Y,p),nrow=p,byrow=T))/h
  numerator<-(1/h)*K(Ymat)%*%d
  denominator<-W(-Ymat)%*%matrix(1,n,1)
  denominator[denominator<=0.005]=0.005
  
  hazard<-as.vector(numerator/denominator)
  }
  ###
  if(type=='interior'){
    Yord<-Y[order(Y)]
    dord<-d[order(Y)]
    Ymat<-(matrix(rep(t,n),ncol=n)-matrix(rep(Yord,p),nrow=p,byrow=T))/h
    dord_divided<-dord/(n-(1:n)+1)
    
    hazard<-as.vector((1/h)*K(Ymat)%*%dord_divided)
  }
  ###
  

  return(hazard) 
}