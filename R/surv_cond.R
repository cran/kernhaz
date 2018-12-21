surv_cond<-function(t,x,Y,X,d,h,kernel = "epan",type.w="nw") {
  ## KERNEL ESTIMATION OF SURVIVAL FUNCTION
  # inputs: t     ...time points in which we want calculate hazard, vector
  #         x     ...vector of covariates we want to calculate hazard
  #         Y     ...observed times, vector
  #         X     ...vector of observed covariates
  #         d     ...censoring indicator, vector, 0..censor, 1..event
  #         h     ...vector of smoothing parameter h=c(ht,hx)
  #           ht    ...smoothing parameter, number
  #           hx    ...parametr for Nadaraya-Watson waights
  #         kernel...one of: 'epan','gauss','quart','rect','trian'
  # output: matrix of hazard values in times t and covariates x
  #         j-th row of matrix is vector of hazard values in times t for j-th covariate in x
  #         i-th column of matrix is vector of hazard values in i-th time in t for vector of covariate x
  #         j-th row and i-th column .. i-th time, j-th covariate
  ##
  
  #predelano ht, hx jako vektor, zkontrolovat vstupy...
  
  ### verification of inputs + omission of NAs
  
  

  
  ###
  hx<-h[2]
  ht<-h[1]
  
  W<-kern(kernel)$W
  n<-length(Y)
  p<-length(t)
  
  Ymat<-(matrix(rep(t,n),ncol=n)-matrix(rep(Y,p),nrow=p,byrow=T))/ht
  
  WM<-wmatrix(x=x,X=X,hx=hx,kernel = kernel,type=type.w)
  
  denominator<-(WM*W(-Ymat)) %*%matrix(1,n,1)
  denominator[denominator<=0.005]=0.005
  
  
  
  cond_surv<-denominator
  
  # one row = all times, one covariate
  return(cond_surv) 
}