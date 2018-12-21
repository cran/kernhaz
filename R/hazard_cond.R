hazard_cond <- function(t,x,Y,X,d,h,kernel = "epanechnikov",type='interior',type.w='nw') {
  ## KERNEL ESTIMATION OF HAZARD FUNCTION
  # inputs: t     ...time points in which we want calculate hazard, vector
  #         x     ...vector of covariates we want to calculate hazard
  #         Y     ...observed times, vector
  #         X     ...vector of observed covariates
  #         d     ...censoring indicator, vector, 0..censor, 1..event
  #         h     ...vector of smoothing parameter h=c(ht,hx)
  #           ht    ...smoothing parameter, number
  #           hx    ...parametr for Nadaraya-Watson waights
  #         kernel...one of: 'epan','gauss','quart','rect','trian'
  #                  'epan'... default
  #         type  ...type of estimate, 
  #                  possible: 'ex' ... exterior
  #                            'in' ... interior, default
  # output: matrix of hazard values in times t and covariates x
  #         j-th row of matrix is vector of hazard values in times t for j-th covariate in x
  #         i-th column of matrix is vector of hazard values in i-th time in t for vector of covariate x
  #         j-th row and i-th column .. i-th time, j-th covariate
  ##
  

 
  ###
  hx<-h[2]
  ht<-h[1]
  
  K<-kern(kernel)$K
  W<-kern(kernel)$W
  n<-length(Y)
  p<-length(t)
  ####################
  ####################
  ###################
  
  
  if(type=='exterior'){
    Ymat<-(matrix(rep(t,n),nrow=n,byrow=T)-matrix(rep(Y,p),ncol=p))/ht
   
    WM<-wmatrix(x=x,X=X,hx=hx,kernel = kernel,type=type.w)
  
    numerator<-(1/ht)* (WM%*%(K(Ymat) *d))
    denominator<-WM%*%W(-Ymat)
    denominator[denominator<=0.005]=0.005
  
    condhaz_t<-numerator/denominator
  }
  
  if(type=='interior'){
    Yord<-Y[order(Y)]
    Xord<-X[order(Y)]
    dord<-d[order(Y)]
    Ymat<-(matrix(rep(t,n),nrow=n,byrow=T)-matrix(rep(Yord,p),ncol=p))/ht
    
    WM<-wmatrix(x=x,X=Xord,hx=hx,kernel = kernel,type=type.w)

    partsums_matrix<-WM%*% upper.tri(matrix(1,n,n))
    partsums_matrix[partsums_matrix==1]<-1-0.00000001
    
    condhaz_t<-(1/ht)*((WM/(1-partsums_matrix))%*%(K(Ymat)*dord))
    
    
  }
  
  
  #poresit NA
  condhaz<-condhaz_t
  
  # one row = all times, one covariate
  return(condhaz) 
}