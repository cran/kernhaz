surv<-function(t,Y,d,h=NA,kernel='epanechnikov',type='kernsmooth'){
  # KERNEL ESTIMATION OF SURVIVAL FUNCTION in times t
  # inputs: t     ...time points in which we want calculate survival function, vector
  #         Y     ...observed times, vector
  #         d     ...censoring indicator, vector, 0..censor, 1..event
  #         h     ...smoothing parameter, number, needed if type='ex'
  #         kernel...one of: 'epan','gauss','quart','rect','trian'
  #                  'epan'... default
  #         type  ...type of estimate, 
  #                  possible: 'empirical'  ... exterior
  #                            'kernsmooth' ... interior, default
  # output: vector of survival function values in times t 
  ##
  
   possible_type <- c('empirical','kernsmooth')
   type2 <- possible_type[pmatch(tolower(type), possible_type)]
   if(is.na(type2)){stop("function surv: type have to be one of: 'empirical','kernsmooth'")}
   
  ## 
  
  
  ###
  n<-length(Y)
  p<-length(t)
  ####
  ####
  ####
  
  
  if(type2=='empirical'){
    I<-matrix(rep(Y,p),nrow=p,byrow=T)<=matrix(rep(t,n),nrow=p)
    surv_t<-1-I%*%matrix(1,nrow=n,ncol=1)/(n+1)
  }
  
  
  if(type2=='kernsmooth'){
    
    if(!(is.numeric(h)  & sum(h>0)==length(h) & length(h)==1)){stop('h have to be positive numeric of length one')}
  
    K<-kern(kernel)$K
    W<-kern(kernel)$W
    
  Ymat<-(matrix(rep(t,n),ncol=n)-matrix(rep(Y,p),nrow=p,byrow=T))/h
  denominator_divided_by_n<-W(-Ymat)%*%matrix(1,n,1)/n
  denominator_divided_by_n[denominator_divided_by_n<=0.005]=0.005
  surv_t<-denominator_divided_by_n
  }
  
  return(surv_t)
}