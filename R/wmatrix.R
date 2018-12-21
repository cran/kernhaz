wmatrix <- function(x,X,hx=1,kernel = "epan",type='nw',p=0) {
  ## WEIGHTS
  # inputs: x     ...vector of covariates we want to calculate
  #         X     ...vector of observed covariate
  #         hx    ...smoothing parameter, number
  #         kernel...one of: 'epan','gauss','quart','rect','trian'    ################ trian neznam coef problem u gm !!!!!!!!!!!!
  #                  default: 'epan'
  #         type  ...type of weigts
  #                  'nw'... Nadaraya-Watson, default
  #                  'll'... local-linear
  #                  'lp'... local-polynomial
  #                             p - order of the local polynomial
  #                                 p=0 - Nadaraya-Watson estimation, default
  #                                 p=1 - local linear estimation
  #                  'gm'... Gasser-Muller
  # output: martix of weights for vector x of covariates
  ##

  ### verification of inputs + omission of NAs
  X2<-na.omit(X)
  x2<-na.omit(x)
  if(!(is.numeric(X2) & sum(-Inf<X2 & X2<Inf )==length(X2) )){stop('X have to be finite numeric')}
  if(!(is.numeric(hx)  & sum(hx>0)==length(hx) & length(hx)==1)){stop('hx have to be positive numeric of length one')}
  ###

  K<-kern(kernel)$K
  n<-length(X2)
  m<-length(x2)
#################################################################################
  if(type=='nw' | (type=='lp' & p==0)){
    Xmat<-(matrix(rep(x2,n),ncol=n)-matrix(rep(X2,m),nrow=m,byrow=T))/hx  # pxn
    # (x1-X1)/hx     ...    (x1-Xn)/hx
    #
    #
    # (xp-X1)/hx     ...    (xp-Xn)/hx
    numerator<-K(Xmat)
    denominator<-K(Xmat)%*%matrix(1,n,1)  # p x 1

    M<-numerator/as.vector(denominator)    # p-th row = p-th covariate in x2

  }

###############################################################################
  if(type=='ll' | (type=='lp' & p>0)){
    if(type=='ll'){pp<-1}
    if(type=='lp'){pp<-p}

    M<-matrix(NA,ncol=m,nrow=n)
    powers<-0:pp

    for(i in 1:m){
      pomx<-x2[i]-X2
      w<-K(pomx/hx)
      Loopy<-lapply(powers,function(s){
        pomx^s
      })
      Xx<-matrix(unlist(Loopy),ncol=pp+1)
      Xxw<-Xx*matrix(w,ncol=pp+1,nrow=n)
      A<-solve(t(Xxw)%*%Xx)%*%t(Xxw)
      M[,i]<-A[1,]
    }
    M<-t(M)
  }
######################################################################
  if(type=='gm'){
    pol<-kern(kernel)$details$coef
    # scaling
    X2<-sort((X2-min(X2))/(max(X2)-min(X2)))
    x2<-(x2-min(x2))/(max(x2)-min(x2))

    pol_int_val<-function(pol,points){
      # value of integral of polynom pol in point
      n_pol<-length(pol)
      int_pol<-c(pol,0)/c(n_pol:1,1)
      values<-int_pol*rep(points,each=n_pol+1)^c(n_pol:0)
      indicator_ofsum<-rep(1:length(points),each=n_pol+1)
      values_final<-tapply(values, list(indicator_ofsum), sum)
      return(values_final)
    }

    M<-matrix(NA,ncol=m,nrow=n)
    X2<-c(0,X2)
    nn<-length(X2)
    for(i in 1:m){
      vecup<-X2[c(nn,2:nn)]
      vecdown<-c(0,X2[1:(nn-1)])
      matup<-cbind(vecup,x2[i]+hx)
      matdown<-cbind(vecdown,x2[i]-hx)

      up<-apply(matup,1,min)
      down<-apply(matdown,1,max)

      cons<-pol_int_val(-pol,(x2[i]-up[1])/hx)-pol_int_val(-pol,(x2[i]-down[1])/hx)
      cons<-as.numeric(cons)

      MM<-(up>down)*(1/cons)*
        (pol_int_val(-pol,(x2[i]-up)/hx)-pol_int_val(-pol,(x2[i]-down)/hx))

      M[,i]<-MM[-1]
    }

    M<-t(M)
  }
#################################################################################

  return(M)
}
