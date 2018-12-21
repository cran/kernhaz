trapezoid<-function(f,a,b,mx,...){
  ##inputs:f  ... function of two variables, type function
  #        a  ... lower border of interval for x
  #        b  ... upper border of interval for x
  #        mx ... number of intervals in direction x
  # output: aproximation of integral of function f on interval 
  #          [a,b] by trapeyoid method
  ##
  hx<-(b-a)/mx
  x<-seq(from=a,to=b,by=hx)
  nx<-length(x)
  
  H<-f(x,...)
  M<-rep(1,nx)
  M[c(1,nx)]<-1/2

  aprox_integral<-H%*%M*hx 
  return(aprox_integral)
}