maximinSLHD <-
function(t,m,k,power=15,nstarts=1,itermax=100,total_iter=1000000){
  
  n=m*t
  
  t00<-Sys.time()
  aaa<-.C("maximinSLHD",as.integer(m),as.integer(k),as.integer(t),as.integer(power),as.integer(nstarts),
          as.integer(itermax),as.integer(total_iter),design=integer(m*t*k),measure=double(1), 
          temp0=double(1),PACKAGE="SLHD")
  t01<-Sys.time()
  
  time_rec=t01-t00
  
  dd=matrix(aaa$design,ncol=k,nrow=n,byrow=TRUE)
  
  Dslice<-list()
  for(i in 1:t){
    Dslice[[i]]<-dd[((i-1)*m+1):(i*m),]
  }
  
  scaled_deisgn<-(apply(dd,2,rank)-0.5)/n
  
  SDslice<-list()
  for(i in 1:t){
    SDslice[[i]]<-scaled_deisgn[((i-1)*m+1):(i*m),]
  }
  
  val<-list(D=dd,DS=Dslice,standD=scaled_deisgn,standDS=SDslice,temp0=aaa$temp0,measure=aaa$measure,time_rec=time_rec)
  
  return(val)
}
