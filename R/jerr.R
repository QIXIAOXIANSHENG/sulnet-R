jerr=function(n,maxit,pmax,family){
  if(n==0) list(n=0,fatal=FALSE,msg="")
  else {
    errlist=switch(family,
              "gaussian"=jerr.elnet(n,maxit,pmax),
              "cox"=jerr.coxnet(n,maxit,pmax)
      )
    names(errlist)=c("n","fatal","msg")
    errlist$msg=paste("from sulnet Fortan code (error code ",n, "); ",errlist$msg,sep="")
    errlist
  }
}
                  
                  
     
