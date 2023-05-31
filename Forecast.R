
ES = function(w,X,alpha){
  return(-mean(sort(X%*%c(w,1-sum(w)))[1:ceiling(alpha*length(X[,1]))]))
}

Forecast = function(steps,Assets,W){
  d = length(Assets[1,])-1;  S = numeric(d);  k = 0;  P = c()
  for(i in 1:length(Assets[,1])){
    S = S + Assets[i,-1];  k = k + 1
    if(k >= steps){
      P = c(P,S%*%c(W$par,1-sum(W$par)));  k = 0;  S = numeric(d)
    }
  }
  return(P)
}


ForecastCombine = function(steps,Assets,which,W){
  P = list()
  for(i in which){ P = append(P,list(Forecast(steps,Assets[[i]],W))) }
  ES = c();  ES.m = 0;  n.m = 0;  M = c();  M.m = 0;  n = 0
  for(i in 1:length(P)){ 
    ES = c(ES,-mean(sort(P[[i]])[1:floor(0.05*length(P[[i]]))]))
    ES.m = ES.m + ES[i]*floor(0.05*length(P[[i]]));  n.m = n.m + floor(0.05*length(P[[i]]))
    M = c(M,mean(P[[i]]));  M.m = M.m + M[i]*length(P[[i]]);  n = n + length(P[[i]])
  }
  ES = list(ES,ES.m/n.m);  M = list(M,M.m/n)
  return(list(P,ES,M))
}

NPCError = function(n,steps,N,NPC){
  Ps = c();  ESs = c()
  for(i in 1:N){
    time = Sys.time()
    W = NPCPortfolio(n,steps,NPC)
    Ps = rbind(Ps,W$par);  ESs = c(ESs,W$value)
    print(paste0(i," out of ",N,". It took ", round(Sys.time() - time,3), " seconds."))
  }
  P.m = rep(1,N)%*%Ps/N;  ES.m = mean(ESs)
  P.v = sqrt(diag(var(Ps)));  ES.v = sqrt(var(ESs))
  return(list(list(Ps,P.m,P.v),list(ESs,ES.m,ES.v)))
}

ESTest = function(Returns,which,W,plot=FALSE){
  X = c()
  for(i in which){ X = c(X,Returns[[i]]) }
  n = length(X);  ES = W$value
  real = c()
  for(i in 1:1000){ Y = sort(sample(X,n,TRUE));  real = c(real,-mean(Y[1:floor(0.05*n)])) }
  I = sort(real)[c(25,975)];  test = FALSE
  if((ES > I[1]) & (ES < I[2])){ test = TRUE }
  if(plot){
    hist(real,40,main="Histogram of Bootstrapped 5% Realized Shortfall",xlab="Shortfall")
    lines(c(ES,ES),c(0,10000),lty=2,col="red")
    lines(c(I[1],I[1]),c(0,10000),lty=2,col="blue")
    lines(c(I[2],I[2]),c(0,10000),lty=2,col="blue")
  }
  return(c(I,test))
}

PortId = function(W,E){
  df = length(W$par)
  stat = sum(((E[[1]][[2]]-W$par)/E[[1]][[3]])^2)
  p = 1-pchisq(stat,df)
  cstat = sum((solve(t(chol(var(E[[1]][[1]]))))%*%t(E[[1]][[2]]-W$par))^2)
  cp = 1-pchisq(cstat,df)
  return(c(df,stat,p,cstat,cp))
}

PortDist = function(E){
  d = length(E[[1]][[1]][1,])+1;  n = length(E[[1]][[1]][,1]);  ta = matrix(0,4,d+1)
  D = cbind(E[[1]][[1]],1-E[[1]][[1]]%*%rep(1,d-1))
  for(i in 1:d){
    jb = jarque.bera.test(D[,i])
    ta[,i] = c(mean(D[,i]),sqrt(var(D[,i])),jb$statistic,jb$p.value)
  }
  jb = jarque.bera.test(E[[2]][[1]])
  ta[,d+1] = c(mean(E[[2]][[1]]),sqrt(var(E[[2]][[1]])),jb$statistic,jb$p.value)
  return(round(ta,3))
}


