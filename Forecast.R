
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


GC = GCCombination(Assets,1:2)
W.GC = GCPortfolio(100000,10,GC)
P.GC = ForecastCombine(10,Assets,3:12,W.GC)


NPC = NPCCombination(Assets,1:2,c(),1:12)
W.NPC = NPCPortfolio(100000,10,NPC)
P.NPC = ForecastCombine(10,Assets,3:12,W.NPC)


NPCAR = NPCCombination(Assets,1:2,c(12),1:12)
W.NPCAR = NPCPortfolio(100000,10,NPCAR)
P.NPCAR = ForecastCombine(10,Assets,3:12,W.NPCAR)

P.NPC[[2]]
P.NPCAR[[2]]

round(1000*P.GC[[2]][[2]],2)
round(1000*P.NPC[[2]][[2]],2)
round(1000*P.NPCAR[[2]][[2]],2)


sum(W.GC$par)
sum(W.NPC$par)
sum(W.NPCAR$par)
GC[[3]]


pbinom(7,10,0.33)




