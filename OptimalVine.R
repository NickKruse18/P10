
OVFind = function(Assets,which,m=100){
  Assets.Full = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ Assets.Full = rbind(Assets.Full,Assets[[which[i]]]) } }
  X = UnifySample(Assets.Full[,-1])
  d = length(X[1,]);  sort = msort = 1:d;  ms = Inf;
  for(i in 1:(d-1)){
    msort = sort
    mind = i
    for(k in (i+(i>1)):d){
      if(i > 1){ msort[i:d] = c(sort[k],sort[-c(1:(i-1),k)]) }
      else{ msort = c(sort[k],sort[-k]) }
      print(msort)
      NPC = NPCCombination(Assets,1:2,c(),msort);  s = 0
      for(j in 1:1){ s = s + NPCPortfolio(100000,1,NPC)$value }
      print(c(s,ms))
      if(s<ms){ ms = s;  mind = k }
    }
    print(mind)
    if(i > 1){ sort[i:d] = c(sort[mind],sort[-c(1:(i-1),mind)]) }
    else{ sort = c(sort[mind],sort[-mind]) }
  }
  return(list(sort,ms))
}

OVFind = function(Assets,which,m=100){
  Assets.Full = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ Assets.Full = rbind(Assets.Full,Assets[[which[i]]]) } }
  Y = UnifySample(Assets.Full[,-1])
  d = length(Y[1,]);  VC = matrix(0,d*m-m,d*m-m)
  sort = 1:d;  ms = Inf
  for(i in 1:(d-1)){
    msort = sort;  mind = i
    for(k in (i+(i>1)):d){
      if(i > 1){ msort[i:d] = c(sort[k],sort[-c(1:(i-1),k)]) }
      else{ msort = c(sort[k],sort[-k]) }
      print(msort)
      X = Y
      for(i1 in i:(d-1)){
        for(j in 1:(d-i1)){
          VC[(i1*m-m+1):(i1*m),(j*m-m+1):(j*m)] = NPCbi(X[,msort[c(i1,j+i1)]],m)
          X[,j+i1] = NPCTransform(X[,msort[j+i1]],X[,msort[i1]],VC[(i1*m-m+1):(i1*m),(j*m-m+1):(j*m)])
        }
      }
      X[X>=1] = 0.9999;  X[X<=0] = 0.0001
      s = abs((sum(qnorm(X)^2)-sum(length(X)))/sqrt(2*sum(length(X))))
      print(c(s,ms))
      if(ms > s){ ms = s;  mind = k }
    }
    if(i > 1){ sort[i:d] = c(sort[mind],sort[-c(1:(i-1),mind)]) }
    else{ sort = c(sort[mind],sort[-mind]) }
    for(j in 1:(d-i)){
      VC[(i*m-m+1):(i*m),(j*m-m+1):(j*m)] = NPCbi(Y[,sort[c(i,j+i)]],m)
      Y[,j+i] = NPCTransform(Y[,sort[j+i]],Y[,sort[i]],VC[(i*m-m+1):(i*m),(j*m-m+1):(j*m)])
    }
  }
  return(sort)
}

OVPort = function(steps,Assets,which){
  Assets.Full = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ Assets.Full = rbind(Assets.Full,Assets[[which[i]]]) } }
  d = length(Assets.Full[1,])-1;  S = numeric(d);  k = 0;  P = matrix(0,0,d)
  for(i in 1:length(Assets.Full[,1])){
    S = S + Assets.Full[i,-1];  k = k + 1
    if(k >= steps){
      P = rbind(P,S);  k = 0;  S = numeric(d)
    }
  }
  W = optim(rep(1/d,d-1),ES,control=list(warn.1d.NelderMead=FALSE,maxit=10000),X=P,alpha=0.05)
  return(W)
}


AssetsUnpackList(list(ForecastCombine(30,Assets,3:12,OVPort(10,Assets,1:2))[[1]]),1:10,c("1"))

OVFind(Assets10,1:2)


NPCPortfolio(100000,1,NPCCombination(Assets10,1:2,c(),c(3,5,1,6,2,4,9,7,8,10,11,12)))

ForecastCombine(1,Assets10,3:12,NPCPortfolio(100000,1,NPCCombination(Assets10,1:2,c(),c(3,5,1,6,2,4,9,7,8,10,11,12))))

