GCVAR = function(Assets){
  D = Marginals(Assets);  ID = InvMarginals(Assets)
  U = UnifySample(Assets);  VC = 2*sin(cor(U)*pi/6)
  return(list(D,ID,VC))
}

GCSample = function(n,steps,GC){
  ID = GC[[2]];  GC = GC[[3]];  d = length(GC[1,])
  Dist = matrix(0,n,d)
  for(i in 1:steps){ dist = pnorm(rmvn(n,mu=rep(0,d),Sigma=GC));  Dist = Dist + SampleMarginals(dist,ID) }
  return(Dist)
}

GCPortfolio = function(n,steps,GC){
  X = GCSample(n,steps,GC);  d = length(X[1,])
  w = optim(rep(1/d,d-1),ES,control=list(warn.1d.NelderMead=FALSE),X=X,alpha=0.05)
  return(w)
}

GCCombination = function(Assets,which){
  Assets.Full = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ Assets.Full = rbind(Assets.Full,Assets[[which[i]]]) } }
  return(GCVAR(Assets.Full[,-1]))
}


