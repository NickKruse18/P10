NPCVARVine = function(Assets,ar=c(),sort=c(),m=100){
  d = length(Assets[1,]);  n = length(Assets[,1])-length(ar)
  if(length(sort)==0){ sort = 1:d };  d = length(sort)
  X = matrix(0,n,d+sum(ar));  ar = c(d,ar);  i1 = 0
  for(j in length(ar):1){ if(ar[j]==0){ next }
    for(i in 1:ar[j]){ i1 = i1 + 1;  X[,i1] = Assets[1:n+length(ar)-j,sort[i]] }
  }
  VC = NPCVine(X,m=m)
  return(VC)
}

NPCVAR = function(Assets,ar=c(),sort=c(),m=100){
  d = length(Assets[1,])
  if(length(sort)==0){ sort = 1:d };  d = length(sort)
  D = Marginals(matrix(Assets[,sort[1:d]],length(Assets[,1]),d))
  ID = InvMarginals(matrix(Assets[,sort[1:d]],length(Assets[,1]),d))
  A.U = UnifySample(Assets)
  VC = NPCVARVine(A.U,ar,sort,m)
  ar = c(d,ar)
  return(list(D,ID,VC,m,ar,sort))
}

NPCSample = function(n,steps,NPC,vals=NULL){
  if(is.null(vals)){ vals = matrix(runif(n),n,1) }
  if(!is.matrix(vals)){ vals = rep(1,n)%*%t(vals) }
  ID = NPC[[2]];  VC = NPC[[3]];  m = NPC[[4]];  ar = NPC[[5]]
  Dist = matrix(0,n,sum(ar));  dist = NPCPredict(n,vals,VC,m);  i1 = 0
  for(j in length(ar):1){
    for(i in 1:ar[j]) { i1 = i1 + 1;  Dist[,i1] = Dist[,i1] + SampleMarginals(matrix(dist[,i1],n,1),ID[,(2*i-1):(2*i)]) } }
  if(steps==1){ return(Dist) }
  for(i in 2:steps){
    if(length(ar)>1){ i1 = 0;  ndist = matrix(0,n,sum(ar[-1]))
    for(j in length(ar):2){ for(i in 1:ar[j]) { i1 = i1 + 1;  ndist[,i1] = dist[,i1+ar[j]] } }
    dist = NPCPredict(n,ndist,VC,m)
    }
    else{ dist = NPCVineGenerate(n,VC,m) }
    i1 = 0
    for(j in length(ar):1){
      for(i in 1:ar[j]) { i1 = i1 + 1;  Dist[,i1] = Dist[,i1] + SampleMarginals(matrix(dist[,i1],n,1),ID[,(2*i-1):(2*i)]) } }
  }
  return(Dist)
}

NPCPortfolio = function(n,steps,NPC,vals=NULL){
  time = Sys.time()
  ar = NPC[[5]];  X = NPCSample(n,steps,NPC,vals)[,(sum(ar[-1])+1):sum(ar)]
  print(Sys.time() - time)
  
  X = X[,match(1:ar[1],NPC[[6]])]
  w = optim(rep(1/ar[1],ar[1]-1),ES,control=list(warn.1d.NelderMead=FALSE,maxit=500),X=X,alpha=0.05)
  return(w)
}

NPCCombination = function(Assets,which,ar=c(),sort=c(),m=100){
  Assets.Full = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ Assets.Full = rbind(Assets.Full,Assets[[which[i]]]) } }
  return(NPCVAR(Assets.Full[,-1],ar,sort,m))
}


