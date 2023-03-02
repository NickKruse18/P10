
AssetSum = function(Assets,which,steps,show){
  X = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ X = rbind(X,Assets[[which[i]]]) } }
  n = length(X[,1]);  s = length(show);  I = diag(ceiling(n/steps))
  I = t(kronecker(I,rep(1,steps))[1:n,]);  X = I%*%X
  par(mfrow=c(floor(sqrt(s)),ceiling(sqrt(s))))
  for(i in 1:s){
    hist(X[,1+show[i]],100)
  }
  par(mfrow=c(1,1))
}

NPCVisDensity = function(Assets,which,show,names,ar=c()){
  X = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ X = rbind(X,Assets[[which[i]]]) } }
  n = length(X[,1]);  X = cbind(X[-1,1],X[-n,-1],X[-1,-1]);  X = UnifySample(X)
  s = length(show)/2;  m = 100
  par(mfrow=c(ceiling(sqrt(s)),2*floor(sqrt(s))));  par(mar = c(4.6, 4.1, 0.6, 0.6))
  for(i in 1:s){
    C = NPCbi(X[,1+c(show[2*i-1],show[2*i])],m)
    C = C-cbind(0,C[,-m]);  C = C-rbind(0,C[-m,])
    plot(X[,1+c(show[2*i-1],show[2*i])],cex=0.5,xlab=names[show[2*i-1]],ylab=names[show[2*i]])
    image(C,col = hcl.colors(100, "terrain"),xlab=names[show[2*i-1]],ylab=names[show[2*i]])
  }
  par(mfrow=c(1,1));  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

#Fig NPCBi
NPCVisDensity(Assets,1:2,c(1,2,8,12,3,7,5,9),I)
#Fig NPCBiAR
NPCVisDensity(Assets,1:2,c(1,13,2,14,5,17,11,23),c(IAR,I),12)

AssetSum(Assets,1:2,10,c(1))

?image

hist(NPCSample(1000,10,NPCAR)[,13],100)
hist(NPCSample(1000,10,NPC)[,1],100)



#Fig NPCSimSPY
{
  par(mfrow=c(2,2))
  X = NPCSample(100000,1,NPCAR)
  hist(X[,13],100,main="Simulated log-returns for SPY",xlab="log-return")
  
  X = NPCSample(100000,1,NPCAR,vals = 0.5)
  hist(X[,13],100,main="Simulated log-returns for SPY, given 0.5",xlab="log-return")
  
  X = NPCSample(100000,1,NPCAR,vals = 0.01)
  hist(X[,13],100,main="Simulated log-returns for SPY, given 0.01",xlab="log-return")
  
  X = NPCSample(100000,1,NPCAR,vals = 0.99)
  hist(X[,13],100,main="Simulated log-returns for SPY, given 0.99",xlab="log-return")
}


#Fig NPCSimTSLA
{
  par(mfrow=c(2,2))
  X = NPCSample(100000,1,NPCCombination(Assets,1:2,c(1),5))
  hist(X[,2],100,main="Simulated log-returns for TSLA",xlab="log-return")
  
  X = NPCSample(100000,1,NPCCombination(Assets,1:2,c(1),5),vals = 0.5)
  hist(X[,2],100,main="Simulated log-returns for TSLA, given 0.5",xlab="log-return")
  
  X = NPCSample(100000,1,NPCCombination(Assets,1:2,c(1),5),vals = 0.01)
  hist(X[,2],100,main="Simulated log-returns for TSLA, given 0.01",xlab="log-return")
  
  X = NPCSample(100000,1,NPCCombination(Assets,1:2,c(1),5),vals = 0.99)
  hist(X[,2],100,main="Simulated log-returns for TSLA, given 0.99",xlab="log-return")
}


















