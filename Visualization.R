
NPCVisDensity = function(Assets,which,show,names,ar=c()){
  X = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ X = rbind(X,Assets[[which[i]]]) } }
  n = length(X[,1]);  X = cbind(X[-1,1],X[-n,-1],X[-1,-1]);  X = UnifySample(X)
  s = length(show)/2;  m = 100;  z = c(-0.0001,0.002)
  par(mfrow=c(ceiling(sqrt(s)),2*floor(sqrt(s))));  par(mar = c(4.6, 4.1, 1.6, 0.6))
  for(i in 1:s){
    C = NPCbi(X[,1+c(show[2*i-1],show[2*i])],m,b=0.02)
    C = C-cbind(0,C[,-m]);  C = C-rbind(0,C[-m,])
    print(c(max(C),min(C)))
    plot(X[,1+c(show[2*i-1],show[2*i])],cex=0.5,xlab=names[show[2*i-1]],ylab=names[show[2*i]])
    image(C,zlim=z,xlab=names[show[2*i-1]],ylab=names[show[2*i]])
    print(sum((C-0.0001)^2))
  }
  par(mfrow=c(1,1));  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

NPCVisSim = function(Assets,which,NPC,show,names){
  X = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ X = rbind(X,Assets[[which[i]]]) } }
  n = length(X[,1]);  X = cbind(X[-1,1],X[-n,-1],X[-1,-1]);  X = UnifySample(X)
  Y = NPCSample(n,1,NPC)[,(sum(NPC[[5]][-1])+1):sum(NPC[[5]])]
  Y = UnifySample(Y);  s = length(show)/2
  par(mfrow=c(ceiling(sqrt(s)),2*floor(sqrt(s))));  par(mar = c(4.6, 4.1, 1.6, 0.6))
  for(i in 1:s){
    plot(X[,1+c(show[2*i-1],show[2*i])],cex=0.5,xlab=names[show[2*i-1]],ylab=names[show[2*i]])
    plot(Y[,c(show[2*i-1],show[2*i])],cex=0.5,xlab=paste0("Sim ",names[show[2*i-1]]),ylab=paste0("Sim ",names[show[2*i]]))
  }
  par(mfrow=c(1,1));  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

AssetsStats = function(Assets,which,names){
  X = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ X = rbind(X,Assets[[which[i]]]) } }
  n = length(X[,1]);  k = length(X[1,-1])
  Stats = matrix(0,4,k);  colnames(Stats) = I;  rownames(Stats) = c("Mean","Std. Dev","Skewness","Kurtosis")
  for(i in 1:k){
    mu = mean(X[,1+i]);  sig = var(X[,1+i])
    Stats[,i] = round(c(10000*mu,10000*sqrt(sig),sum((X[,1+i]-mu)^3)/n/sig^(3/2),sum((X[,1+i]-mu)^4)/n/sig^2),2)
  }
  return(Stats)
}

AssetsUnpackList = function(Assets,which,names){
  k = length(Assets);  X = matrix(0,0,k)
  for(i in which){
    Y = matrix(0,length(Assets[[1]][[i]]),k)
    for(j in 1:k){ Y[,j] = Assets[[j]][[i]] }
    X = rbind(X,Y)
  }
  S = matrix(0,3,k)
  par(mfrow=c(floor(sqrt(k)),max(1,ceiling(sqrt(k)))));  n = length(X[,1])
  for(i in 1:k){
    hist(X[,i],xlab="Log-Return",main=paste0("Histogram of returns for ",names[i]),100)
    mu = mean(X[,i]);  var = sort(X[,i])[floor(0.05*n)];  es = mean(sort(X[,i])[1:floor(0.05*n)])
    lines(c(mu,mu),c(0,10000),lty=2,col="green")
    lines(c(var,var),c(0,10000),lty=2,col="blue")
    lines(c(es,es),c(0,10000),lty=2,col="red")
    S[,i] = c(mu,-var,-es)
  }
  par(mfrow=c(1,1))
  print(length(X[,1]))
  return(round(1000*S,2))
}

KernelVisualization = function(K,res){
  V = matrix(0,res+1,res+1)
  for(i1 in 0:res){
    for(i2 in 0:res){
      V[i1,i2] = K(i1/res,i2/res)
    }
  }
  V[V==0] = NA
  image(V,col = hcl.colors(100, "terrain"),main="Influence of a Single Point",xlab="x",ylab="y")
}

TablePortfolios = function(Ws,Ps){
  n = length(Ws)
  Ports = matrix(0,n,length(Ws[[1]]$par)+2)
  for(i in 1:n){ Ports[i,] = c(Ws[[i]]$par,1-sum(Ws[[i]]$par),1000*Ws[[i]]$value) }
  n = length(Ps)
  Es = matrix(0,n,length(Ps[[1]][[1]])+1)
  for(i in 1:n){ Es[i,] = c(1000*Ps[[i]][[2]][[1]],1) }
  
  return(list(list(round(Ports,2)),list(round(Es,2))))
}

AssetARMissing = function(Assets,which){
  Assets.Full = Assets[[which[1]]]
  if(length(which)>1){ for(i in 2:length(which)){ Assets.Full = rbind(Assets.Full,Assets[[which[i]]]) } }
  Time = diff(Assets.Full[,1])
  m = min(Time)
  return(c(length(Time)+1,sum(Time==m),mean(Time==m)))
}

BootstrapStats = function(Assets,training,test,names,m){
  Stats = AssetsStats(Assets,training,names)
  n2 = length(Assets[[training[1]]][,1])
  if(length(training)>1){ for(i in 2:length(training)){ n2 = n2 + length(Assets[[training[i]]][,1]) } }
  
  X = Assets[[test[1]]]
  if(length(test)>1){ for(i in 2:length(test)){ X = rbind(X,Assets[[test[i]]]) } }
  n = length(X[,1]);  k = length(X[1,-1])
  E = matrix(0,m*length(Stats[,1]),length(Stats[1,]))
  Y = matrix(0,n2,k+1)
  for(i in 1:m){
    for(j in 1:k){ Y[,1+j] = sample(X[,1+j],n2,replace=T) }
    E[4*i-4+1:4,] = AssetsStats(list(Y),1,names)
    print(i)
  }
  Conf = rbind(Stats,Stats)
  Test = Stats;  Test[1:4,1:k] = 0
  for(i in 1:k){
    for(j in 1:4){
      conf = sort(E[4*(1:m)-4+j,i]);  conf = conf[c(ceiling(0.025*m),ceiling(0.975*m))]
      Conf[4*(1:2)-4+j,i] = conf
      if(Stats[j,i]<conf[1]){ Test[j,i] = -1 }
      if(Stats[j,i]>conf[2]){ Test[j,i] = 1 }
    }
  }
  return(list(Conf,Test))
}

GaussToVine = function(){
  set.seed(2)
  X = pnorm(rmvn(10000,mu=rep(0,3),Sigma=matrix(c(1,-0.5,0,-0.5,1,0.8,0,0.8,1),3,3)))
  Y = NPCVineGenerate(10000,NPCVine(X),100)
  par(mfrow=c(2,3));  par(mar = c(4.6, 4.1, 1.6, 0.6))
  plot(X[,c(1,2)],cex=0.5,xlab="Dim 1",ylab="Dim 2")
  plot(X[,c(1,3)],cex=0.5,xlab="Dim 1",ylab="Dim 3")
  plot(X[,c(2,3)],cex=0.5,xlab="Dim 2",ylab="Dim 3")
  plot(Y[,c(1,2)],cex=0.5,xlab="Simulated Dim 1",ylab="Simulated Dim 2")
  plot(Y[,c(1,3)],cex=0.5,xlab="Simulated Dim 1",ylab="Simulated Dim 3")
  plot(Y[,c(2,3)],cex=0.5,xlab="Simulated Dim 2",ylab="Simulated Dim 3")
  par(mfrow=c(1,1));  par(mar = c(5.1, 4.1, 4.1, 2.1))
}


