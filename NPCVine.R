
NPCbi = function(X,m=100,b=0.02){
  n = length(X[,1]);  C = matrix(0.001,m,m);  b = b*m
  for(i in 1:n){
    x = max((round(m*X[i,1])-b),1):min((round(m*X[i,1])+b),m);  y = max((round(m*X[i,2])-b),1):min((round(m*X[i,2])+b),m)
    C[x,y] = C[x,y] + 1
  }
  C[b:1,] = C[b:1,]*((1:b)%*%t(rep(1,m))+b+1)/(b+1);  C[(m-b+1):m,] = C[(m-b+1):m,]*((1:b)%*%t(rep(1,m))+b+1)/(b+1)
  C[,b:1] = C[,b:1]*(rep(1,m)%*%t(1:b)+b+1)/(b+1);  C[,(m-b+1):m] = C[,(m-b+1):m]*(rep(1,m)%*%t(1:b)+b+1)/(b+1)
  for(i in 1:m){ C[,i] = cumsum(C[,i]) }
  for(i in 1:m){ C[i,] = cumsum(C[i,]) }
  
  C = C*C[m,m]*((1:m)/C[,m]/m)%*%t((1:m)/C[m,]/m)
  return(C)
}


NPCTransform = function(X,Y,C){
  X[X>0.9999] = 0.9999;  Y[Y>0.9999] = 0.9999#;  X[X<0.000001] = 0.001;  Y[Y<0.000001] = 0.001
  m = length(C[,1]);  n = length(X);  X.T = numeric(n)
  X.I = floor(m*X)+1;  Y.I = floor(m*Y)+1
  dC = (C-rbind(0,C[1:(length(C[1,])-1),]));  dC = cbind(0,dC/(dC[,m]%*%t(rep(1,m))))
  dC[dC<=0] = 0.0001
  Ind = m*X.I + Y.I - m
  X.T = (X.I-m*X)*dC[Ind] + (1+m*X-X.I)*dC[Ind+m]
  X.T[X.T>0.9999] = 0.9999
  return(X.T)
}


NPCVine = function(X,m=100){
  d = length(X[1,])
  VC = matrix(0,d*m-m,d*m-m)
  for(i in 1:(d-1)){
    for(j in 1:(d-i)){
      VC[(i*m-m+1):(i*m),(j*m-m+1):(j*m)] = NPCbi(X[,c(i,j+i)],m)
      X[,j+i] = NPCTransform(X[,j+i],X[,i],VC[(i*m-m+1):(i*m),(j*m-m+1):(j*m)])
    }
  }
  return(VC)
}

NPCInv = function(C){
  m = length(C[,1]);  dC = (C-rbind(0,C[1:(length(C[1,])-1),]));  dC = cbind(0,dC/(dC[,m]%*%t(rep(1,m))))
  dCInv = matrix(0,m,m)
  for(i1 in 1:m){
    I = 1;  i2 = 1
    while(i2 <= m){
      if(dC[i1,1+I] < i2/m){ I = I + 1;  next }
      t = (i2/m-dC[i1,I])/(dC[i1,1+I]-dC[i1,I])
      dCInv[i1,i2] = t + I - 1;  i2 = i2 + 1
    }
  }
  return(cbind(0,dCInv))
}

NPCBiGenerate = function(n,C,X=0){
  m = length(C[,1]);  dCInv = NPCInv(C);  if(length(X)==1){ X = matrix(runif(2*n),n,2) }
  U = floor(m*X[,2]+1);  Ind = m*U + floor(m*X[,1]+1) - m
  X[,2] = ((U-m*X[,2])*dCInv[Ind] + (1+m*X[,2]-U)*dCInv[Ind+m])/m
  return(X)
}

NPCVineGenerate = function(n,VC,m){
  d = length(VC[,1])/m + 1;  X = matrix(0,n,d);  U = matrix(runif(d*n),n,d);  X[,1] = U[,1]
  for(i in 2:d){
    for(j in (i-1):1){
      U[,i] = NPCBiGenerate(n,VC[(j*m-m+1):(j*m),((i-j)*m-m+1):((i-j)*m)],U[,c(j,i)])[,2]
    }
    X[,i] = U[,i]
    if(i == d){ next }
    for(j in 1:(i-1)){
      U[,i] = NPCTransform(U[,i],U[,j],VC[(j*m-m+1):(j*m),((i-j)*m-m+1):((i-j)*m)])
    }
  }
  return(X)
}

NPCPredict = function(n,vals,VC,m){
  d = length(VC[,1])/m + 1;  k = length(vals[1,]);  X = U = cbind(vals,matrix(runif((d-k)*n),n,d-k))
  for(i in 2:d){
    if(i>k){
      for(j in (i-1):1){
        U[,i] = NPCBiGenerate(n,VC[(j*m-m+1):(j*m),((i-j)*m-m+1):((i-j)*m)],U[,c(j,i)])[,2]
      }
      X[,i] = U[,i]
    }
    if(i == d){ next }
    for(j in 1:(i-1)){
      U[,i] = NPCTransform(U[,i],U[,j],VC[(j*m-m+1):(j*m),((i-j)*m-m+1):((i-j)*m)])
    }
  }
  return(X)
}

NPCVineGenerate = function(n,VC,m){
  d = length(VC[,1])/m + 1;  X = matrix(0,n,d);  U = matrix(runif(d*n),n,d);  X[,1] = U[,1]
  for(i in 2:d){
    X[,i] = U[,i]
    for(j in (i-1):1){
      X[,i] = NPCBiGenerate(n,VC[(j*m-m+1):(j*m),((i-j)*m-m+1):((i-j)*m)],cbind(U[,j],X[,i]))[,2]
    }
  }
  return(X)
}

NPCPredict = function(n,vals,VC,m){
  d = length(VC[,1])/m + 1;  k = length(vals[1,]);  X = U = cbind(vals,matrix(runif((d-k)*n),n,d-k))
  if(k > 1){
    for(i in 2:k){
      for(j in 1:(i-1)){
        U[,i] = NPCTransform(U[,i],U[,j],VC[(j*m-m+1):(j*m),((i-j)*m-m+1):((i-j)*m)])
      }
    }
  }
  for(i in (k+1):d){
    X[,i] = U[,i]
    for(j in (i-1):1){
      X[,i] = NPCBiGenerate(n,VC[(j*m-m+1):(j*m),((i-j)*m-m+1):((i-j)*m)],cbind(U[,j],X[,i]))[,2]
    }
  }
  return(X)
}


image(NPCVine(pnorm(rmvn(10000,mu=rep(0,3),Sigma=matrix(c(1,0.5,0.4,0.5,1,0.7,0.4,0.7,1),3,3)))))

plot(NPCVineGenerate(10000,NPCVine(pnorm(rmvn(10000,mu=rep(0,3),Sigma=matrix(c(1,0.5,0.4,0.5,1,0.9,0.4,0.9,1),3,3)))),100)[,2:3])

plot(pnorm(rmvn(10000,mu=rep(0,3),Sigma=matrix(c(1,0.5,0.4,0.5,1,0.9,0.4,0.9,1),3,3)))[,2:3])


I



