
LinInv = function(x,m=1000){
  y = 1:length(x);  y = (y-y[1]) / (y[length(y)]-y[1])
  C = matrix(0,m+1,2);  C[,1] = (0:m)/m;  i = 0;  j = 2
  while(j <= length(y)){
    if(i/m-0.0000001 >= y[j]){ j = j + 1;  next; }
    C[i+1,2] = (i/m-y[j-1])*(x[j]-x[j-1])/(y[j]-y[j-1]) + x[j-1]
    i = i + 1
  }
  return(C)
}

Lin = function(x,m=1000){
  x = sort(x);  C = matrix(0,m+1,2);  C[,1] = (0:m)/(m+0.0001)*(x[length(x)]-x[1])+x[1];  j = 2;  i = 1
  while(i <= length(C[,1])){
    if(C[i,1] > x[j]){ j = j + 1;  next; }
    C[i,2] = ((C[i,1]-x[j-1])/(x[j]-x[j-1]) + j-1)/length(x)
    i = i + 1
  }
  return(C)
}

UnifySample = function(X){
  for(i in 1:length(X[1,])){
    X[,i] = (match(X[,i],sort(X[,i]))-1)/length(X[,i])
  }
  return(X)
}

Marginals = function(X,m=1000){
  d = length(X[1,]);  D = matrix(0,m+1,2*d)
  for(i in 1:d){
    D[,(2*i-1):(2*i)] = Lin(X[,i],m)
  }
  return(D)
}

InvMarginals = function(X,m=1000){
  d = length(X[1,]);  D = matrix(0,m+1,2*d)
  for(i in 1:d){
    D[,(2*i-1):(2*i)] = LinInv(sort(X[,i]),m)
  }
  return(D)
}

SampleMarginals = function(U,D){
  d = length(U[1,]);  m = length(D[,1])-1;  M = U;  U = m*U+1;  U.I = floor(U)
  for(i in 1:d){
    M[,i] = (U[,i]-U.I[,i])*D[U.I[,i],2*i] + (1-U[,i]+U.I[,i])*D[U.I[,i]+1,2*i]
  }
  return(M)
}

TransformMarginals = function(X,D){
  d = length(X[1,]);  m = length(D[,1])-1;  M = X;  X = m*t((t(X)-D[1,2*(1:d)-1])/(D[m,2*(1:d)-1]-D[1,2*(1:d)-1]))+1
  X = pmax(pmin(X, 1000), 1);  X.I = floor(X)
  for(i in 1:d){
    M[,i] = (X[,i]-X.I[,i])*D[X.I[,i],2*i] + (1-X[,i]+X.I[,i])*D[X.I[,i]+1,2*i]
  }
  return(M)
}

plot(NPCAR[[1]][,1:2],type="l")
hist(Assets[[1]][,2])
TransformMarginals(Assets[[1]][,-1],NPCAR[[1]])[,1]
hist(TransformMarginals(Assets[[2]][,-1],NPCAR[[1]])[,4])

KTau = function(X){
  X = X[order(X[,1]),];  n = length(X[,1]);  tau = 0
  for(i in 2:n){ tau = tau + sum(X[i,2]>X[1:(i-1),2]) }
  return(4*tau/(n*(n-1))-1)
}

CopIndep = function(X){
  n = length(X[,1]);  tau = 3*sqrt(n*(n-1)/(2*(2*n+5)))*abs(KTau(X))
  return(2-2*pnorm(tau))
}




