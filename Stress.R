
StrsTrunc = function(rho,m=100){
  x = y = qnorm(0.999*(1:100)/100+0.0005)
  M = matrix(0,m,m)
  for(i in 1:m){
    M[i,] = pbinorm(x[i],y,cov12 = rho)
  }
  return(M)
}




StrsVine = function(Cor,n,m=100){
  k = length(Cor[1,]);  PCor = Cor;  VC = NPCVAR(rmvn(n,rep(0,k),Cor))[[3]]
  for(i1 in 2:k){
    for(i2 in i1:k){
      if(i1==i2){ next }
      PCor[i1,i2] = Cor[i1,i2]-PCor[i1-1,i1]*PCor[i1-1,i2]
      PCor[i1,i2] = PCor[i1,i2]/sqrt((1-PCor[i1-1,i1]^2)*(1-PCor[i1-1,i2]^2))
      PCor[i2,i1] = PCor[i1,i2]
    }
  }
  if(F){
    Error = 0
    for(i1 in 2:k){
      for(i2 in 1:(i1-1)){
        C1 = rbind(0,cbind(0,StrsTrunc(PCor[i1,i2])));  C1 = C1[1:m,1:m] + C1[2:(m+1),2:(m+1)] - C1[1:m,2:(m+1)] - C1[2:(m+1),1:m]
        C2 = rbind(0,cbind(0,VC[i2*m-m+1:m,(i1-i2)*m-m+1:m]));  C2 = C2[1:m,1:m] + C2[2:(m+1),2:(m+1)] - C2[1:m,2:(m+1)] - C2[2:(m+1),1:m]
        Error = Error + sum(abs(C1-C2))
      }
    }
    Error = 20000*Error/(k^2-k)/m^2
  }
  if(T){
    X = NPCVineGenerate(100000,VC,m)
    Error = 0
    for(i1 in 2:k){
      for(i2 in 1:(i1-1)){
        C1 = rbind(0,cbind(0,StrsTrunc(Cor[i1,i2])));  C1 = C1[1:m,1:m] + C1[2:(m+1),2:(m+1)] - C1[1:m,2:(m+1)] - C1[2:(m+1),1:m]
        C2 = matrix(0,m,m)
        for(i in 1:100000){
          C2[ceiling(m*X[i,i1]),ceiling(m*X[i,i2])] = C2[ceiling(m*X[i,i1]),ceiling(m*X[i,i2])] + 1
        }
        Error = Error + sum(abs(C1-C2/100000))
      }
    }
    Error = 20000*Error/(k^2-k)/m^2
  }
  return(Error)
}

StrsTest = function(K,N,M,m=100){
  Error = Time = matrix(0,length(K)*length(N),M)
  for(i1 in 1:length(K)){
    l = 0.7*runif(K[i1])+0.1;  v = matrix(runif(K[i1]^2),K[i1],K[i1]);  v = eigen(t(v)%*%v)$vectors
    Cor = v%*%diag(l)%*%t(v);  Cor = round(Cor/(sqrt(diag(Cor))%*%t(sqrt(diag(Cor)))),2)
    for(i2 in 1:length(N)){
      for(i3 in 1:M){
        time = Sys.time()
        Error[length(N)*i1-length(N)+i2,i3] = StrsVine(Cor,N[i2],m)
        t = difftime(Sys.time(), time, units = "secs")
        Time[length(N)*i1-length(N)+i2,i3] = t
        print(c(K[i1],N[i2],i3,t))
      }
    }
  }
  return(list(Error,Time))
}

StrsError = StrsTest(c(3,5,10,25,50,100),c(10000),10)




lm(log(StrsError[[1]]%*%rep(1/10,10)) ~ c(3,5,10,25,50,100))
plot(c(3,5,10,25,50,100),StrsError[[1]]%*%rep(1/10,10),pch=16,ylim=c(0.1,0.2),xlab="Dimensions", ylab="Error", main="Estimation Error for Higher Dimensions")

lines(c(3,5,10,25,50,100),exp(-1.8933101+0.0004198*c(3,5,10,25,50,100)),col="red")




OVPort(10,Assets,1:2)



























