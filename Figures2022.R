
#Fi Vin Gau
GaussToVine()

#ta CDF Shape
AssetsStats(Assets,1:2,I)

#ta CDF Shape 10
AssetsStats(Assets10,1:2,I)

#ta CDF Shape Test
AssetsStats(Assets,3:12,I)

#ta CDF Shape Boot
set.seed(2)
BootstrapStats(Assets,1:2,3:12,I,1000)



par(mfrow=c(1,2))
#Fi CDFSPY
hist(c(Assets[[1]][,2],Assets[[2]][,2]),100,main="Minute log-returns for SPY (Training Period)",xlab="Log-Returns")
plot(GC[[1]][,c(1,2)],type="l",main="Sample CDF for SPY (Training Period)",xlab="x",ylab="F(x)")
#Fi CDFNDAQ
hist(c(Assets[[1]][,3],Assets[[2]][,3]),100,main="Minute log-returns for NDAQ (Training Period)",xlab="Log-Returns")
plot(GC[[1]][,c(3,4)],type="l",main="Sample CDF for NDAQ (Training Period)",xlab="x",ylab="F(x)")
par(mfrow=c(1,1))

#ta Con GC
round(GC[[3]][-1,-12],2)*lower.tri(GC[[3]][-1,-12],TRUE)
#Ta Con GC10
round(GC10[[3]][-1,-12],2)*lower.tri(GC10[[3]][-1,-12],TRUE)


#Fi NPCKer
KernelVisualization(function(x,y){ max(0.02-max(abs(x-0.2),abs(y-0.4)),0) },1000)

#Fi NPCBi
NPCVisDensity(Assets,1:2,c(1,2,8,12,3,7,5,9),I)
#Fi NPCBiAR
NPCVisDensity(Assets,1:2,c(1,13,2,14,5,17,11,23),c(IAR,I),12)
#Fi NPCBi 10
NPCVisDensity(Assets10,1:2,c(1,2,8,12,3,7,5,9),I)

#Fi NPCSim
NPCVisSim(Assets,1:2,NPC,c(1,2,8,12,3,7,5,9),I)
#Fi NPCARSim
NPCVisSim(Assets,1:2,NPCAR,c(1,2,8,12,3,7,5,9),I)
#Fi NPCSim10
NPCVisSim(Assets10,1:2,NPC10,c(1,2,8,12,3,7,5,9),I)


#Fi NPCSimSPY
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


#Fi NPCSimTSLA
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

#Ta  Res Port & Ta Res ES
TablePortfolios(list(W.OP,W.GC10,W.NPC10,W.GC,W.NPC,W.NPCAR),list(P.OP,P.GC10,P.NPC10,P.GC,P.NPC,P.NPCAR))

#Ta Res AR Mean
PortDist(E.NPCAR)

#Ta Res Test Rev
PortId(W.RNPCAR,E.NPCAR)

#ta Res Port Rev & ta Res ES Rev
TablePortfolios(list(W.NPCAR,W.RNPCAR),list(P.NPCAR,P.RNPCAR))

#Fi ResHis
AssetsUnpackList(list(P.GC10[[1]],P.NPC10[[1]],P.GC[[1]],P.NPC[[1]],P.NPCAR[[1]],P.RNPCAR[[1]]),
                 1:10,c("GC10","NPC10","GC","NPC","NPC-AR","R-NPC-AR"))

#Ta AR ES & Fi AR ES His
{
  for(i in 1:10){ print(1000*ESTest(P.NPCAR[[1]],i,W.NPCAR)) }
  1000*ESTest(P.NPCAR[[1]],1:10,W.NPCAR,TRUE)
}
1000*ESTest(P.NPC10[[1]],1:10,W.NPC10)


#Ta AR Or ES
AssetsUnpackList(list(P.NPCAR[[1]],P.NPCAR2[[1]],P.NPCAR3[[1]],P.NPCAR4[[1]]),1:10,c("AR 1","AR 2","AR 3","AR 4"))



#NPCAR time horizon
AssetsUnpackList(list(P.NPCARt30[[1]]),1:10,c("NPCAR t30"))
AssetsUnpackList(list(P.NPCARt60[[1]]),1:10,c("NPCAR t60"))
AssetsUnpackList(list(P.NPCARt120[[1]]),1:10,c("NPCAR t120"))

#Ta AR Time Port
TablePortfolios(list(W.NPCAR,W.NPCARt30,W.NPCARt60,W.NPCARt120),list(P.NPCAR,P.NPCARt30,P.NPCARt60,P.NPCARt120))

#ta AR Time Chi
PortId(W.NPCARt30,E.NPCAR)
PortId(W.NPCARt60,E.NPCAR)
PortId(W.NPCARt120,E.NPCAR)

#ta AR Time ES
{
  print(1000*ESTest(P.NPCARt30[[1]],1:10,W.NPCARt30,TRUE))
  print(1000*ESTest(P.NPCARt60[[1]],1:10,W.NPCARt60,TRUE))
  print(1000*ESTest(P.NPCARt120[[1]],1:10,W.NPCARt120,TRUE))
}

#
AssetARMissing(Assets,1:2)
AssetARMissing(Assets,3:12)


chol(matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),3,3))



