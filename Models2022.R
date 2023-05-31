



#One minute data
Assets = AssetsCombine(AssetsRead("22"))
#Ten minute data
Assets10 = AssetsSum(Assets,10)

I = c("SPY","NDAQ","AMD","AAPL","TSLA","NVO","NVDA","BA","LMT","GD","JPM","WFC")
IAR = c("Lag SPY","Lag NDAQ","Lag AMD","Lag AAPL","Lag TSLA","Lag NVO","Lag NVDA","Lag BA","Lag LMT","Lag GD","Lag JPM","Lag WFC")


#Gaussian Copula model
GC = GCCombination(Assets,1:2)
W.GC = GCPortfolio(100000,10,GC)
P.GC = ForecastCombine(10,Assets,3:12,W.GC)

#Non-Parametric Copula model
NPC = NPCCombination(Assets,1:2,c(),1:12)
W.NPC = NPCPortfolio(100000,10,NPC)
P.NPC = ForecastCombine(10,Assets,3:12,W.NPC)

#ten minute Gaussian Copula model
GC10 = GCCombination(Assets10,1:2)
W.GC10 = GCPortfolio(100000,1,GC10)
P.GC10 = ForecastCombine(1,Assets10,3:12,W.GC10)

#ten minute Non-Parametric Copula model
NPC10 = NPCCombination(Assets10,1:2,c(),1:12)
W.NPC10 = NPCPortfolio(100000,1,NPC10)
P.NPC10 = ForecastCombine(1,Assets10,3:12,W.NPC10)

#Non-Parametric Copula model with Assets lag 1 included
NPCAR = NPCCombination(Assets,1:2,c(12),1:12)
W.NPCAR = NPCPortfolio(100000,10,NPCAR)
P.NPCAR = ForecastCombine(10,Assets,3:12,W.NPCAR)

#Control Optimal Portfolio
W.OP = OVPort(10,Assets,3:12)
P.OP = ForecastCombine(10,Assets,3:12,W.OP)

#Reversed Asset order of NPCAR
RNPCAR = NPCCombination(Assets,1:2,c(12),12:1)
W.RNPCAR = NPCPortfolio(100000,10,RNPCAR)
P.RNPCAR = ForecastCombine(10,Assets,3:12,W.RNPCAR)

#Monte Carlo error estimates for NPCAR portfolio
E.NPCAR = NPCError(100000,10,100,NPCAR)

#Ta Res Test Rev
PortDist(E.NPCAR)


#Higher Orders
NPCAR2 = NPCCombination(Assets,1:2,c(12,12),1:12)
W.NPCAR2 = NPCPortfolio(100000,10,NPCAR2)
P.NPCAR2 = ForecastCombine(10,Assets,3:12,W.NPCAR2)

NPCAR3 = NPCCombination(Assets,1:2,c(12,12,12),1:12)
W.NPCAR3 = NPCPortfolio(100000,10,NPCAR3)
P.NPCAR3 = ForecastCombine(10,Assets,3:12,W.NPCAR3)

NPCAR4 = NPCCombination(Assets,1:2,c(12,12,12,12),1:12)
W.NPCAR4 = NPCPortfolio(100000,10,NPCAR4)
P.NPCAR4 = ForecastCombine(10,Assets,3:12,W.NPCAR4)

#Time difference of 58.84316 secs
#Time difference of 1.887792 mins
#Time difference of 2.713561 mins


#Greater Time Horizon
W.NPCARt30 = NPCPortfolio(100000,30,NPCAR)
P.NPCARt30 = ForecastCombine(30,Assets,3:12,W.NPCARt30)

W.NPCARt60 = NPCPortfolio(100000,60,NPCAR)
P.NPCARt60 = ForecastCombine(60,Assets,3:12,W.NPCARt60)

W.NPCARt120 = NPCPortfolio(100000,120,NPCAR)
P.NPCARt120 = ForecastCombine(120,Assets,3:12,W.NPCARt120)

W.NPCARt720 = NPCPortfolio(100000,720,NPCAR)
P.NPCARt720 = ForecastCombine(720,Assets,3:12,W.NPCARt720)

#Time difference of 1.64736 mins
#Time difference of 3.375557 mins
#Time difference of 6.383673 mins

