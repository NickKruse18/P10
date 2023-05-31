

Assets21 = AssetsCombine(AssetsRead("21"))
Assets2110 = AssetsSum(Assets21,10)


NPC2110 = NPCCombination(Assets2110,1:2,c(),1:12)
W.NPC2110 = NPCPortfolio(100000,1,NPC2110)
P.NPC2110 = ForecastCombine(1,Assets2110,3:9,W.NPC2110)


NPC21AR = NPCCombination(Assets21,1:2,c(12),1:12)
W.NPC21AR = NPCPortfolio(100000,10,NPC21AR)
P.NPC21AR = ForecastCombine(10,Assets21,3:9,W.NPC21AR)

#Control Optimal Portfolio
W.OP21 = OVPort(10,Assets21,3:9)
P.OP21 = ForecastCombine(10,Assets21,3:9,W.OP21)


#ta CDF Shape 21
AssetsStats(Assets21,1:2,I)

#Fi NPCBi 21
NPCVisDensity(Assets21,1:2,c(1,2,8,12,3,7,5,9),I)
GCCombination(Assets21,1:2)

#Ta  Res Port 21 & Ta Res ES 21
TablePortfolios(list(W.OP21,W.NPC2110,W.NPC21AR),list(P.OP21,P.NPC2110,P.NPC21AR))
AssetsUnpackList(list(P.NPC2110[[1]],P.NPC21AR[[1]]),1:7,c("NPC-10","NPC-AR"))
AssetsUnpackList(list(P.NPC2110[[1]],P.NPC21AR[[1]]),1:6,c("NPC-10","NPC-AR"))

#Ta Res Test 21
PortId(W.NPC21AR,E.NPCAR)

#ta CDF Shape 21 dec
AssetsStats(Assets21,9,I)

#Ta ES test 21
1000*ESTest(P.NPC2110[[1]],1:6,W.NPC2110)
1000*ESTest(P.NPC21AR[[1]],1:6,W.NPC21AR)




