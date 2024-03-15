library(ggplot2)
repro <- read.csv("D:/A Novel Neural Network by Topological Dynamic Space Transform Method/KAILIN/Database/program2k_result.csv")
ggplot(repro,aes(x=Var1,y=Var2,color=Var3))+geom_point()
ggplot(repro,aes(x=Var1,y=Var2,color=Var3))+geom_point()