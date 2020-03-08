#! /urs/bin/env Rscript
# 计算基于otu表的anosim


library("vegan")

distance.bray<-vegdist(FishBio,method = 'bray')

hclust.fish<-hclust(distance.bray,method = "average")

plot(hclust.fish)

distance.bray<-vegdist(FishBio,method = 'bray')

anosim.result<-anosim(distance.bray,FishEnv$Group,permutations = 999)

summary(anosim.result)