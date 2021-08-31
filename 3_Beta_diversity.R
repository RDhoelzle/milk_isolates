### Beta Diversity

source("Functions.R")


#Q1 Does extraction method result in significantly different MB profiles (HFD)
adonis(sqrt(otu1) ~ GroupLabel, data=env1, method='euc')
  otu1.rda <- rda(otu1 ~ env1$GroupLabel)
    palette(c("black","grey"))  
    custom.plot.pca(ord = otu1.rda, group = factor(env1$GroupLabel), legend.position = "bottomleft",
                    title = "Q1", plottype = "n", sd.val = 10, scaling.val = 3)
    
#Q2 Does extraction method result in significantly different MB profiles (LFD)
adonis(sqrt(otu2) ~ GroupLabel, data=env2, method='euc')
  otu2.rda <- rda(otu2 ~ env2$GroupLabel)
    palette(c("red","orange"))  
    custom.plot.pca(ord = otu2.rda, group = factor(env2$GroupLabel), legend.position = "bottomleft",
                    title = "Q2", plottype = "n", sd.val = 10, scaling.val = 3)
    
#Q3 Does extraction method result in significantly different MB profiles (all)
#   What is the microbiome profile of the May 2021 raw fibre model milk samples
adonis(sqrt(otu3) ~ AdultFDiet * OtherInfo, data=env3, method='euc')
  otu3.rda <- rda(sqrt(otu3) ~ env3$GroupLabel)
    palette(c("black","grey","red","orange"))  
    custom.plot.rda(ord = otu3.rda, group = factor(env3$GroupLabel), legend.position = "bottomleft",
                    title = "Q3", plottype = "n", sd.val = 10, scaling.val = 3)
    