### Alpha Diversity

source("Functions.R")

## PND15

#Q1 Does extraction method result in significantly different MB profiles (HFD)
anova(lm(Sobs ~ factor(GroupLabel), data=env1))
bargraph.CI(x.factor = GroupLabel, response = Sobs, legend=TRUE, data=env1)
anova(lm(Chao1 ~ factor(GroupLabel), data=env1))
bargraph.CI(x.factor = GroupLabel, response = Chao1, legend=TRUE, data=env1)
anova(lm(Shan ~ factor(GroupLabel), data=env1))
bargraph.CI(x.factor = GroupLabel, response = Shan, legend=TRUE, data=env1)
anova(lm(PD ~ factor(GroupLabel), data=env1))
bargraph.CI(x.factor = GroupLabel, response = PD, legend=TRUE, data=env1)

#Q2 Does extraction method result in significantly different MB profiles (LFD)
anova(lm(Sobs ~ factor(GroupLabel), data=env2))
bargraph.CI(x.factor = GroupLabel, response = Sobs, legend=TRUE, data=env2)
anova(lm(Chao1 ~ factor(GroupLabel), data=env2))
bargraph.CI(x.factor = GroupLabel, response = Chao1, legend=TRUE, data=env2)
anova(lm(Shan ~ factor(GroupLabel), data=env2))
bargraph.CI(x.factor = GroupLabel, response = Shan, legend=TRUE, data=env2)
anova(lm(PD ~ factor(GroupLabel), data=env2))
bargraph.CI(x.factor = GroupLabel, response = PD, legend=TRUE, data=env2)

#Q3 Does extraction method result in significantly different MB profiles (all)
#   What is the microbiome profile of the May 2021 raw fibre model milk samples
anova(lm(Sobs ~ factor(GroupLabel), data=env3))
bargraph.CI(x.factor = GroupLabel, response = Sobs, legend=TRUE, data=env3)
anova(lm(Chao1 ~ factor(GroupLabel), data=env3))
bargraph.CI(x.factor = GroupLabel, response = Chao1, legend=TRUE, data=env3)
anova(lm(Shan ~ factor(GroupLabel), data=env3))
bargraph.CI(x.factor = GroupLabel, response = Shan, legend=TRUE, data=env3)
anova(lm(PD ~ factor(GroupLabel), data=env3))
bargraph.CI(x.factor = GroupLabel, response = PD, legend=TRUE, data=env3)


## Additional analyses

# All GroupLabels in both time points


bargraph.CI(x.factor = GroupLabel, response = Sobs, GroupLabel = PND, legend=TRUE, data=env)
bargraph.CI(x.factor = GroupLabel, response = Chao1, GroupLabel = PND, legend=TRUE, data=env)
bargraph.CI(x.factor = GroupLabel, response = Shan, GroupLabel = PND, legend=TRUE, data=env)
bargraph.CI(x.factor = GroupLabel, response = PD, GroupLabel = PND, legend=TRUE, data=env)

for(i in colnames(env[,10:13])) {
  print(i)
  print(summary(aov(env[,i] ~ env$GroupLabel * env$PND)))
}


for(i in colnames(env[,10:13])) {
  print(i)
  print(cld(lsmeans(aov(env[,i] ~ GroupLabel * PND, data= env), 
                    ~ GroupLabel * PND, 
                    Letters = letters, 
                    adjust = "tukey")))
}




