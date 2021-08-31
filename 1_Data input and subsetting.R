### Nestle 16S analyses 

# Paul Dennis and Rob Hoelzle Aug 2021

# Data input and subsetting

# Get data

# Get environmental metadata

env.all <- read.csv2('../data_may2021_samples_799/env.csv', header = T, sep = ',', row.names = 1)

# Get OTU table

otu.tmp <- read.table('../data_may2021_samples_799/otu_with_tax_10650.csv', header=TRUE, sep=',', row.names=1)
otu.all <- as.data.frame(t(otu.tmp[,-25]/10650)) # transpose and divide by 4700 to get relative abundances

# Get Taxonomy

taxonomy <- otu.tmp[,25] # Make a list of OTUs with the respective taxonomy
taxonomy <- as.data.frame(taxonomy)
row.names(taxonomy) <- row.names(otu.tmp)
taxonomy$OTU = row.names(taxonomy)
taxonomy$OTU.tax = paste(taxonomy$OTU,taxonomy$taxonomy,sep='.')

# Get weighted unifrac

w_unifrac.all <- read.table('../data_may2021_samples_799/weighted_unifrac.csv', header=TRUE, sep=',', row.names=1)

# Check the samples are in the same order etc...

row.names(otu.all) == row.names(env.all)
row.names(otu.all) == row.names(w_unifrac.all)
rm(otu.tmp) # tidy up

#Subset to questions

env1 <- env.all[env.all$Q1 == 1,]
env2 <- env.all[env.all$Q2 == 1,]
env3 <- env.all[env.all$Q3 == 1,]
env4 <- env.all[env.all$Q4 == 1,]
env5 <- env.all[env.all$Q5 == 1,]

otu1 <- otu.all[env.all$Q1 == 1,]
otu2 <- otu.all[env.all$Q2 == 1,]
otu3 <- otu.all[env.all$Q3 == 1,]
otu4 <- otu.all[env.all$Q4 == 1,]
otu5 <- otu.all[env.all$Q5 == 1,]

w_unifrac1 <- w_unifrac.all[env.all$Q1 == 1,env.all$Q1 == 1]
w_unifrac2 <- w_unifrac.all[env.all$Q2 == 1,env.all$Q2 == 1]
w_unifrac3 <- w_unifrac.all[env.all$Q3 == 1,env.all$Q3 == 1]
