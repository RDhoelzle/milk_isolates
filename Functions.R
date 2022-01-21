#  _____                 _   _                 
# |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#
# Paul Dennis & Robert Hoelzle

library(vegan)
library(sciplot)
library(multcomp)
library(lsmeans)
library(multcompView)
library(corrplot)
library(mvabund)
library(pheatmap)
library(RColorBrewer)
library(genoPlotR)

## Axis percent 
# Returns the percentage varation on ordination axes

axis.percent <- function(ordination){
  round((100*eigenvals(ordination)[1:2]/ordination$tot.chi[[1]]),digits=2)
}

# END



## pcoa
# Makes a PCoA ordination object
pcoa <- function(ord){
  cmdscale(dist(ord,'euc'), k=(nrow(ord)-1), eig=TRUE)
}

## axis.percent.pcoa 
# Returns the percentage varation on ordination axes
axis.percent.pcoa <- function(ordination){
  round(as.vector(eigenvals(ordination)/sum(eigenvals(ordination)))[1:2]*100,digits=2)
}

## custom.plot.pcoa
# Plots a PCoA ordination with group enrichment

custom.plot.pcoa <- function(ord,group,title,plottype){
  ordiplot(scores(ord, choices=c(1,2)), type=plottype, 
           xlab=paste("PCoA1 (",axis.percent.pcoa(
             ord)[[1]],"%)",sep=""),
           ylab=paste("PCoA2 (",axis.percent.pcoa(ord)[[2]],"%)",sep=""),
           main=title)
  ordiellipse(ord,group,
              kind='sd')
  #conf=0.95)
  points(ord$points, pch=21, cex=3, bg=factor(group))    
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  legend("topright",legend=unique(factor(group)),pch=19,col=unique(factor(group)))
}

## custom.plot.II.pcoa
# Plots a PCoA ordination with group enrichment

custom.plot.II.pcoa <- function(ord,colours,ellipses,size,title,plottype, legend.position){
  ordiplot(scores(ord, choices=c(1,2)), type=plottype, 
           xlab=paste("PCoA1 (",axis.percent.pcoa(
             ord)[[1]],"%)",sep=""),
           ylab=paste("PCoA2 (",axis.percent.pcoa(ord)[[2]],"%)",sep=""),
           main=title)
  ordiellipse(ord,ellipses,
              kind='sd', lty=4)
  #conf=0.95)
  points(ord$points, pch=21, cex=size, bg=factor(colours))    
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  legend(legend.position,legend=unique(factor(colours)),pch=19,col=unique(factor(colours)))
}

## plot soil PCA - plots a pregenerated PCA ordination with groups for ellipses	

custom.plot.pca <- function(ord,group,title,plottype,legend.position,sd.val,scaling.val){
  plot(ord, type=plottype, scaling=scaling.val, 
       #xlim=c(-0.3,0.35), ylim=c(-0.4,0.4),
       xlab=paste("PC1 (",axis.percent(
         ord)[[1]],"%)",sep=""),
       ylab=paste("PC2 (",axis.percent(
         ord)[[2]],"%)",sep=""),
       main=title)
  #ordispider(ord, col='darkgrey',lty=3, scaling=3)
  ordiellipse(ord, group, col=group, lty=3, kind='sd', scaling=scaling.val)
  points(ord,dis='sp',pch=4,col='grey',cex=0.6,scaling=scaling.val)
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),],pch=4,col='blue',cex=1)	
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),],pch=4,col='blue',cex=1) 
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),],pch=4,col='blue',cex=1)
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),],pch=4,col='blue',cex=1)
  points(ord, dis='sites',pch=21,bg=group,cex=3,scaling=scaling.val)
  legend(legend.position,legend=unique(factor(group)),pch=19,col=unique(factor(group)))
}

# END


custom.plot.rda <- function(ord,group,title,plottype,sd.val,scaling.val,legend.position){
  plot(ord, type=plottype, scaling=scaling.val, 
       #xlim=c(-0.4,0.5), ylim=c(-0.4,0.4),
       xlab=paste("RDA1 (",axis.percent(
         ord)[[1]],"%)",sep=""),
       ylab=paste("RDA2 (",axis.percent(
         ord)[[2]],"%)",sep=""),
       main=title)
  #ordispider(ord, col=group,lty=3, scaling=scaling.val)
  ordiellipse(ord, group, col=group, lty=3, kind='sd', scaling=scaling.val)
  points(ord,dis='sp',pch=4,col='grey',cex=0.6,scaling=scaling.val)
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1])),],pch=4,col='blue',cex=1)	
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1]))),],pch=4,col='blue',cex=1) 
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2])),],pch=4,col='blue',cex=1)
  points(scores(ord,scaling=scaling.val)$sp[which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2]))),],pch=4,col='blue',cex=1) 
  #text(scores(ord,scaling=scaling.val)$sp[names(which(scores(ord)$sp[,1] > sd.val * sd(scores(ord)$sp[,1]))),],cex=0.8)
  #text(scores(ord,scaling=scaling.val)$sp[names(which(scores(ord)$sp[,1] < 0 - (sd.val * sd(scores(ord)$sp[,1])))),],cex=0.8) 
  #text(scores(ord,scaling=scaling.val)$sp[names(which(scores(ord)$sp[,2] > sd.val * sd(scores(ord)$sp[,2]))),],cex=0.8)
  #text(scores(ord,scaling=scaling.val)$sp[names(which(scores(ord)$sp[,2] < 0 - (sd.val * sd(scores(ord)$sp[,2])))),],cex=0.8) 
  points(ord, dis='sites',pch=21,bg=group,cex=3,scaling=scaling.val)
  legend(legend.position,legend=unique(factor(group)),pch=19,col=unique(factor(group)))
}

# END

#### Calculate the p values within a correlogram
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

