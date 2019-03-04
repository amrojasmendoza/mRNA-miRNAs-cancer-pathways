# this script must be run after "analysis.R"
# questions about this script to Sergio Alonso (2019). Credit to salonsou@igtp.cat


par(mfrow=c(3,5))
for(i in cancers) {

# extract the r calculated from the ANOVA analyses
  
cors <- sapply(general$anovas[[i]],function(x) {
  if(is.data.frame(x)) x["mirnaExp","r"] else NA 
})

# calculate partial correlations

pcors <- sapply(general$models[[i]],function(x) {
  if(is.list(x)) {
    geneColumn <- which(colnames(x$model)=="geneExp")
    mirnaColumn <- which(colnames(x$model)=="mirnaExp")
  
    res1 <- lm(geneExp ~ .,x$model[,-mirnaColumn])$residuals
    res2 <- lm(mirnaExp ~ .,x$model[,-geneColumn])$residuals
    cor(res1,res2)} else NA
})

# plot r vs partial correlation

range0 <- range(c(cors,pcors),na.rm=T)

plot(pcors,cors,
     xlab="Partial correlations",
     ylab="signed sqrt(r2)",
     xlim=range0,
     ylim=range0,
     pch=19)
abline(0,1)
abline(h=0,v=0,lty=2)
title(i)
mtext(sprintf("r2=%1.2f",cor(pcors,cors,use="complete")^2),
      3,-2,cex=.9)

}
