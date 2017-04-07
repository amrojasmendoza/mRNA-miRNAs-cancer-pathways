# Analysis of correlations between gene expression and miRNA expressio
# by Sergio Alonso
# Jan 2017

sink(file="log.txt",split=T)

# Load required libraries. Install them if necessary  ----

if(!require(corrplot)) {install.packages("corrplot");library(corrplot)}
if(!require(survival)) {install.packages("survival");library(survival)}

# Create output directories, if needed ----

if(!file.exists("figures")) dir.create("figures")
if(!file.exists("tables")) dir.create("tables")

# Auxiliary functions ----

# Substitute infinite values by either NAs or the minimum value != -Infinite
subInf <- function(x,substitute=c("na","min")) {
  x[is.infinite(x)] <- NA
  if(match.arg(substitute)=="min") x[is.na(x)] <- min(x,na.rm=T)
  return(x)
}

# Apply padjust to a matrix of Pvalues
# It can be applied considering all the values in the matrix (all)
# or by row (row)
# or by column (col)

p.adjust2 <- function(m,by=c("all","row","col"),method="fdr") {
  p <- m
  by <- match.arg(by,c("all","row","col"))
  
  if(by=="all") {
    p <- matrix(p.adjust(unlist(p),method=method),ncol=ncol(p))
    colnames(p) <- colnames(m)
    rownames(p) <- rownames(m)
  }
  if(by=="col") {
    for(i in 1:ncol(p)) {
      p[,i] <- p.adjust(p[,i],method=method)
    }
  }
  if(by=="row") {
    for(i in 1:nrow(p)) {
      p[i,] <- p.adjust(p[i,],method=method)
    }
  }
  return(p)
}

# Functions to convert ß-values and M-values (methylation)

BtoM <- function(beta) log(beta/(1-beta))
MtoB <- function(m) exp(m)/(1+exp(m))

# Lineal regression analysis functions ----

# Create a lineal model with gene expression as dependent variable and Tumor
# type (if applicable), CNAs, Methylation and miRNA expression as explanatory
# variables. CNAs can enter the function as categorical variables (as a factor)
# or as a continuos numeric variable based on the GISTIC score
# see (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3218867/)

generateLM <- function(gene,mirna,dataset,TumorAsFactor=T,cna=c("Categorical","Continous")) {
  gene <- as.character(gene)
  mirna <- gsub("-",".",mirna,fixed = T)
  
  geneExp <- subInf(log(dataset[,gene]))
  mirnaExp <- subInf(log(dataset[,mirna]))
  
  cna <- match.arg(cna)
  
  if(cna=="Categorical") CNA <- dataset[,paste0(gene,".CNA")]
  if(cna=="Continous") CNA <- as.numeric(dataset[,paste0(gene,".CNA")])
  
  metprobes <- grep(paste0(gene,".cg"),colnames(dataset),value = T)
  MET <- data.frame(dataset[,metprobes])
  # Remove methylation probes with data for less than 15 patients
  
  valid <- apply(!is.na(MET),2,sum) >= 15
  MET <- MET[,valid]
  metprobes <- colnames(MET) <- gsub(gene,"MET",colnames(MET))
  
  d0 <- data.frame(Tumor=factor(dataset$Tumor),geneExp=geneExp,mirnaExp=mirnaExp,CNA=CNA)
  d0 <- cbind(d0,MET)
  
  # To avoid collapsing models because 
  # the number of factors is larger than the number
  # of valid observations, the models will include a 
  # a subsset of best-correlating methylation probes
  
  n <- nrow(na.omit(d0)) # complete observations
  
  if(n < 5) return(NA)
  
  n <- n-5-1 # maximum number of meth probes in the model ( degress of freedom in CNA and 1 in mirnaExp)
  
  if(length(metprobes)>n) {
    met <- coef(summary(lm(d0$geneExp ~ .,MET)))[-1,] 
    met <- met[order(met[,4]),] # sort by p-value
    metprobes <- rownames(met)[1:n]
  }
  
  # Construct the appropriate formula
  
  metfactors <- paste(c("",metprobes),collapse=" + ")
  
  if(length(levels(d0$Tumor))<2) TumorAsFactor <- F
  
  formula0 <- ifelse(TumorAsFactor,
                     "geneExp ~ Tumor + CNA %s + mirnaExp + Tumor:mirnaExp",
                     "geneExp ~ CNA %s + mirnaExp")
  
  formula0 <- formula(sprintf(formula0,metfactors))
  
  lm0 <- NA
  try(lm0 <- lm(formula0,d0),silent = F)
  return(lm0)
  
}

# Calculate an ANOVA table from a lineal model, calculating also correlation
# coeffiicnets for every independent variable. Methylation probes are collapsed
# into a single MET value.

anovaLM <- function(lm0) {
  a <- NA
  if(!is.na(lm0)[1]) {
    a <- anova(lm0)
    colnames(a)[5] <- "pVal"
    
    # Combine the MET probes
    metprobes <- which(substr(rownames(a),1,3)=="MET")
    m <- colSums(a[metprobes,],na.rm=T)
    
    m["Mean Sq"] <- m["Sum Sq"]/m["Df"]
    m["F value"] <- m["Mean Sq"]/a["Residuals","Mean Sq"]
    m["pVal"] <- pf(m["F value"],m["Df"],a["Residuals","Df"],lower.tail = F)
    if(m["Df"]==0) m <- rep(NA,5)
    
    variable <- grep("MET",rownames(a),value = T,invert = T)
    
    cna <- which(variable=="CNA")
    variable <- append(variable,"MET",after = cna)
    a <- a[variable,]
    rownames(a) <- variable
    a["MET",] <- m
    
    a$r2 <- a$`Sum Sq`/sum(a$`Sum Sq`,na.rm=T)
    
    a$r <- sqrt(a$r2) * sign(coef(lm0)[c(rownames(a))])
  }
  return(a)
}

# Wrapper function to generate the lineal models, the anova tables and populate
# the correlation and Pvalues matrices

generateLModels <- function(env) {
  with(env,{
    
    models <- list()
    anovas <- list()
    
    correlations <- matrix(NA,ncol=length(cancers)+1,nrow=nrow(pairs))
    colnames(correlations) <- c(cancers,"ALL")
    rownames(correlations) <- sprintf("%10s - %16s",pairs$Gen,pairs$miRNA)
    pVals <- correlations
    
    for(cancer in c(cancers,"ALL")) {
      for(pair in 1:nrow(pairs)) {
        pairname <- sprintf("Pair%02i",pair)
        
        gene <- as.character(pairs$Gen[pair])
        mirna <- as.character(pairs$miRNA[pair])
        
        if(cancer=="ALL") {
          m <- tryCatch(generateLM(gene,mirna,rpkm,TumorAsFactor = T),finally=NA)
        } else {
          m <- tryCatch(generateLM(gene,mirna,subset(rpkm,Tumor==cancer),TumorAsFactor = F),finally=NA) 
        }
        
        a <- tryCatch(anovaLM(m),finally=NA)
        
        cat("\n\nLineal regression analysis",cancer,gene,mirna,fill=T)
        print(a)
        
        try({
          correlations[pair,cancer] <- a["mirnaExp","r"]
          pVals[pair,cancer] <- a["mirnaExp","pVal"]
        },silent = T)
        
        
        models[[cancer]][[pairname]] <- m
        anovas[[cancer]][[pairname]] <- a
      }
    }
    
    pValsAdj <- p.adjust2(pVals,by = "row")
    
    select <- correlations[,"ALL"] < 0 & pValsAdj[,"ALL"] < 0.05
    
    rm(cancer,pair,pairname,gene,mirna,m,a)
    
  })
}

# Coxph regression analysis functions ---- Create a coxph model including Stage
# and Gene expression as explanatory variables. Gene expression can be
# introduced as a categorical semiquantitative variable with nquantile>2
# categories, or as a continous variable if nquantiles is set to 0

coxphStageGene <- function(cancer,gene,mirna,dataset,nquantiles=4,parameters,substituteInfinite=c("min","na")) {
  
  substituteInfinite <- match.arg(substituteInfinite)
  
  qcut <- function(x,groups=4,verbose=F) {
    cuts <- seq(0,1,l=groups+1)
    if(verbose) cat(cuts)
    try(x <- as.numeric(cut(x,quantile(x,cuts,na.rm=T),include.lowest = T)))
    return(x)
  }
  
  gene <- as.character(gene)
  mirna <- as.character(mirna)
  d1 <- subset(dataset,Tumor==cancer)[,c("stage2",gene,mirna,"vital_status","time")]
  
  d1$geneExp <- NA
  d1$mirnaExp <- NA
  
  if(nquantiles>=1) {
    d1$geneExp <- qcut(d1[,gene],nquantiles)
    d1$mirnaExp <- qcut(d1[,mirna],nquantiles)
  }
  if(nquantiles<1) {
    d1$geneExp <- subInf(log(d1[,gene]),substitute = substituteInfinite)
    d1$mirnaExp <- subInf(log(d1[,gene]),substitute = substituteInfinite)
  }
  
  s1 <- Surv(d1$time,d1$vital_status=="Dead")
  f1 <- formula(sprintf("s1 ~ %s",parameters))
  
  print(c1 <- coxph(f1,d1))
  
  return(coef(summary(c1)))
} 



# Graphic functions and variables ----

emptyPlot <- function(xlim=c(-1,1),ylim=c(-1,1)) {
  plot(NA,xlim=xlim,ylim=ylim,xaxt="n",yaxt="n",
       xlab=NA,ylab=NA,bty="n")
}

# cMat: function to return a color matrix represented the values of the matrix m

cMat <- function(m,palette=colorRamp(c("green","black","red")),na.color="yellow",min=NULL,max=NULL) {
  
  if(!is.null(min)) m[m<min] <- min
  if(!is.null(max)) m[m>max] <- max
  
  # scale between 0 and 1
  max <- max(abs(m),na.rm=T)
  m <- m/max/2+.5
  nas <- is.na(m)
  m[nas] <- 0
  for(i in 1:nrow(m)) m[i,] <- rgb(palette(m[i,]),maxColorValue = 255)
  m[nas] <- na.color
  return(m)
}

# defaultCorrPlot: Wrapper function with the default parameters for corrplot

defaultCorrPlot <- function(cor,pval,cl.ratio=.2) {
  corrplot(cor,mar=c(4,1,4,1),
           addgrid.col = NA,
           na.label = "_",
           p.mat = pval,
           pch.cex = 1,
           pch=4,
           sig.level=0.05,
           tl.col="black",
           cl.align.text="l",
           cl.offset=.1,
           cl.ratio = cl.ratio,
           add=F)
  
  grid2(ncol(cor),nrow(cor),col="grey",lwd=2)
}

# GoPlot: function to plot a matrix of values as dots. 
# xcoordinates and ycoordinates can be explicity specified.

GoPlot <- function(size,color,pch=19,k=1,minSize=.2,cex.text=1,add=F,xcoord=NULL,ycoord=NULL,rownames=NULL,gridcolor="lightgrey") {
  if(is.null(rownames)) rownames <- rownames(size)
  x <- ncol(size)
  y <- nrow(size)
  if(is.null(xcoord)) xcoord=1:x
  if(is.null(ycoord)) ycoord=1:y
  xscale <- diff(range(xcoord)) / x
  if(!add) plot(NA,xlim=range(xcoord)+c(-xscale*10,+xscale),ylim=rev(range(ycoord))+c(2,-2),xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n")
  for(col in xcoord) segments(col,min(ycoord),col,max(ycoord),col=gridcolor,lwd=2)
  for(row in ycoord) segments(min(xcoord),row,max(xcoord),row,col=gridcolor,lwd=2)
  for(row in 1:y) {
    points(xcoord,rep(ycoord[row],x),cex=size[row,]*k+minSize,pch=pch,col=as.character(color[row,]))
  }
  #text(xcoord[1]/2,ycoord,rownames,pos=2,cex=cex.text) 
}

# grid2: plots a grid of nx vertical lines and ny horizontal lines

grid2 <- function(nx,ny,...) {
  segments((0:nx)+.5,0.5,(0:nx)+.5,ny+.5,...)
  segments(0.5,(0:ny)+.5,nx+.5,(0:ny)+.5,...)
}


# correlationPalette: color palette copied from library corrplot

correlationPalette <- c("#67001F", "#B2182B", "#D6604D", 
                        "#F4A582", "#FDDBC7", "#FFFFFF", 
                        "#D1E5F0", "#92C5DE", "#4393C3", 
                        "#2166AC", "#053061") 


# plotCorrSurvival: Complete function to plot the correlations and the coxph 
# coefficients in a single plot, using two overlapping GoPlots

plotCorrSurvival <- function(env1) {
  
  with(env1, {

    x1 <- as.matrix(correlations[select,cancers])
    x1 <- cMat(x1,colorRamp(correlationPalette),min=-1,max=+1)
    
    x2 <- as.matrix(HR.coef[select,cancers])
    x2 <- cMat(x2,colorRamp(rev(correlationPalette)),min=-log(2),max=log(2))
    
    s1 <- as.matrix(pVals[,cancers][select,])
    if(adjustFDR) s1 <- p.adjust2(s1,by="row")
    s1[s1>.05] <- NA
    s1[s1<0.001] <- 1.5
    s1[s1<0.01] <- 1
    s1[s1<0.05] <- .5
    
    s2 <- as.matrix(HR.pval[select,cancers])
    if(adjustFDR) s2 <- p.adjust2(s2,by="row")
    s2[s2>.05] <- NA
    s2[s2<0.001] <- 1.5
    s2[s2<0.01] <- 1
    s2[s2<0.05] <- .5
    
    bothSig <- !(is.na(s1) | is.na(s2))
    
    for(i in 1:nrow(s1)) {
      rect(1:ncol(s1)-.4,
           i-.4,
           1:ncol(s1)+.4,
           i+.4,
           col="white",
           lwd=ifelse(bothSig[i,],2,.5))}
    
    GoPlot(s1,x1,pch=16,
           xcoord=1:ncol(s1)-.15,
           add=T,gridcolor = NA)
    GoPlot(s2,x2,pch=18,
           xcoord=1:ncol(s2)+.15,
           add=T,gridcolor = NA)
    
    text(1:ncol(x1),.4,cancers,srt=90,adj=0,cex=1.5)
    text(.4,1:nrow(x1),rownames(x1),adj=1)
    
    
  })
}


# legend2: Legend for plotCorrSurvival

legend2 <- function(x0,y0,k=.2) {
  points(rep(x0-.5,4),y0+c(0,3:5),pch=16,cex=c(2.5,1.5,1,.5)+k)
  points(rep(x0+.5,4),y0+c(0,3:5),pch=18,cex=c(2.5,1.5,1,.5)+k)
  rect(x0-1,y0-1,x0+1,y0+1)
  text(x0+c(-.5,+.5),y0-1,c(" miRNA-Gene Corr"," Gene-Survival"),srt=90,adj=0)
  
  text(x0-1,y0+3:5,sprintf("P < %1.2g",c(.001,.01,.05)),adj=1)
  
  y1 <- y0+7+1:21
  
  rect(x0-.3,y1-.5,x0+.3,y1+.5,
       col=colorRampPalette(rev(correlationPalette))(21),
       border=NA)
  rect(x0-.3,y0+7.5,x0+.3,y0+28.5)
  text(x0-.4,y0+c(8,18,28),c(1,0,-1),adj=1)
  text(x0+.4,y0+c(8,18,28),c(-1,0,1),adj=0)
  text(x0-.4,y0+7,"Correlation",adj=1)
  text(x0+.4,y0+7,"log2(HR)",adj=0)
  
  
}




# Load data ----
# data is loaded into separate environments

general <- new.env()
with(general,{
  rpkm <- read.csv("csvdata/general.pairs.data.csv")
  pairs <- read.csv("csvdata/general.pairs.csv")
  cancers <- levels(rpkm$Tumor)
})


lung <- new.env()
with(lung,{
  rpkm <- read.csv("csvdata/lung.pairs.data.csv")
  pairs <- read.csv("csvdata/lung.pairs.csv")
  cancers <- levels(rpkm$Tumor)
})


# Correlation analysis ------------

generateLModels(general)
generateLModels(lung)

with(general,{
    pdf("figures/generalPairsCorrelations01.pdf",15,10)
    par(family="mono",cex=1)
    defaultCorrPlot(correlations,pValsAdj)
    title("Gene-miRNA correlation\nafter correcting for CNAs and DNA methylation")
    dev.off()
    
    pdf("figures/generalPairsCorrelations02.pdf",15,10)
    par(family="mono",cex=1)
    defaultCorrPlot(correlations[select,],pValsAdj[select,])
    title("Gene-miRNA correlation\nafter correcting for CNAs and DNA methylation")
    dev.off()
})
    
with(lung,{
  pdf("figures/lungPairsCorrelations01.pdf",10,15)
  par(family="mono",cex=1)
  defaultCorrPlot(correlations,pValsAdj,cl.ratio = .8)
  title("Lung-exclusive Gene-miRNA correlation\nafter correcting for CNAs and DNA methylation")
  dev.off()
  
  pdf("figures/lungPairsCorrelations02.pdf",10,10)
  par(family="mono",cex=1)
  defaultCorrPlot(correlations[select,],pValsAdj[select,],cl.ratio = .8)
  title("Lung-exclusive Gene-miRNA correlation\nafter correcting for CNAs and DNA methylation")
  dev.off()
})  

with(general,{
  sapply(anovas,function(cancer){
    sapply(cancer, function(pair) {
      x <- rep(NA,4)
      try(x <- pair[c("Tumor","CNA","MET","mirnaExp"),"r2"],silent=T)
      names(x) <- c("Tumor","CNA","MET","miRNA")
      return(x)
    })
  }) -> contributions
  
  contributions <- round(contributions*100,2)
  rownames(contributions) <- paste(rep(pairs$Gen,each=4),rep(pairs$miRNA,each=4),c("Tumor type","CNA","MET","miRNA"),sep=" - ")
})


# Write correlation analyses results tables 

write.csv(general$correlations,quote=F,file="tables/general.correlation.csv")
write.csv(general$pVals,quote=F,file="tables/general.correlation.pval.csv")
write.csv(general$contributions,quote=F,file="tables/general.correlation.contribution.csv")

write.csv(lung$correlations,quote=F,file="tables/lung.correlation.csv")
write.csv(lung$pVals,quote=F,file="tables/lung.correlation.pval.csv")


# Survival analysis ----

with(general,{
  sapply(cancers,function(cancer) {
    sapply(1:nrow(pairs),function(pair) {
      gene <- as.character(pairs$Gen[pair])
      mirna <- as.character(pairs$miRNA[pair])
      cat("\nSurvival analysis",cancer,gene,mirna,"\n")
      if(length(unique(subset(rpkm,Tumor==cancer)[,gene]))==1) return(c(NA,NA))
      tryCatch(coxphStageGene(cancer,gene,mirna,rpkm,nquantiles = 0,"stage2 + geneExp",substituteInfinite = "min")[2,c(1,5)],finally=c(NA,NA))
    })
  }) -> HR.mat
  
  HR.coef <- HR.mat[seq(1,nrow(HR.mat),2),]
  HR.pval <- HR.mat[seq(1,nrow(HR.mat),2)+1,]
  row.names(HR.coef) <- row.names(HR.pval) <- sprintf("%10s - %16s",pairs$Gen,pairs$miRNA)
})

with(lung,{
  sapply(cancers,function(cancer) {
    sapply(1:nrow(pairs),function(pair) {
      gene <- as.character(pairs$Gen[pair])
      mirna <- as.character(pairs$miRNA[pair])
      cat("\nSurvival analysis",cancer,gene,mirna,"\n")
      if(length(unique(subset(rpkm,Tumor==cancer)[,gene]))==1) return(c(NA,NA))
      tryCatch(coxphStageGene(cancer,gene,mirna,rpkm,nquantiles = 0,"stage2 + geneExp",substituteInfinite = "min")[2,c(1,5)],finally=c(NA,NA))    })
  }) -> HR.mat
  
  HR.coef <- HR.mat[seq(1,nrow(HR.mat),2),]
  HR.pval <- HR.mat[seq(1,nrow(HR.mat),2)+1,]
  row.names(HR.coef) <- row.names(HR.pval) <- sprintf("%10s - %16s",pairs$Gen,pairs$miRNA)
  
  })


general$adjustFDR=T
pdf("figures/generalPairsSurvival.pdf",12,12)
par(family="mono")
emptyPlot(xlim=c(-12,16),ylim=c(36,-5))
plotCorrSurvival(general)
legend2(-9.5,7,.3)
dev.off()

lung$adjustFDR=T
pdf("figures/lungPairsSurvival.pdf",12,13)
par(family="mono")
emptyPlot(xlim=c(-12,16),ylim=c(40,-5))
plotCorrSurvival(lung)
legend2(-9.5,9,.3)
dev.off()

write.csv(general$HR.coef,quote=F,file = "tables/general.coxph.coefs.csv")
write.csv(general$HR.pval,quote=F,file = "tables/general.coxph.pvals.csv")
write.csv(lung$HR.coef,quote=F,file = "tables/lung.coxph.coefs.csv")
write.csv(lung$HR.pval,quote=F,file = "tables/lung.coxph.pvals.csv")

sink()








