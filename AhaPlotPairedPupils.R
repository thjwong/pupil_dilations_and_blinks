## Bonferroni Correction
## http://www.brainvoyager.com/bvqx/doc/UsersGuide/StatisticalAnalysis/TheMultipleComparisonsProblem.html
## Adjusted p-value thresholds = threshold / N, so that:
BonfCorr005 <- 0.05 / length(bin) #20001
BonfCorr001 <- 0.01 / length(bin) #20001
    
## Use this to find out the critical clsuter size
## "the [α × N]+1 th largest value over the sampling distribution"
## https://stats.stackexchange.com/questions/196769/multiple-comparison-correction-for-temporally-correlated-tests
## There simply were not as many clusters:
## > sort(brunlen$lengths[which(brunlen$values==1)])
##  [1]   1   1   1   2   2  14  43  76 175 357
## > sort(prunlen$lengths[which(prunlen$values==1)])
## [1]   42 2035

#source("AhaPupil.R")

##***Assume there are only 2 list in global environment that both end with .lst***
##pnames <- ls(pattern=".lst")
##pdata <- list()
##for(i in pnames) pdata[[i]] <- get(i)

##create 2 matrix lists to store pupil & pupil maxmin (peak & minimum) matrices of two groups of subjects
mdata <- list()
mmdata <- list()

for(i in 1:2) {
    tmp.lst <- comp.lst[[i]] ##*Assume comp.lst already stored the two list for comparison!
    print(paste("Processing comp.lst[[",i,"]]",sep=""),quote=FALSE)
    subj <- unique(g[which(tmp.lst),,drop=FALSE]$Subject)
    Psubjmatrix <- NULL
    Psubjmaxmin <- NULL
    session <- g[which(tmp.lst),,drop=FALSE]$Session[1]
    for(j in 1:length(subj)) {
        subjindex <- which(tmp.lst & g$Subject==subj[j])
        print(paste("Start working on Subject ",subj[j],sep=""),quote=FALSE)
        ##Compute trial matrix and then add rowmean to subject pupil matrix
        Ptrialmatrix <- Ppmastermatrix2[subjindex,,drop=FALSE] ##HERE!!!
        Ptrialmatrix <- t(Ptrialmatrix)
        Psubjmatrix <- cbind(Psubjmatrix, rowMeans(Ptrialmatrix, na.rm=TRUE, dims=1))
        ##Compute max and min pupilsize for each trial and then add to subject max matrix
        Pmaxmin <- NULL
        for(k in 1:ncol(Ptrialmatrix)) {
            Pmaxmin <-cbind(Pmaxmin,c(max(Ptrialmatrix[,k],na.rm=TRUE),min(Ptrialmatrix[,k],na.rm=TRUE)))
        }
        Psubjmaxmin <- cbind(Psubjmaxmin,rowMeans(Pmaxmin, na.rm=TRUE, dims=1))
    }
    mdata[[i]] <- Psubjmatrix
    ##or mdata[[names(pdata)[i]]] <- Psubjmatrix
    mmdata[[i]] <- Psubjmaxmin
}

Psubjmin <- min(min(rowMeans(mdata[[1]],na.rm=TRUE,dims=1)),min(rowMeans(mdata[[2]],na.rm=TRUE,dims=1)))
Psubjmax <- max(max(rowMeans(mdata[[1]],na.rm=TRUE,dims=1)),max(rowMeans(mdata[[2]],na.rm=TRUE,dims=1)))

##Define plot area
plot(NA,NA,xlim=c(bin[1],bin[length(bin)]),ylim=c(Psubjmin-.05,Psubjmax+.05),xlab="",ylab="",axes=FALSE)

##Compute p-values of ttest comparisons on each msec & then plot on minimal Y
Pttest <- rep(NA,length(bin))
##handle exception if all values are the same (equal) or are all missing
t.test.2sam.pv <- function(s1,s2) {
    obj<-try(t.test(s1,s2,paired=T), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
##Plot t.test results on the bottom
for(i in 1:length(bin)) {
    Pttest[i] <- t.test.2sam.pv((mdata[[1]])[i,],(mdata[[2]])[i,])
}

PttestC <- rep(NA,length(bin))
for(i in 1:length(bin)) {
    if(!is.na(Pttest[i])) {
        #if(Pttest[i] > 0.05) {PttestC[i]<-0}
        if(Pttest[i] > BonfCorr005) {PttestC[i]<-0}
        else {
            #if(Pttest[i] <= 0.05 & Pttest[i] > 0.01) {PttestC[i]<-1}
            if(Pttest[i] <= BonfCorr005 & Pttest[i] > BonfCorr001) {PttestC[i]<-1}
            #else {if(Pttest[i] <= 0.01) {PttestC[i]<-1}}
            else {if(Pttest[i] <= BonfCorr001) {PttestC[i]<-1}}
        }
    }
    else {PttestC[i]<-0}
}

PttestC2 <- rep(NA,length(bin))
prunlen<-rle(PttestC)
#pcon<-which(prunlen$lengths>=500 & prunlen$values==1)
pcon<-which(prunlen$lengths>=100 & prunlen$values==1)
#pcon<-which(prunlen$lengths>=50 & prunlen$values==1)
if(length(pcon)>0) {
    ##fill in the head of series if not begin by long1
    if(pcon[1]!=1) {
        psum <- sum(prunlen$lengths[1:(pcon[1]-1)])
        PttestC2[1:psum]<-0
    }
    for(i in 1:length(pcon)) {
        #fill from the first segment
        PttestC2[(sum(prunlen$lengths[1:(pcon[i]-1)])+1):sum(prunlen$lengths[1:pcon[i]])]<-1
        if(i < length(pcon)) {
            PttestC2[(sum(prunlen$lengths[1:pcon[i]])+1):sum(prunlen$lengths[1:(pcon[i+1]-1)])]<-0
        } else {if (sum(prunlen$lengths[1:pcon[i]]) < length(bin)) {
            PttestC2[(sum(prunlen$lengths[1:pcon[i]])+1):sum(prunlen$lengths)]<-0
        }}
    }
} else {
    #print("No >500 significant segment")
    print("No >100 significant segment")
    #print("No >50 significant segment")
    PttestC2<-rep(0,length(bin))
}

for(i in 1:length(bin)) {
    if(PttestC2[i]==1) {points(bin[i],Psubjmin,pch="|",col="pink",cex=1.5)}#cyan, magenta
}

PsubjPmeans <- rowMeans(mdata[[1]],na.rm=TRUE,dims=1)
PsubjMmeans <- rowMeans(mdata[[1]],na.rm=TRUE,dims=1)
lines(bin,PsubjPmeans,col="black",type="l",lty=1,lwd=6)

PsubjPmeans <- rowMeans(mdata[[2]],na.rm=TRUE,dims=1)
PsubjMmeans <- rowMeans(mdata[[2]],na.rm=TRUE,dims=1)
lines(bin,PsubjPmeans,col="darkgrey",type="l",lty=4,lwd=6)
if(FALSE){
    PsubjPsd <- sd(t(Psubjmatrix),na.rm=TRUE)
    PsubjYlower <- PsubjPmeans-(1.96*PsubjPsd)
    PsubjYupper <- PsubjPmeans+(1.96*PsubjPsd)
    library("Hmisc")
    errbar(bin,PsubjPmeans,PsubjYupper,PsubjYlower,add=TRUE,col="dark red")
}
lines(x=c(0,0),y=c(Psubjmin-0.01,Psubjmax+0.01),col="black",lty=3)
axis(2,tick=T,pos=bin[1],las=2)
##X-axis
axis(1,tick=T,tck=-0.01,pos=Psubjmin-0.01)
