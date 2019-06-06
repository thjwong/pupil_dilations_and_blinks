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

Bsubjmin <- min(min(rowMeans(fdata[[1]],na.rm=TRUE,dims=1)),min(rowMeans(fdata[[2]],na.rm=TRUE,dims=1)))
Bsubjmax <- max(max(rowMeans(fdata[[1]],na.rm=TRUE,dims=1)),max(rowMeans(fdata[[2]],na.rm=TRUE,dims=1)))

plot(NA,NA,xlim=c(bin[1]/1000,bin[length(bin)]/1000),ylim=c(0,1),xlab="", ylab="",axes=FALSE)

##Compute p-values of ttest comparisons on each msec & then plot on minimal Y
Bttest <- rep(NA,length(bin))
##handle exception if all values are the same (equal) or are all missing
t.test.2sam.pv <- function(s1,s2) {
    obj<-try(t.test(s1,s2,paired=T), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
##Plot t.test results on the bottom
for(i in 1:length(bin)) {
    Bttest[i] <- t.test.2sam.pv((fdata[[1]])[i,],(fdata[[2]])[i,])
}

BttestC <- rep(NA,length(bin))
for(i in 1:length(bin)) {
    if(!is.na(Bttest[i])) {
        #if(Bttest[i] > 0.05) {BttestC[i]<-0}
        if(Bttest[i] > BonfCorr005) {BttestC[i]<-0}
        else {
            #if(Bttest[i] <= 0.05 & Bttest[i] > 0.01) {BttestC[i]<-1}
            if(Bttest[i] <= BonfCorr005 & Bttest[i] > BonfCorr001) {BttestC[i]<-1}
            #else {if(Bttest[i] <= 0.01) {BttestC[i]<-1}}
            else {if(Bttest[i] <= BonfCorr001) {BttestC[i]<-1}}
        }
    }
    else {BttestC[i]<-0}
}

BttestC2 <- rep(NA,length(bin))
brunlen<-rle(BttestC)
#bcon<-which(brunlen$lengths>=500 & brunlen$values==1)
bcon<-which(brunlen$lengths>=100 & brunlen$values==1)
if(length(bcon)>0) {
    ##fill in the head of series if not begin by long1
    if(bcon[1]!=1) {
        bsum <- sum(brunlen$lengths[1:(bcon[1]-1)])
        BttestC2[1:bsum]<-0
    }
    else {
        bsum <- brunlen$lengths[1]
        BttestC2[1:bsum]<-0
    }
    for(i in 1:length(bcon)) {
        #fill from the first segment
        BttestC2[(sum(brunlen$lengths[1:(bcon[i]-1)])+1):sum(brunlen$lengths[1:bcon[i]])]<-1
        if(i < length(bcon)) {
            BttestC2[(sum(brunlen$lengths[1:bcon[i]])+1):sum(brunlen$lengths[1:(bcon[i+1]-1)])]<-0
        } else {if (sum(brunlen$lengths[1:bcon[i]]) < length(bin)) { #end of list
            BttestC2[(sum(brunlen$lengths[1:bcon[i]])+1):sum(brunlen$lengths)]<-0
        }}
    }
} else {
    #print("No >500 significant segment")
    print("No >100 significant segment")
    BttestC2<-rep(0,length(bin))
}

for(i in 1:length(bin)) {
    if(BttestC2[i]==1) {points(bin[i]/1000,0.02,pch="|",col="pink",cex=1.5)}#cyan, magenta
}

BsubjBmeans <- rowMeans(fdata[[1]],na.rm=TRUE,dims=1)
BsubjMmeans <- rowMeans(fdata[[1]],na.rm=TRUE,dims=1)
lines(bin/1000,BsubjBmeans,col="black",type="l",lty=1,lwd=3)

BsubjBmeans <- rowMeans(fdata[[2]],na.rm=TRUE,dims=1)
BsubjMmeans <- rowMeans(fdata[[2]],na.rm=TRUE,dims=1)
lines(bin/1000,BsubjBmeans,col="darkgrey",type="l",lty=4,lwd=3)
if(FALSE){
    BsubjBsd <- sd(t(Bsubjmatrix),na.rm=TRUE)
    BsubjYlower <- BsubjBmeans-(1.96*BsubjBsd)
    BsubjYupper <- BsubjBmeans+(1.96*BsubjBsd)
    library("Hmisc")
    errbar(bin,BsubjBmeans,BsubjYupper,BsubjYlower,add=TRUE,col="dark red")
}
lines(x=c(0,0),y=c(0,1),col="black",lty=3)
axis(2,tick=T,pos=bin[1]/1000,las=2,cex.axis=1.5)
##X-axis
axis(1,tick=T,tck=-0.01,pos=0,cex.axis=1.5) ##X-axis at the lowest limit

##Check the continuous regions:
rleC2<-rle(BttestC2)
TRegion<-NULL
TEndtime<-NULL
for(i in 1:length(rle(BttestC2)[1]$lengths)) {
    TRegion <- c(TRegion,rle(BttestC2)[2]$values[i])
    TEndtime <- c(TEndtime,((bin[1]-1)+sum(rle(BttestC2)[1]$lengths[1:i])))
}

##Store the region info
BRegions<-data.frame(TRegion,TEndtime)
