logit<-function(x)log(x/(1-x))
convertOddsToProp<-function(z){z/(1+z)}

qnorm(10^-10, mean = 0, sd = 1, lower.tail = F)

stdErrBinom<-function(n1,n2,p1,p2) {
  (p1-p2)/sqrt(((p1*(1-p1))/n1)+((p2*(1-p2))/n2))
}

stdErrOR<-function(n1,n2,p1,p2) {
  sqrt(
    (1/(p1*n1))+
    (1/((1-p1)*n1))+
    (1/(p2*n2))+
    (1/((1-p2)*n2))
  )
}

getTwoWeekCombo<-function(subset1,subset2,filterAlpha=0.05){
  #addReference(subset1)
  qAlpha<-qnorm(filterAlpha/2, mean = 0, sd = 1, lower.tail = F)
  
  tempCombo<-merge((subset1),
                   (subset2),
                   by="posBase",all=T
  )
#  tempCombo<-imputeNonMatched(tempCombo)
  
  tempCombo$stderr<-sqrt((1/tempCombo$nBase.x)+
                           (1/(tempCombo$nSeq.x-tempCombo$nBase.x))+
                           (1/tempCombo$nBase.y)+
                           (1/(tempCombo$nSeq.y-tempCombo$nBase.y))
  )
  tempCombo$filtered
  
  tempCombo$lower<-(logit(tempCombo$pBase.y)-logit(tempCombo$pBase.x))-(qAlpha*tempCombo$stderr)
  tempCombo$upper<-(logit(tempCombo$pBase.y)-logit(tempCombo$pBase.x))+(qAlpha*tempCombo$stderr)
  tempCombo$incl<-(tempCombo$lower>0)|(tempCombo$upper<0)
  tempCombo[is.na(tempCombo$incl),"incl"]<-F
  return(tempCombo)
}

addLtoSNPfreq<-function(sfDF,deltaT=1){
  
  stdErrBinom<-function(n1,n2,p1,p2) {
    (p2-p1)/sqrt(((p1*(1-p1))/n1)+((p2*(1-p2))/n2))
  }
  
  sfDF<-sfDF%>%group_by(posBase)%>%
    arrange(weeksince)%>%
    mutate(
      #L1=logit(lead(propPosBase,n=deltaT)),
      #L2=logit(propPosBase),
      currentL=logit(lead(pBase,n=deltaT))-logit(pBase),
      nextPBase=lead(pBase,n=deltaT),
      nextNBase=lead(nBase,n=deltaT),
      nextNSeq=lead(nSeq,n=deltaT)#,
      #posDif=getciDF(nSeq,nextNSeq,propPosBase,nextPropPosBase)
    )%>%
    ungroup()%>%
    mutate(
      stdErrBinom=stdErrBinom(nSeq,nextNSeq,pBase,nextPBase),
      stdErrOR=stdErrOR(nSeq,nextNSeq,pBase,nextPBase)
    )%>%
    as.data.frame()
  
  return(sfDF)
}


spaceInvaderPlot<-function(snpFreqs,tFirst,tSecond,examplePosBase="39465T",filterAlpha=0.05){
  #allPropsDE[allPropsDE$weeksince==104&allPropsDE$pBase>0.1&allPropsDE$pBase<0.7,]
  snpFreqs$posBase<-paste0(snpFreqs$pos,snpFreqs$base)
  
  twoWeekDF<-getTwoWeekCombo(snpFreqs[snpFreqs$weeksince==tFirst,],snpFreqs[snpFreqs$weeksince==tSecond,],filterAlpha=filterAlpha)
  exampleSNPDF<-twoWeekDF[twoWeekDF$posBase==examplePosBase,]
  print(exampleSNPDF)
  ggplot(twoWeekDF)+
    geom_abline(slope=1,intercept = (c(c(1,2,3)*(5/7),(c(-1,-2,-3)*(5/7)))),linetype=2)+
    geom_abline(slope=1,intercept = 0,linetype=1)+
    geom_point(aes(x=logit(pBase.x),y=logit(pBase.y),color=incl))+
    geom_point(data=exampleSNPDF[exampleSNPDF$weeksince.x==tFirst,],aes(x=logit(pBase.x),y=logit(pBase.y)),
             size=5,shape=21,color="black",fill="#00000040" )+
    geom_point(data=exampleSNPDF[exampleSNPDF$weeksince.x==tFirst,],aes(x=logit(pBase.x),y=logit(pBase.y)),
             size=1,color="black")+
    geom_vline(xintercept = -log10(1),linetype=3)+
    geom_hline(yintercept = -log10(1),linetype=3)+
    scale_x_continuous(limits=c(-7,7))+
    scale_y_continuous(limits=c(-7,7))+
    scale_color_manual(limits=c(T,F),values=c("#0000FF","#FF000020"))+
    theme_bw()+
    theme(legend.position = "none")+
    labs(title="",x="logit(SNP proportion previous week)",y="logit(SNP proportion current week)")+
    coord_fixed()
}


snpProportionTimeline<-function(snpFreqs,tFirst,tLast,tExample,examplePosBase="39465T"){
  snpFreqs$posBase<-paste0(snpFreqs$pos,snpFreqs$base)
  exampleSNPDF<-as.data.frame(snpFreqs[snpFreqs$posBase==examplePosBase,])
  print(exampleSNPDF)
  Lexample<-as.numeric(logit(exampleSNPDF[exampleSNPDF$weeksince==(tExample+1),"pBase"])-
    logit(exampleSNPDF[exampleSNPDF$weeksince==tExample,"pBase"]))
  print(Lexample)
  projectionDF=data.frame(
    t=(tFirst:tLast),
    projectedProp=convertOddsToProp(exp(logit(exampleSNPDF[exampleSNPDF$weeksince==tExample,"pBase"])+(((tFirst:tLast)-tExample)*Lexample)))
    )
  print(projectionDF)
  
  ggplot(snpFreqs[(snpFreqs$weeksince>=tFirst)&(snpFreqs$weeksince<=tLast),])+
  #ggplot(addReference(snpFreqs[snpFreqs$posBase%in%twoWeekDF[twoWeekDF$incl,"posBase"],]))+
    geom_line(aes(x=fConvWtoDate(weeksince),y=pBase,group=posBase),color="#00000030")+
    
      geom_line(data=projectionDF,
        aes(x=fConvWtoDate(t),y=projectedProp),color="#0000FFA0",size=2
      )+  
      geom_point(data=exampleSNPDF[exampleSNPDF$weeksince==tExample,],aes(x=fConvWtoDate(weeksince),y=pBase),size=5,shape=21,
                color="black",fill="#0000FF70")+
    scale_x_date(limits=c(fConvWtoDate(tFirst),fConvWtoDate(tLast)),
               minor_breaks = as.Date("2020-12-01")+months(1:24),#breaks=c(),
               #sec.axis = dup_axis()
    )+
    labs(x="Collection date",y="SNP Proportion")+
    theme_bw()
}


L_overTimePlot<-function(snpFreqs,tFirst,tLast,tExample,filterAlpha){

  snpFreqs$posBase<-paste0(snpFreqs$pos,snpFreqs$base)
  snpFreqs<-addLtoSNPfreq(snpFreqs)
  tempSNPfreqs<-snpFreqs
  qAlpha<-qnorm(filterAlpha, mean = 0, sd = 1, lower.tail = F)
  
  tempSNPfreqs$lower<-(tempSNPfreqs$currentL)-(qAlpha*tempSNPfreqs$stdErrOR)
  tempSNPfreqs$upper<-(tempSNPfreqs$currentL)+(qAlpha*tempSNPfreqs$stdErrOR)
  tempSNPfreqs<-tempSNPfreqs[(!is.na(tempSNPfreqs$lower)|is.infinite(tempSNPfreqs$lower))&
                               (!is.na(tempSNPfreqs$upper)|is.infinite(tempSNPfreqs$upper))
                             ,]
  tempSNPfreqs<-tempSNPfreqs[tempSNPfreqs$lower>0|tempSNPfreqs$upper<0,]
  
  timerange<-min(tempSNPfreqs$weeksince):max(tempSNPfreqs$weeksince)
  tempSNPfreqs<-merge(
    tempSNPfreqs,
    bind_rows(lapply(unique(tempSNPfreqs$posBase),function(x)data.frame(posBase=x,dummy="-",weeksince=timerange))),
    by=c("posBase","weeksince"),all=T
  )
  
  ggplot(tempSNPfreqs)+
    geom_point(data=snpFreqs,aes(x=fConvWtoDate(weeksince),y=currentL),color="#FF000001")+
    geom_point(aes(x=fConvWtoDate(weeksince),y=currentL),color="#0000FF20")+
    geom_line(aes(x=fConvWtoDate(weeksince),y=currentL,group=posBase),color="#0000FF20")+
    scale_y_continuous(limits=c(-3,3))+
    scale_x_date(limits=c(as.Date("2021-01-01"),as.Date("2022-11-01")))+
    labs(x="Collection date",y="Logit difference in SNP proportion (L)")+
    theme_bw()
}
