
calculateAllProportions<-function(dataDF,nSeqSet=NULL){
  localProps<-dataDF%>%
    filter(!is.na(weeksince))%>%
    group_by(weeksince,pos,base)%>%
    summarise(
      nBase=n(),
      nSeq=nSeq[1]
    )%>%
   # mutate(nSeq=sum(nBase))%>%
    ungroup%>%
    add_row(
      group_by(.,pos,base)%>%
      reframe(
            weeksince=min(weeksince):max(weeksince),
            nBase=0,
            nSeq=0
      )
    )%>%ungroup()%>%
    group_by(pos,base,weeksince)%>%
    summarise(
      nBase=sum(nBase),
      nSeq=sum(nSeq)
    )%>%ungroup()%>%
    group_by(weeksince)%>%
    mutate(
      nSeq=max(nSeq,na.rm = T),
      pBase=nBase/nSeq
    )%>%
    ungroup()%>%
    add_row(#Add reference positions
      group_by(.,pos,weeksince)%>%
        summarise(
          base="Reference",#sREFch[pos],
          nBase=nSeq[1]-sum(nBase),
          pBase=nBase/nSeq[1],
          nSeq=nSeq[1])
    )%>%
  #  arrange(base,pos,weeksince)%>% #TODO: delete the arrange and grouping. 
    ungroup()%>%
   # group_by(pos,base)%>%
    mutate(
      oneMinuspB=1-pBase,
      Hx=pBase*oneMinuspB
    )%>%ungroup()
  return(localProps)
}

calculateAllProportionsNSeq<-function(dataDF,nSeqSet=10){
  localProps<-dataDF%>%
    filter(!is.na(weeksince))%>%
    group_by(weeksince,pos,base)%>%
    summarise(
      nBase=n(),
      nSeq=nSeqSet
    )%>%
    mutate(nSeq=sum(nBase))%>%
    ungroup%>%
    add_row(
      group_by(.,pos,base)%>%
        reframe(
          weeksince=min(weeksince):max(weeksince),
          nBase=0,
          nSeq=0
        )
    )%>%ungroup()%>%
    group_by(pos,base,weeksince)%>%
    summarise(
      nBase=sum(nBase),
      nSeq=(nSeqSet)
    )%>%ungroup()%>%
    group_by(weeksince)%>%
    mutate(
      nSeq=nSeqSet,#max(nSeq,na.rm = T),
      pBase=nBase/nSeqSet
    )%>%
    ungroup()%>%
    add_row(#Add reference positions
      group_by(.,pos,weeksince)%>%
        summarise(
          base="Reference",#sREFch[pos],
          nBase=nSeqSet-sum(nBase),
          pBase=nBase/nSeqSet,
          nSeq=nSeqSet)
    )%>%
    #  arrange(base,pos,weeksince)%>% #TODO: delete the arrange and grouping. 
    ungroup()%>%
    # group_by(pos,base)%>%
    mutate(
      oneMinuspB=1-pBase,
      Hx=pBase*oneMinuspB
    )%>%ungroup()
  return(localProps)
}


calculateLogKsSingle<-function(propDF,focalWeek){

  focalWeek<-propDF[propDF$weeksince==focalWeek,]%>%
    mutate(Hy=Hx,
           pBasey=pBase
           )%>%
    select(pos,base,Hy,pBasey)
  
  logKdf1<-(propDF%>%
              arrange(base,pos,weeksince)%>%
             # group_by(base,pos)%>% #TODO: delete?
              left_join(.,focalWeek,by=c("pos","base")
              )%>%group_by(base)%>%
             # mutate(pBasey = replace(pBasey, base!="Reference"&is.na(pBasey),0))%>%
            #  mutate(pBasey = replace(pBasey, base=="Reference"&is.na(pBasey),1))%>%
            #  replace_na(list(Hy = 0))%>%
              mutate(
                Hp=Hx*Hy,
                meanProp=(pBase+pBasey)/2,
                Hs=(meanProp*(1-meanProp))*(meanProp*(1-meanProp))#,#^2
              ))%>%
    group_by(weeksince)%>%
    summarise(naHp=sum(is.na(Hp)),
              sumHp=sum(Hp,na.rm=T),
              naHs=sum(is.na(Hs)),
              sumHs=sum(Hs,na.rm=T),
              logK=log(sumHs/sumHp)
    )
  return(logKdf1) 
}

calculateLogKs<-function(propDF,lagVec=2:9){
  propDF2<-
    data.frame(
      weeksince=rep(propDF$weeksince,each=length(lagVec)),
      pos=rep(propDF$pos,each=length(lagVec)),
      base=rep(propDF$base,each=length(lagVec)),
      nBase=rep(propDF$nBase,each=length(lagVec)),
      pBase=rep(propDF$pBase,each=length(lagVec)),
      nSeq=rep(propDF$nSeq,each=length(lagVec)),
      oneMinuspB=rep(propDF$oneMinuspB,each=length(lagVec)),
      Hx=rep(propDF$Hx,each=length(lagVec)),
      lagT=rep(lagVec,nrow(propDF))
    )
  
  logKdf1<-(propDF2%>%
             arrange(lagT,base,pos,weeksince)%>%
             group_by(lagT)%>% #TODO: delete?
             mutate(
               Hy=lag(Hx,n=lagT[1]),
               Hp=Hx*Hy,
               meanProp=(pBase+lag(pBase,n=lagT[1]))/2,
               Hs=(meanProp*(1-meanProp))*(meanProp*(1-meanProp)),#^2
               bool=(pos==lag(pos,n=lagT[1]))&(base==lag(base,n=lagT[1]))&(lagT==lag(lagT,n=lagT[1]))&(weeksince==(lag(weeksince,n=lagT[1])+lagT[1])),
               emptyPrevRec=(lagT==lag(lagT,n=lagT[1]))&(weeksince!=(lag(weeksince,n=lagT[1])+lagT[1]))
             )%>%filter(bool)%>%
             group_by(weeksince,lagT)%>%
             summarise(naHp=sum(is.na(Hp)),
                       sumHp=sum(Hp,na.rm=T),
                       naHs=sum(is.na(Hs)),
                       sumHs=sum(Hs,na.rm=T),
                       logK=log(sumHs/sumHp)
             )
            )
  return(logKdf1) 
}

estimDeltaLogK<-function(SNPdataDF,lagVec=0:9,splitLogK=F){
  print("calculating proportions")
  allProps<-calculateAllProportions(SNPdataDF)

  
  print("calculating logK")
  if(splitLogK){
    print("going through each lag separately...")
    logKdf<-bind_rows(lapply(lagVec,function(x){
      print(paste0("Lag ",x))
      calculateLogKs(allProps,lagVec=x)
      }))
  } else {
    logKdf<-calculateLogKs(allProps,lagVec=lagVec)
  }
  print("calculating delta logK")
  dLk<-logKdf%>%group_by(weeksince)%>%
    summarise(meanLk=mean(diff(logK,lag=1)),
            maxLk=max(diff(logK,lag=1)),
            meanHp=mean(sumHp),
            meanHs=mean(sumHs)
    )%>%as.data.frame()
  return(dLk)
}

imputeNonMatched<-function(tempCombo){
  unMatchedx<-is.na(tempCombo$nPosBase.x)
  unMatchedRefs<-is.na(tempCombo$nPosBase.x)&(tempCombo$nOthers.y>0)
  unMatchedNonRefs<-is.na(tempCombo$nPosBase.x)&(tempCombo$nOthers.y==0)
  
  unMatchedy<-is.na(tempCombo$nPosBase.y)
  unMatchedRefsy<-is.na(tempCombo$nPosBase.y)&(tempCombo$nOthers.x>0)
  unMatchedNonRefsy<-is.na(tempCombo$nPosBase.y)&(tempCombo$nOthers.x==0)
  
  tempCombo[unMatchedx,"nSeq.x"]<-mean(tempCombo$nSeq.x,na.rm=T)
  tempCombo[unMatchedx,"t.x"]<-mean(tempCombo$t.x,na.rm=T)
  tempCombo[unMatchedx,"pos.x"]<-tempCombo[unMatchedx,"pos.y"]
  tempCombo[unMatchedx,"base.x"]<-tempCombo[unMatchedx,"base.y"]
  tempCombo[unMatchedx,"nOthers.x"]<-0
  
  tempCombo[unMatchedRefs,"nPosBase.x"]<-mean(tempCombo$nSeq.x,na.rm=T)
  tempCombo[unMatchedRefs,"propPosBase.x"]<-1
  tempCombo[unMatchedNonRefs,"nPosBase.x"]<-0
  tempCombo[unMatchedNonRefs,"propPosBase.x"]<-0
  
  tempCombo[unMatchedy,"nSeq.y"]<-mean(tempCombo$nSeq.y,na.rm=T)
  tempCombo[unMatchedy,"t.y"]<-mean(tempCombo$t.y,na.rm=T)
  tempCombo[unMatchedy,"pos.y"]<-tempCombo[unMatchedy,"pos.x"]
  tempCombo[unMatchedy,"base.y"]<-tempCombo[unMatchedy,"base.x"]
  tempCombo[unMatchedy,"nOthers.y"]<-0
  
  tempCombo[unMatchedRefsy,"nPosBase.y"]<-mean(tempCombo$nSeq.y,na.rm=T)
  tempCombo[unMatchedRefsy,"propPosBase.y"]<-1
  tempCombo[unMatchedNonRefsy,"nPosBase.y"]<-0
  tempCombo[unMatchedNonRefsy,"propPosBase.y"]<-0
  return(tempCombo)
}

calcLogK<-function(subset1,subset2){
  tempCombo<-merge(addReference(subset1),
                   addReference(subset2),
                   by="posBase",all=T
  )
  tempCombo<-imputeNonMatched(tempCombo)
  
  tempCombo$propPosBase.both<-(tempCombo$propPosBase.x+tempCombo$propPosBase.y)/2 #Assuming the isolates are a representative sample of two equally sized populations.
  
  tempCombo$px<-tempCombo$propPosBase.x*(1-tempCombo$propPosBase.x)
  tempCombo$py<-tempCombo$propPosBase.y*(1-tempCombo$propPosBase.y)
  tempCombo$pb<-tempCombo$propPosBase.both*(1-tempCombo$propPosBase.both)
  
  Hp<-sum(tempCombo$px*tempCombo$py)
  Hs<-sum(tempCombo$pb*tempCombo$pb)
  return(log(Hs/Hp))
}

#system.time(
createNxNlogKmatrix<-function(snpFreqs,startW=50,endW=146){
  timelogK<-bind_rows(
    lapply(startW:endW,function(y){
        bind_rows(
        #mclapply(startW:endW, function(x){
#        parLapply(cl, startW:endW,function(x){
         lapply(startW:endW, function(x){
            print(paste0(x,"->",y))
            return(data.frame(x=x,y=y,fsp=calcLogK(snpFreqs[snpFreqs$t==x,],snpFreqs[snpFreqs$t==y,])))
          }
        )
      )}
    )
  )
  
  timelogK<-
    timelogK%>%
    group_by(y)%>%
    #arrange(x,by_group=T)%>%
    mutate(logitFSP=fsp,
           diffTime=x-y,
           diffTimeSign=diffTime/abs(diffTime),
           diffLfsp=diffTimeSign*((logitFSP)-lag(logitFSP,1))/(1),
           diffLfsp2=diffTimeSign*((logitFSP)-lag(logitFSP,2))/(2),
           diffLfsp3=diffTimeSign*((logitFSP)-lag(logitFSP,3))/(3),
           diffLfspRev=-diffTimeSign*((logitFSP)-lead(logitFSP,1))/(1),
           diffLfsp2Rev=-diffTimeSign*((logitFSP)-lead(logitFSP,2))/(2),
           diffLfsp3Rev=-diffTimeSign*((logitFSP)-lead(logitFSP,3))/(3)
    )%>%as.data.frame()
  
  return(timelogK)
}

fitnessEstimate<-function(logKmatrix,serialInterval=5, estiStart=-7, estiEnd=-2,threshold=0.125){
  averageL<-logKmatrix%>%group_by(y)%>%
    filter(diffTime>=(estiStart),diffTime<=(estiEnd))%>%
    summarise(maxdiff=max(diffLfspRev,na.rm=T),
              meandiff=median(diffLfspRev,na.rm=T)
              
              )%>%
    as.data.frame()
  
  averageL<-averageL[(3):nrow(averageL),]
  averageL$fitEstiMax<-(exp(averageL$maxdiff*serialInterval/7)-1)
  averageL$fitEsti<-(exp(averageL$meandiff*serialInterval/7)-1)
  averageL$above05<-(exp(averageL$meandiff*serialInterval/7)-1)>threshold
  averageL$lastAbove<-(lead(averageL$above05)!=averageL$above05)&(averageL$above05)
  averageL$firstAbove<-(lag(averageL$above05)!=averageL$above05)&(averageL$above05)
  averageL$date<-fConvWtoDate(averageL$y)
  
  averageL[which(is.na(averageL$lastAbove)),"lastAbove"]<-TRUE
  averageL[which(is.na(averageL$firstAbove)),"firstAbove"]<-averageL[which(is.na(averageL$firstAbove)),"above05"]
  return(averageL)
}

onlyFitnessEstimate<-function(snpFreqs,startW=50,endW=146,serialInterval=5, estiStart=-7, estiEnd=-2,threshold=0.125){
  timelogK<-bind_rows(
    lapply(startW:endW,function(y){
      bind_rows(
        #mclapply(startW:endW, function(x){
        #        parLapply(cl, startW:endW,function(x){
        lapply((y+estiStart):(y+estiEnd), function(x){
          print(paste0(x,"->",y))
          return(data.frame(x=x,y=y,fsp=calcLogK(snpFreqs[snpFreqs$t==x,],snpFreqs[snpFreqs$t==y,])))
        }
        )
      )}
    )
  )
  
  timelogK<-
    timelogK%>%
    group_by(y)%>%
    #arrange(x,by_group=T)%>%
    mutate(logitFSP=fsp,
           diffTime=x-y,
           diffTimeSign=diffTime/abs(diffTime),
           diffLfsp=diffTimeSign*((logitFSP)-lag(logitFSP,1))/(1),
           diffLfsp2=diffTimeSign*((logitFSP)-lag(logitFSP,2))/(2),
           diffLfsp3=diffTimeSign*((logitFSP)-lag(logitFSP,3))/(3),
           diffLfspRev=-diffTimeSign*((logitFSP)-lead(logitFSP,1))/(1),
           diffLfsp2Rev=-diffTimeSign*((logitFSP)-lead(logitFSP,2))/(2),
           diffLfsp3Rev=-diffTimeSign*((logitFSP)-lead(logitFSP,3))/(3)
    )%>%as.data.frame()
  
  averageL<-timelogK%>%group_by(y)%>%
    filter(diffTime>=(estiStart),diffTime<=(estiEnd))%>%
    summarise(maxdiff=max(diffLfspRev,na.rm=T),
              meandiff=median(diffLfspRev,na.rm=T)
              
    )%>%
    as.data.frame()
  
  averageL<-averageL[(3):nrow(averageL),]
  averageL$fitEstiMax<-(exp(averageL$maxdiff*serialInterval/7)-1)
  averageL$fitEsti<-(exp(averageL$meandiff*serialInterval/7)-1)
  averageL$above05<-(exp(averageL$meandiff*serialInterval/7)-1)>threshold
  averageL$lastAbove<-(lead(averageL$above05)!=averageL$above05)&(averageL$above05)
  averageL$firstAbove<-(lag(averageL$above05)!=averageL$above05)&(averageL$above05)
  averageL$date<-fConvWtoDate(averageL$y)
  
  averageL[which(is.na(averageL$lastAbove)),"lastAbove"]<-TRUE
  averageL[which(is.na(averageL$firstAbove)),"firstAbove"]<-averageL[which(is.na(averageL$firstAbove)),"above05"]
  return(averageL)
}

extractPeriods<-function(estiDF,threshold=0.125,thresholdOnMax=T){
  estiDF$above05<-estiDF$maxLk>threshold
  if(!thresholdOnMax) estiDF$above05<-estiDF$meanLk>threshold
  estiDF$lastAbove<-(lead(estiDF$above05)!=estiDF$above05)&(estiDF$above05)
  estiDF$firstAbove<-(lag(estiDF$above05)!=estiDF$above05)&(estiDF$above05)
  estiDF$date<-fConvWtoDate(estiDF$weeksince)
  print(estiDF)
  print(estiDF[estiDF$firstAbove,])
  print(estiDF[estiDF$lastAbove,])
  data.frame(
    firsts=estiDF[estiDF$firstAbove,"date"],
    lasts=estiDF[estiDF$lastAbove,"date"],
    firstWk=estiDF[estiDF$firstAbove,"weeksince"],
    lastWk=estiDF[estiDF$lastAbove,"weeksince"]
  )  
}

#Plot of the estimated fitness advantage
fitnessPlot<-function(estiDF,startTF=1,endTF=6,serialInterval=5,inclMaxLine=F,threshold=0.3,thresholdOnMax=T){
  estiDF$aboveThreshold<-estiDF$maxLk>threshold
  if(!thresholdOnMax) estiDF$aboveThreshold<-estiDF$meanLk>threshold
  estiDF$lastAbove<-(lead(estiDF$aboveThreshold)!=estiDF$aboveThreshold)&(estiDF$aboveThreshold)
  estiDF$firstAbove<-(lag(estiDF$aboveThreshold)!=estiDF$aboveThreshold)&(estiDF$aboveThreshold)
  estiDF$date<-fConvWtoDate(estiDF$weeksince)
  ymaxEsti<-max(1,max(estiDF$maxLk[!(estiDF$maxLk%in%c(Inf,-Inf,NA,NaN))])*1.1)
  print(ymaxEsti)
  ggplot(estiDF)+
    geom_line(aes(x=date,y=(meanLk)),color="blue",size=0.8)+
    geom_line(aes(x=date,y=(maxLk)),color="red",size=0.8)+
  #  geom_line(data=NxNdf[NxNdf$diffTime>=(startTF)&NxNdf$diffTime<=(endTF),],
  #            aes(x=fConvWtoDate(weeksince),y=(diffLfsp2Rev),group=x),color="#00000040")+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = 0.25*7/5.5,linetype=2)+
    geom_hline(yintercept = 0.5*7/5.5,linetype=2)+
    geom_hline(yintercept = 0.75*7/5.5,linetype=2)+
   # geom_hline(yintercept = 1.0*7/5.5,linetype=2)+
    #  geom_hline(yintercept = log(1+2)*7/5,linetype=2)+
    
    
    annotate(geom="label", x=as.Date("2021-10-21"), y=0.25*7/5.5, label="+25%",
             color="black",fill="white",label.size = NA)+
    annotate(geom="label", x=as.Date("2021-10-21"), y=0.5*7/5.5, label="+50%",
             color="black",fill="white",label.size = NA)+
    annotate(geom="label", x=as.Date("2021-10-21"), y=0.75*7/5.5, label="+75%",
             color="black",fill="white",label.size = NA)+
    (if(ymaxEsti>1*7/5.5){geom_hline(yintercept = 1*7/5.5,linetype=2)})+
    (if(ymaxEsti>1*7/5.5){annotate(geom="label", x=as.Date("2021-10-21"), y=7/5.5, label="+100%",
                               color="black",fill="white",label.size =NA)})+
    (if(ymaxEsti>1.25*7/5.5){geom_hline(yintercept = 1.25*7/5.5,linetype=2)})+
    (if(ymaxEsti>1.25*7/5.5){annotate(geom="label", x=as.Date("2021-10-21"), y=1.25*7/5.5, label="+125%",
                             color="black",fill="white",label.size =NA)})+
      #annotate(geom="label", x=as.Date("2021-10-21"), y=0.75*7/5.5, label="+75%",
      #                     color="black",fill="white",label.size = NA)
    #} else {theme_bw()})+
  #  annotate(geom="label", x=as.Date("2021-10-21"), y=5/7, label="+100%",
  #           color="black",fill="white",label.size = NA)+
    #annotate(geom="label", x=as.Date("2021-10-01"), y=log(1+2.0)*7/5, label="+200%",
    #         color="black",fill="white",label.size = NA)+
    scale_x_date(limits=c(as.Date("2021-01-01")-4,as.Date("2022-10-01")),
                 minor_breaks = as.Date("2020-12-01")+months(1:24),breaks=c(),
                 sec.axis = dup_axis())+
    theme_bw()+
    geom_rect(data=data.frame(
      firsts=estiDF[estiDF$firstAbove,"date"]-days(3),
      lasts=estiDF[estiDF$lastAbove,"date"]+days(3)),
      aes(xmin=firsts,xmax=lasts,ymin=-.25,ymax=ymaxEsti),fill="#00000020",color="#00000000"
    )+
    # theme(panel.grid.minor = element_line(colour="blue", size=0.5))+
    #theme(axis.line = element_rect(colour = "black"))+
    labs(x="",y=expression(Delta*log(k)))
}

timelineStiExamplePlot<-function(estiDF,startT=-9,endT=-2){
  refWeek<-estiDF[estiDF$logK==0,"weeksince"]
  refWeek<-refWeek[!is.na(refWeek)]
  print(refWeek)
  subdata=as.data.frame(estiDF[estiDF$weeksince>=(refWeek+startT)& estiDF$weeksince<=(refWeek+endT),])
  subdata$deltaL<-subdata$logK-lag(subdata$logK)
  print(subdata)
  steepSlopeID<-(which(subdata$deltaL==min(subdata$deltaL,na.rm = T)))[1]
  print(as.numeric(fConvWtoDate(subdata[steepSlopeID,"weeksince"])))
  print(subdata[steepSlopeID,"deltaL"])
  
  
  maxLine<-data.frame(weeksince=c(subdata[steepSlopeID,"weeksince"]-3,subdata[steepSlopeID,"weeksince"]+3),
             logK=c(subdata[steepSlopeID,"logK"]-(3*subdata[steepSlopeID,"deltaL"]),
                    subdata[steepSlopeID,"logK"]+(3*subdata[steepSlopeID,"deltaL"]))
             )
  meanLine<-data.frame(weeksince=c(refWeek+startT-1,refWeek+endT+1),
                       logK=c(subdata[1,"logK"]-(1*mean(subdata$deltaL,na.rm=T)),
                              subdata[1,"logK"]+(((endT-startT)+1)*mean(subdata$deltaL,na.rm=T)))
                       )
  print(meanLine)
  
ggplot(estiDF)+
  geom_line(aes(x=fConvWtoDate(weeksince),y=logK),color="#AAAAAA")+
  geom_line(data=subdata,aes(x=fConvWtoDate(weeksince),y=logK),color="black")+
  geom_point(data=estiDF[estiDF$weeksince>=(refWeek+startT)& estiDF$weeksince<=(refWeek+endT),],aes(x=fConvWtoDate(weeksince),y=logK))+
  geom_point(data=estiDF[estiDF$weeksince<(refWeek+startT)|estiDF$weeksince>(refWeek+endT),],aes(x=fConvWtoDate(weeksince),y=logK),shape=21,color="#999999",fill="white")+
  geom_line(data=maxLine,aes(x=fConvWtoDate(weeksince),y=logK),color="red")+
  geom_line(data=meanLine,aes(x=fConvWtoDate(weeksince),y=logK),color="blue")+
  #geom_point(data=estiDF[estiDF$y==exampleT&NxNdf$diffTime>0,],aes(x=fConvWtoDate(x),y=fsp),color="#999999")+
  geom_vline(xintercept = as.numeric(fConvWtoDate(refWeek)),linetype=2)+
  geom_vline(xintercept = as.numeric(fConvWtoDate(refWeek+endT)),linetype=1)+
  geom_vline(xintercept = as.numeric(fConvWtoDate(refWeek+startT)),linetype=1)+
  scale_x_date(limits=c(fConvWtoDate(refWeek-26),fConvWtoDate(refWeek)))+
  labs(x="Date",y="log(k)")+
  theme_bw()+
  theme(aspect.ratio = 1)
}

fTheory<-function(p,L){ ((p* exp(L))/(1 + (-1 + exp(L))*p)) }

fLogitFSPtheory0Alexis<-function(p,L) {
  log((
    ((((p+fTheory(p,L))/2)*(1-((p+fTheory(p,L))/2)))^2)/
      (
        p*(1-p)*fTheory(p,L)*(1-fTheory(p,L))
      )
  ))
}

flogKApprox<-function(p,L) {
  res<-L+ (2*log(1/4))+(2*log(p+1))
  res[L<=0]<-NA
  return(res)
}

flogKApprox2<-function(p,L) {
  #-L+log(((1/4)^2))+log((2-p)^2)
  res<-(-L+ (2*log(1/4))+(2*log(2-p)))
  res[L>=0]<-NA
  return(res)
}

simulateFSP<-function(p=0.5,varvar=0.0,varwt=0.0,beginL=-10,endL=10,stepsize=0.1){
  L=seq(beginL,endL+stepsize,stepsize)
  pNew=convertOddsToProp(exp(logit(p)+L))
  pSNP=(pNew*(1-varvar))+((1-pNew)*varwt)
  
  print(pSNP)
  zeropoint<-round((abs(beginL)/stepsize)+1)
  print(zeropoint)
  simDF<-data.frame(L=L,
             p=pNew,
             pSNP=pSNP)
  
  simDF$Hx<-simDF$pSNP*(1-simDF$pSNP)
  simDF$Hy<-simDF[zeropoint,"Hx"]
  simDF$meanProp<-(simDF$pSNP+simDF[zeropoint,"pSNP"])/2
  simDF$Hp=simDF$Hx*simDF$Hy
  simDF$Hs=(simDF$meanProp*(1-simDF$meanProp))*(simDF$meanProp*(1-simDF$meanProp))#,#^2
  simDF$logK=log(simDF$Hs/simDF$Hp)
  return(simDF)
}

theoryPlotVar<-function(p=0.5,beginL=-10,endL=10,stepSize=0.1){
  theoryData<-data.frame(
    L=rep(seq(beginL,endL+stepSize,0.1),length(p)),
    fsp0=unlist(lapply(p,function(x)fLogitFSPtheory0Alexis(x,seq(beginL,endL+stepSize,0.1)))),
    fsp1=unlist(lapply(p,function(x)flogKApprox(x,seq(beginL,endL+stepSize,0.1)))),
    fsp2=unlist(lapply(p,function(x)flogKApprox2(x,seq(beginL,endL+stepSize,0.1)))),
    p=as.factor(rep(p,each=(((endL+stepSize)-beginL)/stepSize)+1))
  )
  
  ggplot(theoryData)+
    geom_line(data=simulateFSP(beginL=-10,endL=10,varvar = 0.02) ,aes(x=L,y=logK,color="Partial Variant",linetype=as.factor(0.02)))+
    geom_line(data=simulateFSP(beginL=-10,endL=10,varvar = 0.01),aes(x=L,y=logK,color="Partial Variant",linetype=as.factor(0.01)))+
    geom_line(data=simulateFSP(beginL=-10,endL=10,varwt = 0.01),aes(x=L,y=logK,color="Partial part of wildtype",linetype=as.factor(0.01)))+
    geom_line(data=simulateFSP(beginL=-10,endL=10,varwt = 0.02),aes(x=L,y=logK,color="Partial part of wildtype",linetype=as.factor(0.02)))+
    geom_line(aes(x=L,y=fsp1,color="Approximation",group=p),linetype=2)+
    geom_line(aes(x=L,y=fsp2,color="Approximation",group=p),linetype=2)+
    scale_y_continuous(limits=c(-2,8))+
    scale_color_manual(limits=c(
      "Partial Variant",
      "Partial part of wildtype",
      "Approximation"
    ),values=c(
      "red","blue","black"
    ))+
    labs(y="ln(k)",x="L (fitness advantage x t)",linetype="Proportion",colour="Type of variation")+
    theme_bw()+
    theme(aspect.ratio = 1)
}



theoryPlot<-function(p=0.5,beginL=-10,endL=10){
  stepSize=0.1
  print(rep(p,each=(((endL+stepSize)-beginL)/stepSize)+1))
  
  theoryData<-data.frame(
    L=rep(seq(beginL,endL+stepSize,0.1),length(p)),
    fsp0=unlist(lapply(p,function(x)fLogitFSPtheory0Alexis(x,seq(beginL,endL+stepSize,0.1)))),
    fsp1=unlist(lapply(p,function(x)flogKApprox(x,seq(beginL,endL+stepSize,0.1)))),
    fsp2=unlist(lapply(p,function(x)flogKApprox2(x,seq(beginL,endL+stepSize,0.1)))),
    p=as.factor(rep(p,each=(((endL+stepSize)-beginL)/stepSize)+1))
  )
    
  ggplot(theoryData)+
  #  geom_line(aes(x=L,y=fsp1),color="gray")+
#  geom_line(aes(x=L,y=fsp,color="logit(Fsp)"))+
  geom_line(aes(x=L,y=fsp0,color="log(k)",group=p,linetype=p))+
  geom_line(aes(x=L,y=fsp1,color="Approximation",group=p,linetype=p))+
  geom_line(aes(x=L,y=fsp2,color="Approximation",group=p,linetype=p))+
    
  #geom_line(data=data.frame(L=(-101:99)/10,
  #                          fsp=calculatedFST),
  #          aes(x=L,y=fsp,color="Simulated"),linetype=2,size=1)+
  scale_color_manual(limits=c("Approximation","log(k)"),
                     values=c("red","black")
  )+
  scale_y_continuous(limits=c(-2,8))+
  labs(y="ln(k)",x="L (fitness advantage x t)")+
  theme_bw()+
  theme(aspect.ratio = 1)
}

resampledLK<-function(snpDF,sampleSize,repeats=100,verbose=F){
lapply(1:repeats,function(x){
  if(verbose)print(paste0("loop ",x))
  subdata<-samplePerWeek(snpDF,size=sampleSize)
  tim<-system.time(
    dLk<-estimDeltaLogK(subdata,nSeqSet=sampleSize,lagVec = 2:8)
  )
  dLk$loop<-x
  dLk$sampleSize<-sampleSize
  dLk$timing<-tim[3]
  return(dLk)
})%>% bind_rows()
}
