getPosNegLineages<-function(startW,endW,snpdata,propdata,metaData,SeqMeta,scoreCut=0){
  #rownames(SeqMeta)<-SeqMeta$seqID
  startDate=fConvWtoDate(startW)
  endDate=fConvWtoDate(endW+1)
  
  selectedSeqIds<-metaData[(metaData$date>=startDate&metaData$date<endDate),"seqID"]
  selectedSeqIds<-selectedSeqIds[!is.na(selectedSeqIds)]

  actualSeqs <- snpdata[snpdata$seqID %in% selectedSeqIds, ]
  actualSeqs$posBase<-paste0(actualSeqs$pos,actualSeqs$base)
  propdata$posBase<-paste0(propdata$pos,propdata$base)
  
  w1=propdata%>%
    filter(weeksince==startW)%>%
    mutate(pBase1=pBase)%>%
    select(posBase,pBase1)
    
  thisTwoWks=propdata%>%
    filter(weeksince==endW)%>%
    mutate(pBase2=pBase)%>%
    select(posBase,pBase2)%>%
    full_join(w1)%>%
    filter((!is.na(pBase1)&(!is.na(pBase2))))%>%
    filter(pBase1>0)%>%filter(pBase2>0)%>%
    filter(pBase1<1)%>%filter(pBase2<1)
  
  print("Done with the merger")
  
#  thisTwoWks<-getTwoWeekCombo(snpFreqs[snpFreqs$t==startW,],snpFreqs[snpFreqs$t==(endW),])
  #actualSeqs%>%group_by(posbase,)
  print("Calculating L score per posbase")
  thisTwoWks$L<-logit(thisTwoWks$pBase2)-logit(thisTwoWks$pBase1)
  thisTwoWks$logit1<-logit(thisTwoWks$pBase1)
  thisTwoWks$logit2<-logit(thisTwoWks$pBase2)
  
  print("Calculating the score of the reference")
  #refScore<-merge(data.frame(
  #  pos=1:length(sREFch),
  #  posBase=paste0(1:length(sREFch),sREFch)
  ##),
  #print( paste0(1:length(sREFch),sREFch))
  print( thisTwoWks$posBase)
#  refScoreEntries<-(thisTwoWks[thisTwoWks$posBase %in% paste0(1:length(sREFch),sREFch),"L"])
  refScoreEntries<-(thisTwoWks[thisTwoWks$posBase %in% paste0(1:length(sREFch),"Reference"),])
  
  refScore<-sum(refScoreEntries$L)
  print(refScoreEntries)
  #thisTwoWks[,c("posBase","L")],by="posBase",all.x=F)
  print(paste0("Refscore is ",refScore, "from ",nrow(refScoreEntries)," entries"))
  #thisTwoWks[thisTwoWks$incl,c("posBase","L")],by="posBase",all.x=F)
  #actualSeqsPlusScores<-merge(actualSeqs,thisTwoWks[thisTwoWks$incl,c("posBase","L")],by="posBase",all.x=T)
  
  print("Calculating the score of the sequences")
  
 # snpdata<-snpdata[snpdata$seqID %in% selectedSeqIds$seqID,]
  getScoreForSNPs<-function(posBases){
    sum(thisTwoWks[thisTwoWks$posBase %in% posBases,"L"])
  }
  rownames(SeqMeta)<-SeqMeta$Accession.ID
  scoredSeqs<-actualSeqs%>%
    group_by(seqID)%>%
    summarise(
      n=n(),
      score=getScoreForSNPs(posBase)#
      #Pango=SeqMeta[metaData[seqID[1],"EPI"],"Pango.lineage"]
    )%>%ungroup()%>%
    mutate(Pango=SeqMeta[metaData[seqID,"EPI"],"Pango.lineage"])%>%
    as.data.frame()
  
  print(scoredSeqs)
  print("Calculating the score of the sequences, step 2")
  
  if(is.null(scoreCut))scoreCut<-mean(scoredSeqs$score)
  allPosiEPI<-metaData[metaData$seqID %in% ((scoredSeqs[scoredSeqs$score>scoreCut,"seqID"])),"EPI"]
  allNegiEPI<-metaData[metaData$seqID %in% ((scoredSeqs[scoredSeqs$score<scoreCut,"seqID"])),"EPI"]
  print("Exporting")
  
  return(
    list(
      refScore=refScore,
      scoredSeqs=scoredSeqs,
      allPosiEPI=allPosiEPI,
      posiPango=SeqMeta[SeqMeta$Accession.ID %in% allPosiEPI,"Pango.lineage"],
      allNegiEPI=allNegiEPI,
      negiPango=SeqMeta[SeqMeta$Accession.ID %in% allNegiEPI,"Pango.lineage"]
    )
  )
}

#allPeriods

