getPosNegLineages_old<-function(db,startW,endW,snpFreqs,metaData,isolateMetaData,scoreCut=0){
  
  startDate=as.Date("2020-01-01")+(7*startW)
  endDate=as.Date("2020-01-01")+(7*(endW+1))
  
  selectedSeqIds<-metaData[(metaData$date>=startDate&metaData$date<endDate),"seqID"]
  selectedSeqIds<-selectedSeqIds[!is.na(selectedSeqIds)]
  #selectedSeqIds
  listIDs<-paste0(c(selectedSeqIds),collapse = ", ")
  myQuery<-paste0("SELECT * FROM seqdata WHERE seqID IN (",listIDs,")")
  actualSeqs <- dbGetQuery(db, myQuery)
  actualSeqs$posBase<-paste0(actualSeqs$pos,actualSeqs$base)
  
  thisTwoWks<-getTwoWeekCombo(snpFreqs[snpFreqs$t==startW,],snpFreqs[snpFreqs$t==(endW),])
  thisTwoWks$L<-logit(thisTwoWks$propPosBase.y)-logit(thisTwoWks$propPosBase.x)
  
  refScore<-merge(data.frame(
    pos=1:length(sREFch),
    posBase=paste0(1:length(sREFch),sREFch)
  ),
  thisTwoWks[thisTwoWks$incl,c("posBase","L")],by="posBase",all.x=F)
  
  actualSeqsPlusScores<-merge(actualSeqs,thisTwoWks[thisTwoWks$incl,c("posBase","L")],by="posBase",all.x=T)
  actualSeqsPlusScores[is.na(actualSeqsPlusScores$L),"L"]<-0
  
  rownames(metaData)<-metaData$seqID
  
  scoredSeqs<-actualSeqsPlusScores%>%group_by(seqID)%>%
    summarise(
      n=n(),
      refScore=sum(refScore[!(refScore$pos %in% pos),"L"]),
      seqScore=sum(L),
      score=refScore+seqScore
    )%>%#mutate(
    #  EPI=metaData[seqID,"EPI"]#,
     #Pango=isolateMetaData[EPI,"Pango.lineage"]
    #%>%
    as.data.frame()
  
  #scoredSeqs$Pango<-isolateMetaData$Accession.ID %in% metaData[scoredSeqs$EPI,"Pango.lineage"]
  scoredSeqs<-merge(scoredSeqs,metaData,by="seqID",all.x=T)
  scoredSeqs<-merge(scoredSeqs,isolateMetaData,by.x="EPI",by.y="Accession.ID",all.x=T)
  
  if(is.null(scoreCut))scoreCut<-mean(scoredSeqs$score)
  allPosiEPI<-metaData[metaData$seqID %in% ((scoredSeqs[scoredSeqs$score>scoreCut,"seqID"])),"EPI"]
  allNegiEPI<-metaData[metaData$seqID %in% ((scoredSeqs[scoredSeqs$score<scoreCut,"seqID"])),"EPI"]
  
  return(
    list(
      refScore=refScore,
      scoredSeqs=scoredSeqs,
      allPosiEPI=allPosiEPI,
      posiPango=isolateMetaData[isolateMetaData$Accession.ID %in% allPosiEPI,"Pango.lineage"],
      posiPango2=scoredSeqs[scoredSeqs$score>scoreCut,"Pango"],
      allNegiEPI=allNegiEPI,
      negiPango=isolateMetaData[isolateMetaData$Accession.ID %in% allNegiEPI,"Pango.lineage"]
    )
  )
}

# getScorePlot<-function(p,scoreResults,timeperiods,lineageClustering){
#   bwidth<-round((max(scoreResults[[p]]$scoredSeq$score)-min(scoreResults[[p]]$scoredSeq$score)),digits = 0)/50
#   print(bwidth)
#   ggplot(scoreResults[[p]]$scoredSeqs, aes(x=score,fill=Pango.lineage))+
#     geom_bar(stat="bin",binwidth=bwidth)+
#     scale_fill_manual(
#       limits=lineageClustering$fullName,
#       values = lineageClustering$variantCol
#     )+
#     labs(x="Sequence total SNP advantage score",y="Count",title=paste0(timeperiods[p,"firsts"]," until ",timeperiods[p,"lasts"]))+
#     theme_bw()+
#     theme(legend.position = "None",
#           aspect.ratio = 1
#     )
# }

get2x2Table<-function(p,scoreResults,timeperiods,lineageClustering,showN=1){
  negiCount<-apply(lineageClustering[scoreResults[[p]]$negiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")],2,sum)
  posiCount<-apply(lineageClustering[scoreResults[[p]]$posiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")],2,sum)
  
  varNames<-names(posiCount[order(-posiCount)])[1:showN]
  negiCount[order(-posiCount)]
  
  temp1<-lineageClustering[scoreResults[[p]]$posiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")]
  temp1$SNP_score<-"Positive"
  temp2<-lineageClustering[scoreResults[[p]]$negiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")]
  temp2$SNP_score<-"Negative"
  
  temp<-rbind(temp1,temp2)
  temp$Variant<-"Other"#paste0("Non-",varName)
  
  lapply(1:showN,function(n){temp[temp[,varNames[[n]]],"Variant"]<<-varNames[[n]]})
  #temp[temp[,varName],"Variant"]<-varName
  TP<-(sum(temp$Variant!="Other"&temp$SNP_score=="Positive"))
  FN<-(sum(temp$Variant!="Other"&temp$SNP_score=="Negative"))
  
  FP<-(sum(temp$Variant=="Other"&temp$SNP_score=="Positive"))
  TN<-(sum(temp$Variant=="Other"&temp$SNP_score=="Negative"))
  sens<-TP/(TP+FN)
  spec<-TN/(TN+FP)
  print(paste0(timeperiods[p,"firsts"]," until ",timeperiods[p,"lasts"]))
  print(paste0("Sensitivity: ",round(sens*100,digits = 2),"%, specificity: ",round(spec*100,digits = 2),"%."))
  print(table(temp[, c("Variant", "SNP_score")]))
  
  #  print(table(rbind(temp1,temp2)[, c(names(posiCount[order(-posiCount)])[[1]], "SNP_score")]))
  contTable(temp[, c("Variant", "SNP_score")])
}
