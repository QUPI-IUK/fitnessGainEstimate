fConvWtoDate<-function(x)as.Date("2020-01-06")+weeks(x)

weeksSince<-function(w){
  data.frame(
    weeksSince=w+1,
    weekPad=paste0(isoyear(ymd("2020-01-06")+days(7*w)),"-",
                   str_pad(isoweek(ymd("2020-01-06")+days(7*w)),width=2,side="left",pad="0")
    ),
    weekNonPad=paste0(isoyear(ymd("2020-01-06")+days(7*w)),"-",
                      isoweek(ymd("2020-01-06")+days(7*w))
    )
  )
}

weeksinceToDate<-function(w){
  return(as.Date("2020-01-06")+days(7*w))
}

convertMetaTSVtoRDS<-function(metaFilename,blockSize=1000000){
  system(paste("split","-l", format(blockSize, scientific=FALSE),# 100000,",#blockSize*2,
               metaFilename,
               paste0(dirname(metaFilename),"/tempSplit-Seg")))
  segmentFileList=list.files(dirname(metaFilename),"tempSplit-Seg")
  print(paste0("Found ",length(segmentFileList)," split files"))
  nBlock<-length(segmentFileList)
  print(segmentFileList)
  detectedColNames<-NULL
  allData<-lapply(1:nBlock, function(bl){
    #pbmclapply(1:2, function(bl){
    # lapply(1:nBlock, function(bl){
    print(bl)
    if(is.null(detectedColNames)){
      metaTemp<-read.csv(paste0(dirname(metaFilename),"/",segmentFileList[bl]),sep="\t")
      detectedColNames<<-colnames(metaTemp)
    } else {
      metaTemp<-read.csv(paste0(dirname(metaFilename),"/",segmentFileList[bl]),sep="\t")
      colnames(metaTemp) <- detectedColNames
    }
    return(metaTemp[,c("Virus.name","Accession.ID","Collection.date","Clade","Pango.lineage","Submission.date")])
  })%>%bind_rows()
  print(paste0("Loaded ",nrow(allData)," sequence meta data records"))
  saveRDS(allData,"Resources/allMetaData.rds")
  #file.remove(paste0(dirname(msaFilename),"/",segmentFileList))
}

bind_rows(lapply(c("bagbaga"),function(x)data.frame(x=1:2,y=4:5)))

loadSeqMetaData<-function(geoString=""){ #TODO: This is still based on German data, we now have the world.
  print("Loading metadata file...")
  metaData<-readRDS("Resources/allMetaData.rds")
  if(geoString!=""){
     metaData<-bind_rows(lapply(geoString,function(x){
       print(x)
       metaData[str_detect(metaData$Virus.name,x),]
     }
    )
  )
  }
  
  print(nrow(metaData))
  rownames(metaData)<-metaData$Accession.ID
  metaData$runningWeek<-ceiling(as.numeric(as.Date(metaData$Collection.date)-as.Date("2020-01-06"))/7)
  print("Creating color scheme...")
  colorSchemeLineages<-createLineageColorscheme(metaData)
  metaData<-metaData%>%
    mutate(weeksince=as.numeric(floor((as.Date(Collection.date)-as.Date("2020-01-06"))/7)),
           weeksinceSubmit=as.numeric(floor((as.Date(Submission.date)-as.Date("2020-01-06"))/7))
    )%>%
    left_join(colorSchemeLineages, by = join_by(Pango.lineage == fullName))
  
  return(metaData)
}

addReference<-function(tframe){
  #  tframe=snpFreqs[snpFreqs$t==50,]
  tframe$pos<-as.numeric(str_sub(tframe[,"posBase"],end=-2))
  tframe$base<-str_sub(tframe[,"posBase"],start=-1,end=-1)
  tframe$nOthers<-0
  
  tempJoin<-rbind(
    tframe%>% 
      select(t,nOthers,pos,nPosBase,propPosBase,base,posBase,nSeq)
    ,
    tframe%>% 
      select(t,nOthers,pos,nPosBase,propPosBase,base,posBase,nSeq)%>%
      group_by(pos,t)%>%
      summarise(
        #t=t,
        nOthers=n(),
        nPosBase=nSeq[1]-sum(nPosBase),
        propPosBase=nPosBase/nSeq[1],
        pos=pos,
        base=sREFch[pos],
        posBase=paste0(pos,base),
        nSeq=nSeq[1],
        .groups = "drop"
      )%>%as.data.frame()
  )
  return(as.data.frame(tempJoin[order(tempJoin$pos),]))
}

getSeqsBeforeSubmDate<-function(dataDF,metaData,seqMetaData,submDate){
  selectID<-metaData[
    metaData$EPI %in% seqMetaData[seqMetaData$Submission.date<as.Date(submDate),"Accession.ID"],"seqID"]
  
  subDataDF<-dataDF[dataDF$seqID %in% selectID,]
  
  existingIDs=(subDataDF%>%group_by(weeksince,seqID)%>%
                 summarise(y=1))
  weekCounts<-existingIDs%>%group_by(weeksince)%>%
    summarise(nSeq=n())
  
  subDataDF<-subDataDF%>%select(!nSeq)%>%left_join(weekCounts,by="weeksince")
  return(subDataDF)
}

getSeqsBeforeSubmWeek<-function(dataDF,metaData,seqMetaData,submWeek){
  selectID<-metaData[
    metaData$EPI %in% seqMetaData[seqMetaData$weeksinceSubmit<=submWeek,"Accession.ID"],"seqID"]
  subDataDF<-dataDF[dataDF$seqID %in% selectID,]
  
  existingIDs=(subDataDF%>%group_by(weeksince,seqID)%>%
                 summarise(y=1))
  weekCounts<-existingIDs%>%group_by(weeksince)%>%
    summarise(nSeq=n())
  
  subDataDF<-subDataDF%>%select(!nSeq)%>%left_join(weekCounts,by="weeksince")
  return(subDataDF)
}

samplePerWeek<-function(dataDF,size=10){
  existingIDs=(dataDF%>%group_by(weeksince,seqID)%>%
                 summarise(y=1))
  resIDs<-existingIDs%>%
    group_by(weeksince)%>%
    slice_sample(n=size,replace = T)
  
  subDataDF<-dataDF[dataDF$seqID %in% resIDs$seqID,]
  weekCounts<-resIDs%>%group_by(weeksince)%>%
    summarise(nSeq=n())
  print(weekCounts)
  subDataDF<-subDataDF%>%select(!nSeq)%>%left_join(weekCounts,by="weeksince")
  
  return(subDataDF)
}

tlLkQunatiles<-function(subSampled){
  subSampled%>%group_by(weeksince)%>%
    summarise(lowerRange=quantile(maxLk,0.0,na.rm=T),
              lowerIQR=quantile(maxLk,0.25,na.rm=T),
              median=quantile(maxLk,0.5,na.rm=T),
              upperIQR=quantile(maxLk,0.75,na.rm=T),
              upperRange=quantile(maxLk,1.0,na.rm=T)
    )
}
# system.time({
#   myFunction<-function(x){
#     print(paste0("Submission Date cut-off:  ",x))
#     allProps<-calculateAllProportions(getSeqsBeforeSubmDate(allSNPDataDump,allMetaData,fullSeqMetaData,x))
#     logKdf<-calculateLogKs(allProps,lagVec=0:9)
#     dLk<-logKdf%>%group_by(weeksince)%>%
#       summarise(meanLk=mean(diff(logK,lag=1)),
#                 maxLk=max(diff(logK,lag=1)),
#                 meanHp=mean(sumHp),
#                 meanHs=mean(sumHs),
#                 cutIOffDate=x
#       )%>%as.data.frame()
#     return(dLk)
#   }
#   timeSteps<-lapply(as.Date("2021-01-03")+seq(0,730,7),myFunction)#,mc.cores=6)
#   timeSteps<-timeSteps%>% bind_rows()
# }
# )

#sREF = readDNAStringSet("/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta")[1]
#sREFch = unlist(strsplit(as.character(sREF), split = ""))
#sREFch[[383]]