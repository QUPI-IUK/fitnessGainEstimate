
createLineageColorscheme<-function(lineages){

  lineageClustering<-data.frame(
    fullName=unique(lineages$Pango.lineage),
    orig=F,
    alpha=F,
    Delta=F,
    Omicron1=F,
    Omicron2=F,
    Omicron3=F,
    Omicron5=F
  )

  rownames(lineageClustering)<-lineageClustering$fullName
  lineageClustering$Omicron1<-lineageClustering$Omicron1|str_detect(lineageClustering$fullName,"BA.1")
  lineageClustering$Omicron1<-lineageClustering$Omicron1|str_detect(lineageClustering$fullName,"B.1.1.529")
  lineageClustering$Omicron2<-lineageClustering$Omicron2|str_detect(lineageClustering$fullName,"BA.2")
  lineageClustering$Omicron3<-lineageClustering$Omicron3|str_detect(lineageClustering$fullName,"BA.3")
  lineageClustering$Omicron5<-lineageClustering$Omicron5|str_detect(lineageClustering$fullName,"BA.5")
  lineageClustering$Omicron5<-lineageClustering$Omicron5|str_detect(lineageClustering$fullName,"BE.")
  lineageClustering$Omicron5<-lineageClustering$Omicron5|str_detect(lineageClustering$fullName,"BF.")
  lineageClustering$Delta<-lineageClustering$Delta|str_detect(lineageClustering$fullName,"AY")
  lineageClustering$Delta<-lineageClustering$Delta|str_detect(lineageClustering$fullName,"B.1.617.2")
  lineageClustering$alpha<-lineageClustering$alpha|str_detect(lineageClustering$fullName,"B.1.1.7")
  lineageClustering$any<-(lineageClustering$alpha|lineageClustering$Delta|lineageClustering$Omicron1|lineageClustering$Omicron2|lineageClustering$Omicron3|lineageClustering$Omicron5)
  lineageClustering$none<-!lineageClustering$any

  lineageClustering<-lineageClustering[order(lineageClustering$fullName),]
  lineageClustering$variantCol<-"#888888"
  lineageClustering[lineageClustering$alpha,"variantCol"]<- darken(brewer.pal(7, "Paired")[1],amount=(((1:sum(lineageClustering$alpha))/sum(lineageClustering$alpha))/2)-0.25)
  lineageClustering[lineageClustering$Delta,"variantCol"]<- darken(brewer.pal(7, "Paired")[2],amount=(((1:sum(lineageClustering$Delta))/sum(lineageClustering$Delta))/2)-0.25)
  lineageClustering[lineageClustering$Omicron1,"variantCol"]<- darken(brewer.pal(7, "Paired")[3],amount=(((1:sum(lineageClustering$Omicron1))/sum(lineageClustering$Omicron1))/2)-0.25)
  lineageClustering[lineageClustering$Omicron2,"variantCol"]<- darken(brewer.pal(7, "Paired")[4],amount=(((1:sum(lineageClustering$Omicron2))/sum(lineageClustering$Omicron2))/2)-0.25)
  lineageClustering[lineageClustering$Omicron3,"variantCol"]<- darken(brewer.pal(7, "Paired")[5],amount=(((1:sum(lineageClustering$Omicron3))/sum(lineageClustering$Omicron3))/2)-0.25)
  lineageClustering[lineageClustering$Omicron5,"variantCol"]<- darken(brewer.pal(7, "Paired")[6],amount=(((1:sum(lineageClustering$Omicron5))/sum(lineageClustering$Omicron5))/2)-0.25)

  lineageClustering$variantOrder<-paste0("0_",lineageClustering$fullName)
  lineageClustering[lineageClustering$alpha,"variantOrder"]<-paste0("1_",lineageClustering[lineageClustering$alpha,"fullName"])
  lineageClustering[lineageClustering$Delta,"variantOrder"]<-paste0("2_",lineageClustering[lineageClustering$Delta,"fullName"])
  lineageClustering[lineageClustering$Omicron1,"variantOrder"]<-paste0("3_",lineageClustering[lineageClustering$Omicron1,"fullName"])
  lineageClustering[lineageClustering$Omicron2,"variantOrder"]<-paste0("4_",lineageClustering[lineageClustering$Omicron2,"fullName"])
  lineageClustering[lineageClustering$Omicron3,"variantOrder"]<-paste0("5_",lineageClustering[lineageClustering$Omicron3,"fullName"])
  lineageClustering[lineageClustering$Omicron5,"variantOrder"]<-paste0("6_",lineageClustering[lineageClustering$Omicron5,"fullName"])
  return(lineageClustering)
}

lineagePerWeekPlot<-function(metaData,estiDF,threshold=0.125,thresholdOnMax=T){
  estiDF$above05<-estiDF$maxLk>threshold
  if(!thresholdOnMax) estiDF$above05<-estiDF$meanLk>threshold
  estiDF$lastAbove<-(lead(estiDF$above05)!=estiDF$above05)&(estiDF$above05)
  estiDF$firstAbove<-(lag(estiDF$above05)!=estiDF$above05)&(estiDF$above05)
  estiDF$date<-fConvWtoDate(estiDF$weeksince)
  
  linCountPerWeek<-metaData%>%group_by(weeksince,Pango.lineage)%>%
    summarise(nLineage=n())%>%
    as.data.frame()

  linColorScheme<-unique(metaData[,c("Pango.lineage","variantCol","variantOrder")])
  row.names(linColorScheme)<-linColorScheme$Pango.lineage
  
  linCountPerWeek<-rbind(linCountPerWeek,
                       data.frame(weeksince=rep(8:max(metaData$weeksince,na.rm=T),each=length(unique(metaData$Pango.lineage))),
                                  Pango.lineage=rep(unique(metaData$Pango.lineage),max(metaData$weeksince,na.rm=T)-7),
                                  nLineage=0
                       ))

  linCountPerWeek<-linCountPerWeek%>%group_by(weeksince,Pango.lineage)%>%
    summarise(nLineage=sum(nLineage))%>%
    group_by(weeksince)%>%
    mutate(nIso=sum(nLineage),
           propLin=nLineage/nIso,
           orderName=linColorScheme[Pango.lineage,"variantOrder"]
    )%>%
    as.data.frame()


  ggplot(linCountPerWeek)+
    #geom_line(aes(x=date,y=percent_variant,group=variant))
    geom_area(aes(x=fConvWtoDate(weeksince),y=propLin*100,fill=orderName),color="#00000005")+
     geom_rect(data=data.frame(
       firsts=estiDF[estiDF$firstAbove,"date"]-3,
       lasts=estiDF[estiDF$lastAbove,"date"]+3),
       aes(xmin=firsts,xmax=lasts,ymin=0,ymax=100),fill="#00000020",color="#FFFFFF20"
     )+
    scale_x_date(limits=c(as.Date("2021-01-01"),as.Date("2022-10-01")),
                 minor_breaks = as.Date("2020-12-01")+months(1:24),#breaks=c(),
                 #sec.axis = dup_axis()
    )+
    scale_fill_manual(
      limits=linColorScheme$variantOrder,
      values = linColorScheme$variantCol
    )+
    theme_bw()+
    theme(legend.position = "None",
          aspect.ratio=1/1
    )+
    labs(x="Date",y="Percentage VoC")
}

getScorePlot<-function(scoredSeqs,linColorScheme,allPeriods,p){
  ggplot(scoredSeqs, aes(x=score,fill=Pango))+
    geom_bar(stat="bin",binwidth=5)+
    scale_fill_manual(
      limits=linColorScheme$fullName,
      values = linColorScheme$variantCol
    )+
    labs(x="Sequence total SNP advantage score",y="Count",title=paste0(allPeriods[p,"firsts"]," until ",allPeriods[p,"lasts"]))+
    theme_bw()+
    theme(legend.position = "None",
          aspect.ratio = 1
    )
}

getScorePlotList<-function(resAllLin,linColorScheme,allPeriods,p){
  ggplot(resAllLin[[p]]$scoredSeqs, aes(x=score,fill=Pango))+
    geom_bar(stat="bin",binwidth=5)+
    scale_fill_manual(
      limits=linColorScheme$fullName,
      values = linColorScheme$variantCol
    )+
    labs(x="Sequence total SNP advantage score",y="Count",title=paste0(allPeriods[p,"firsts"]," until ",allPeriods[p,"lasts"]))+
    theme_bw()+
    theme(legend.position = "None",
          aspect.ratio = 1
    )
}




### Construct 2x2 tables # NOTE: resPNLineagesAll needs own name.
get2x2Table<-function(linColorScheme,p,showN=1){
  negiCount<-apply(linColorScheme[resPNLineagesAll[[p]]$negiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")],2,sum)
  posiCount<-apply(linColorScheme[resPNLineagesAll[[p]]$posiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")],2,sum)
  
  varNames<-names(posiCount[order(-posiCount)])[1:showN]
  negiCount[order(-posiCount)]
  
  temp1<-linColorScheme[resPNLineagesAll[[p]]$posiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")]
  temp1$SNP_score<-"Positive"
  temp2<-linColorScheme[resPNLineagesAll[[p]]$negiPango,c("alpha","Delta","Omicron1","Omicron2","Omicron3","Omicron5","none")]
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
  print(paste0(allPeriods[p,"firsts"]," until ",allPeriods[p,"lasts"]))
  print(paste0("Sensitivity: ",round(sens*100,digits = 2),"%, specificity: ",round(spec*100,digits = 2),"%."))
  print(table(temp[, c("Variant", "SNP_score")]))
  
  #  print(table(rbind(temp1,temp2)[, c(names(posiCount[order(-posiCount)])[[1]], "SNP_score")]))
#  contTable(temp[, c("Variant", "SNP_score")])
}

effectiveRplot<-function(RTdata,estiDF,threshold=0.3,thresholdOnMax=T){
    estiDF$aboveThreshold<-estiDF$maxLk>threshold
    if(!thresholdOnMax) estiDF$aboveThreshold<-estiDF$meanLk>threshold
    estiDF$lastAbove<-(lead(estiDF$aboveThreshold)!=estiDF$aboveThreshold)&(estiDF$aboveThreshold)
    estiDF$firstAbove<-(lag(estiDF$aboveThreshold)!=estiDF$aboveThreshold)&(estiDF$aboveThreshold)
    estiDF$date<-fConvWtoDate(estiDF$weeksince)
    
  ggplot(RTdata)+
    geom_ribbon(aes(x=as.Date(Date),ymin=(R_Q025_raw),ymax=(R_Q975_raw)),fill="#0000FF40")+
    geom_line(aes(x=as.Date(Date),y=(R_mean_raw)),color="blue")+
    geom_hline(yintercept = 1,linetype=2)+
    scale_x_date(limits=c(as.Date("2021-01-01"),as.Date("2022-10-01")),minor_breaks = as.Date("2020-12-01")+months(1:24))+
    geom_rect(data=data.frame(
      firsts=estiDF[estiDF$firstAbove,"date"]-3,
      lasts=estiDF[estiDF$lastAbove,"date"]+3),
      aes(xmin=firsts,xmax=lasts,ymin=0.5,ymax=1.6),fill="#00000020",color="#00000000"
    )+
    scale_y_continuous(limits=c(0.5,1.6))+
    theme_bw()+
    #theme(axis.line = element_rect(colour = "black"))+
    labs(x="Date",y=expression(R_e))
}

incidencePlot<-function(RTdata, estiDF,threshold=0.3,thresholdOnMax=T){
  estiDF$aboveThreshold<-estiDF$maxLk>threshold
  if(!thresholdOnMax) estiDF$aboveThreshold<-estiDF$meanLk>threshold
  estiDF$lastAbove<-(lead(estiDF$aboveThreshold)!=estiDF$aboveThreshold)&(estiDF$aboveThreshold)
  estiDF$firstAbove<-(lag(estiDF$aboveThreshold)!=estiDF$aboveThreshold)&(estiDF$aboveThreshold)
  estiDF$date<-fConvWtoDate(estiDF$weeksince)
  
  ggplot(RTdata)+
    geom_line(aes(x=as.Date(Date),y=(epicurveRaw)),color="darkgreen")+
    scale_x_date(limits=c(as.Date("2021-01-01"),as.Date("2022-10-01")),minor_breaks = as.Date("2020-12-01")+months(1:24))+
    geom_rect(data=data.frame(
      firsts=estiDF[estiDF$firstAbove,"date"]-3,
      lasts=estiDF[estiDF$lastAbove,"date"]+3),
      aes(xmin=firsts,xmax=lasts,ymin=1,ymax=50000),fill="#00000020",color="#00000000"
    )+
    theme_bw()+
    labs(x="Date",y="reported cases")
}
