px<-function(t,px0,a){
  1/(1+(((1-px0)/px0)*exp(-r*t)))
}
test<-function(x,C1,a){1 / (1+ ((1-C1)/C1)*exp(-r*x))}


getVOCcounts<-function(data,week){
subData<-data%>%filter(weeksinceSubmit<=week)%>%
  group_by(weeksince)%>%
  summarise(nsample=n(),
            nVOC=sum(Omicron1),
            propVOC=nVOC/n(),
            week=(weeksince[[1]]),
            date=weeksinceToDate(weeksince[[1]]),
            cutoff=weeksinceToDate(week)
  )%>%as.data.frame()
return(subData[!is.na(subData$week),])
}

getVOCfit<-function(data,week,incSub){
  VOCcounts<-getVOCcounts(data,week)
  VOCcounts$time<-(VOCcounts$week)*7
  startProp<-(VOCcounts[VOCcounts$time>0,][1,"propVOC"])
  startProp<-VOCcounts[VOCcounts$propVOC>0.01,][1,"propVOC"]
  VOCcounts$time<-VOCcounts$time-(VOCcounts[VOCcounts$propVOC>0.01,][1,"time"])
  print(startProp)
  finalResDF<-tryCatch({
    doublefit<-nls(nVOC ~ nsample / (1+ (((1-C1)/C1)*exp(-a*time))),data=VOCcounts,start=list(C1=startProp,a=0.74))
    coefs<-as.data.frame(summary(doublefit)[10]$coefficients)
    resDF<-data.frame(
      date=week,
      #refDate=refDate,
      totVOC=sum(VOCcounts$nVOC),
      pointa=coefs["a","Estimate"],
      lowera=coefs["a","Estimate"]-1.96*coefs["a","Std. Error"],
      uppera=coefs["a","Estimate"]+1.96*coefs["a","Std. Error"],
      pointC=coefs["C1","Estimate"],
      lowerC=coefs["C1","Estimate"]-1.96*coefs["C1","Std. Error"],
      upperC=coefs["C1","Estimate"]+1.96*coefs["C1","Std. Error"]
    )
    return(list(resDF=resDF,VOCcounts=VOCcounts))
  },
  error=function(cond) {
    resDF<-data.frame(
      date=week,
     #refDate=refDate,
      totVOC=sum(VOCcounts$nVOC),
      pointa=NA,
      lowera=NA,
      uppera=NA,
      pointC=NA,
      lowerC=NA,
      upperC=NA
    )
    return(list(resDF=resDF,VOCcounts=VOCcounts))
    
  }
  )
  #eturn(Reduce(function(l){rbind(l$resDF)},finalResDF))
  return(finalResDF)
}


getVOCfit1<-function(VOCcounts,dum){
  VOCcounts<-as.data.frame(VOCcounts)
   # VOCcounts<-getVOCcounts(data,week)
  VOCcounts$time1<-(VOCcounts$week)*7
#  startProp<-(VOCcounts[VOCcounts$time>0,][1,"propVOC"])
  startProp<-VOCcounts[VOCcounts$propVOC>0.01,][1,"propVOC"]
#  print(VOCcounts[VOCcounts$propVOC>0.01,][1,"time1"])

  VOCcounts$time<-(VOCcounts$time1)-(VOCcounts[VOCcounts$propVOC>0.01,][1,"time1"])
  #print(startProp)
  #print(VOCcounts)
  finalResDF<-tryCatch({
#    doublefit<-glm(nVOC ~ nsample / (1+ (((1-C1)/C1)*exp(-a*time))),data=VOCcounts,start=list(C1=startProp,a=0.14))
 #   print(doublefit)
    doublefit<-nls(nVOC ~ nsample / (1+ (((1-C1)/C1)*exp(-a*time))),data=VOCcounts,start=list(C1=startProp,a=0.14))
    
    coefs<-as.data.frame(summary(doublefit)[10]$coefficients)
    resDF<-data.frame(
  #    date=week,
      #refDate=refDate,
      startProp=startProp,
      totSample=sum(VOCcounts$nsample),
      totVOC=sum(VOCcounts$nVOC),
      pointa=coefs["a","Estimate"],
      lowera=coefs["a","Estimate"]-1.96*coefs["a","Std. Error"],
      uppera=coefs["a","Estimate"]+1.96*coefs["a","Std. Error"],
      pointC=coefs["C1","Estimate"],
      lowerC=coefs["C1","Estimate"]-1.96*coefs["C1","Std. Error"],
      upperC=coefs["C1","Estimate"]+1.96*coefs["C1","Std. Error"]
    )
    return(resDF)
    #return(list(resDF=resDF,VOCcounts=VOCcounts))
  },
  error=function(cond) {
    resDF<-data.frame(
     # date=week,
      #refDate=refDate,
      startProp=startProp,
      totSample=sum(VOCcounts$nsample),
      totVOC=sum(VOCcounts$nVOC),
      pointa=NA,
      lowera=NA,
      uppera=NA,
      pointC=NA,
      lowerC=NA,
      upperC=NA
    )
    return(resDF)
    #return(list(resDF=resDF,VOCcounts=VOCcounts))
    
  }
  )
  #eturn(Reduce(function(l){rbind(l$resDF)},finalResDF))
  return(finalResDF)
}

# Calculate the variant proportions per week.
getVarProps<-function(seqMetaData){
  seqMetaData%>%
  group_by(weeksince)%>%
  summarise(
    n=n(),
    propAlpha=sum(alpha)/n(),
    propDelta=sum(Delta)/n(),
    propOmicron1=sum(Omicron1)/n(),
    propOmicron2=sum(Omicron2)/n(),
    propOmicron3=sum(Omicron3)/n(),
    propOmicron5=sum(Omicron5)/n()
  )%>%
  mutate(
    alphaGrowth=(propAlpha>0.05&propAlpha<lead(propAlpha,4)),
    deltaGrowth=(propDelta>0.05&propDelta<lead(propDelta,4)),
    omicron1Growth=(propOmicron1>0.05&propOmicron1<lead(propOmicron1,4)),
    omicron2Growth=(propOmicron2>0.05&propOmicron2<lead(propOmicron2,4)),
    omicron3Growth=(propOmicron2>0.05&propOmicron2<lead(propOmicron3,4)),
    omicron5Growth=(propOmicron2>0.05&propOmicron2<lead(propOmicron5,4)),
  )
}