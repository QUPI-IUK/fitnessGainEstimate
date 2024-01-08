# Load the source code and data ####
source("dependencies.R")
source("dataLoading.R")
source("SNPlogitDiff.R")
source("lineagePlots.r")
source("logKmethod.r")
source("BuildOwnDb.R")
source("sequenceScoring.r")
source("pangoEstimates.r")

#convertMetaTSVtoRDS("/media/Data/CountryData/Genomic_Data/metadata.tsv",blockSize=1000000)
sREFch<-readRDS("Resources/sREFch.rds")
#TODO: load the correct Rt data


## Build the databases for each geographical area ####
#### Baden-Württemberg ####
system.time(
  BWdata<-buildSNPDF(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_BW_full.fasta",
                   refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                   maskingFilename="Resources/masking_deMaio.tsv"
)
)
saveRDS(BWdata,"BWDataSNPs.rds")

system.time(
  BWdata<-buildSNPDF_blocks(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_BW_full.fasta",
                            refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                            maskingFilename="Resources/masking_deMaio.tsv", outputPrefix="BW2",blockSize = 5000,cores=8
  )
)


#### Germany ####
DEdata<-buildSNPDF(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_DE_full_try.fasta",
                   refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                   maskingFilename="Resources/masking_deMaio.tsv"
)
saveRDS(DEdata,"DEDataSNPs.rds")

DEdata<-buildSNPDF_blocks(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_DE_full_try.fasta",
                   refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                   maskingFilename="Resources/masking_deMaio.tsv", outputPrefix="DE",blockSize = 5000,cores=8
)

#### Denmark ####
DKdata<-buildSNPDF_blocks(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_DK_full_try.fasta",
                          refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                          maskingFilename="Resources/masking_deMaio.tsv", outputPrefix="DK",blockSize = 5000,cores=8
)

#### South Africa ####
system.time(
  RSAdata<-buildSNPDF_blocks(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_RSA_full_try.fasta",
                            refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                            maskingFilename="Resources/masking_deMaio.tsv", outputPrefix="RSA",blockSize = 5000,cores=8
  )
)

#### United States ####
system.time(
  USAdata<-buildSNPDF_blocks(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_USA_full_try.fasta",
                             refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                             maskingFilename="Resources/masking_deMaio.tsv", outputPrefix="USA",blockSize = 5000, cores=8
  )
)

#### United Kingdom ####
system.time(
  UKdata<-buildSNPDF_blocks(msaFilename="/media/Data/CountryData/Genomic_Data/msa_1101/MSA_UK_full.fasta",
                            refFilename="/media/Data/CountryData/Genomic_Data/msa_1101/head.fasta",
                            maskingFilename="Resources/masking_deMaio.tsv", outputPrefix="UK"
  )
)



# Load the sequence meta data (incl. Pango-lineage assignment)

# The main analysis #####
## Loading the German data ####
DEsnps<-readRDS("DESNPs-blocked.rds")
DEmeta<-readRDS("DEmetaData-blocked.rds")
DEfullSeqMetaData<-loadSeqMetaData("Germany")
allPropsDE<-calculateAllProportions(DEsnps)  
RTdata<-read.csv("Resources/Exported_Rt_2023-04-14.csv")

## General statistics on the database ####
# Number of unique position/base combinations
perPosBase<-
  DEsnps%>%group_by(pos,base)%>%
  summarise(orig=sREFch[pos[[1]]],
            #sum=sum(nPosBase),
            n=n()
  )%>%as.data.frame()
nrow(perPosBase)

# Number of changed positions:
perPos<-DEsnps%>%group_by(pos)%>%
  count(base)%>%pivot_wider(names_from = base, values_from = n)%>%
  mutate(nBaseMut = 4-(is.na(A)+is.na(T)+is.na(G)+is.na(C)),
         nTotal = sum(c(A,T,C,G),na.rm=T)
         )%>%
  ungroup()%>%
  arrange(desc(nTotal))%>%
  mutate(index=1:n())  
table(perPos$nBaseMut)
nrow(perPos)

refPositions<-sum(sREFch%in% c("A","T","C","G"))
insertPositions<-sum(perPosBase$orig=="-")

genomeLength<-(refPositions+insertPositions)
changedPosition<-nrow(perPos)
refPositions
insertPositions
genomeLength
# The actual stats
#number of sequences:
nrow(DEmeta)
#number of positions in reference
refPositions
#number of insertions
insertPositions
#percentage of positions with changes
nrow(perPos)/genomeLength
nrow(perPos)
table(perPos$orig)

posBaseCount<-snpFreqs%>%group_by(posBase)%>%
  summarise(nTotal=sum(nPosBase)
  )%>%arrange(desc(nTotal))%>%as.data.frame()
posBaseCount$index<-1:nrow(posBaseCount)

#number of unique SNPs:
nrow(perPosBase)
sum(perPos$nBaseMut==1)
sum(perPos$nBaseMut==1)/ nrow(perPosBase)
(nrow(DEmeta)/2)
#number of singeltons
sum(perPosBase$n==1)/ nrow(perPosBase)

#sum(perPos$nTotal10000)/nrow(DEmeta)

#number of SNPs present in over 50% of the sequence
sum(perPosBase$n>=(nrow(DEmeta)/2))

100*sum(perPosBase$n>=(nrow(DEmeta)/2))/ nrow(DEmeta)

ggplot(perPos)+
  geom_point(aes(x=index,y=nTotal/nrow(DEmeta)))+
  scale_x_log10()+
  scale_y_log10(labels = scales::percent_format())+
  theme_bw()+
  labs(x="SNP index",y="Percentage of sequences",title="Presence of SNPs in isolates")
ggsave("SNPpresence.svg")

DEmeta[1,]
min(as.Date(DEmeta$date),na.rm=T)
max(as.Date(DEmeta$date),na.rm=T)
### Figure S1: Plot of isolates per week #####
weeklyIsolates<-as.data.frame((table(DEmeta$weeksince)))%>%
  mutate(date=fConvWtoDate(as.numeric(as.character(Var1))))

ggplot(weeklyIsolates)+
  geom_ribbon(aes(x=date,ymin=0,ymax=Freq),fill="#0000FF90")+
  theme_bw()+
  labs(x="Isolation date",y="Number of sequences")
ggsave("isolatesPerWeek-DE.svg")

## SNP based estimate ####
### Figure 1: Space invader plots ####
#filtering parameters:
filterAlpha<-0.01/(genomeLength*146*4)
filterAlphaNonConserv<-0.0000000001/nrow(allPropsDE)
spaceInvaderPlot(allPropsDE,102,103,examplePosBase="11836Reference",filterAlpha = filterAlphaNonConserv)
ggsave("SpaceInvaderExample-DE.svg")

snpProportionTimeline(allPropsDE,96,118,102,examplePosBase="11836Reference")+theme(aspect.ratio = 1)
ggsave("snpProportionTimeline-DE.svg")

L_overTimePlot(allPropsDE,0,146,filterAlpha = filterAlphaNonConserv)
ggsave("L_overTimePlot-DE.svg")
ggsave("L_overTimePlot-DE.png")

## log K over time ####
### Figure 2: Theory on log(k) ####

theoryPlot(p=c(0.5,0.25,0.01))
ggsave("theoryplot.svg")

theoryPlotVar()
ggsave("theoryplotWithVar.svg")

timelineExample<-calculateLogKsSingle(calculateAllProportions(DEsnps),focalWeek = 106)
timelineStiExamplePlot(timelineExample)
ggsave("exampleEstiTL-1.svg")
fConvWtoDate(106)

timelineExample<-calculateLogKsSingle(calculateAllProportions(DEsnps),focalWeek = 118)
timelineStiExamplePlot(timelineExample)
ggsave("exampleEstiTL-2.svg")
fConvWtoDate(118)

## Delta Log(k) #####
# The MAIN thing: 
system.time(
  dLk_DE<-estimDeltaLogK(DEsnps)
)
increasedPeriods<-extractPeriods(dLk_DE[dLk_DE$weeksince>50,],threshold=0.125)
saveRDS(dLk_DE,"delta_log_K_DE.rds")

increasedPeriods

#### Figure 3: delta log k timeline ####

#trial of the plots:
  # Fitness over time #
    fitnessPlot(dLk_DE[dLk_DE$weeksince>50,],startTF=(-8),endTF=(-2),serialInterval = 5,threshold=0.125,thresholdOnMax=T)
  # Lineages over time #
    lineagePerWeekPlot(DEfullSeqMetaData,dLk_DE,threshold=0.125)
  # Effective R #
    effectiveRplot(RTdata,dLk_DE,threshold=0.125)
  # Incidence #
    incidencePlot(RTdata,dLk_DE,threshold=0.125)

  #Combined figure
    g1<-ggplotGrob(fitnessPlot(dLk_DE[dLk_DE$weeksince>50,],startTF=(-8),endTF=(-2),serialInterval = 5,threshold=0.125,thresholdOnMax=T))
    g2<-ggplotGrob(lineagePerWeekPlot(DEfullSeqMetaData,dLk_DE,threshold=0.125))
    g3<-ggplotGrob(effectiveRplot(RTdata,dLk_DE,threshold=0.125))
    g4<-ggplotGrob(incidencePlot(RTdata,dLk_DE,threshold=0.125))

    gCombi <- rbind(rbind(rbind(g1, g2, size = "first"), g3,size = "first"),g4,size = "first")

    gCombi$widths <-unit.pmax(g1$widths, g2$widths)
    grid.newpage()
    grid.draw(gCombi)

    ggsave("DE-FitnessAndVariantPeriods.svg",gCombi,width=8,height=12)
    ggsave("DE-FitnessAndVariantPeriods.pdf",gCombi,width=8,height=12)
    
    increasedPeriods
    dLk_DE

    colorSchemeLineages<-createLineageColorscheme(DEfullSeqMetaData)
    increasedPeriods
    DEfullSeqMetaData%>%filter(weeksince==55)%>%
      group_by(alpha)%>%
      summarise(n())
    1851/(1851+1529)
    #
## Score isolates for aligning with increasing SNPs for each time frame ####

#a single period:
scoresP3<-getPosNegLineages(increasedPeriods[3,"firstWk"]-1,
                  increasedPeriods[3,"lastWk"],
                  DEsnps,
                  allPropsDE,
                  DEmeta,
                  DEfullSeqMetaData,
                  scoreCut=0)
getScorePlot(scoresP3$scoredSeqs,colorSchemeLineages,increasedPeriods,3)

seqScoresAllPeriods<-lapply(1:nrow(increasedPeriods),
                         function(x){
                           print(x)
                           getPosNegLineages(increasedPeriods[x,"firstWk"]-1,
                             increasedPeriods[x,"lastWk"],
                             DEsnps,
                             allPropsDE,
                             DEmeta,
                             DEfullSeqMetaData,
                             scoreCut=0)})

#try one plot or table

getScorePlot(3,seqScoresAllPeriods,increasedPeriods,colorSchemeLineages)
get2x2Table(5,seqScoresAllPeriods,increasedPeriods,colorSchemeLineages,showN=1)

getScorePlotList(seqScoresAllPeriods,colorSchemeLineages,increasedPeriods,2)

#produce all
allScorePlots<-lapply(1:nrow(increasedPeriods),function(x)getScorePlotList(seqScoresAllPeriods,colorSchemeLineages,increasedPeriods,x))

ggarrange(plotlist=allScorePlots[2:5])
ggsave("SeqScorePeriods-DE.svg",width=10,height=8)
ggsave("SeqScorePeriods-DE.pdf",width=10,height=8)

# Sensitivity analysis ####
## Time of detection #####
gc()
system.time({
  timeSteps<-
    lapply(50:146,function(w){
      gc()
      print(paste0("Week ",w))
      dLk_temp<-estimDeltaLogK(getSeqsBeforeSubmWeek(DEsnps,DEmeta,DEfullSeqMetaData,w),splitLogK = T)
      dLk_temp$week<-w
      saveRDS(dLk_temp,paste0("dLk_temp",w,".RDS"))
      return(dLk_temp)
    })%>% bind_rows()
})

timeSteps<-lapply(56:146,function(w)
       readRDS(paste0("dLk_temp",w,".RDS"))
       )%>% bind_rows()

saveRDS(timeSteps,"timeSteps-DE.rds")
timeSteps[timeSteps$weeksince==timeSteps$week-3,]

timeSteps[timeSteps$week==102,]
fConvWtoDate(101)

ggplot(timeSteps[timeSteps$week==103,])+
  geom_step(aes(x=weeksinceToDate(weeksince),y=maxLk),color="red")+
  geom_hline(yintercept = 0.125,linetype=2)+
  theme_bw()

ggplot(timeSteps[timeSteps$weeksince==(timeSteps$week-2),])+
  geom_step(aes(x=weeksinceToDate(week),y=maxLk),color="red")+
  geom_hline(yintercept = 0.125,linetype=2)+
  theme_bw()+
  scale_y_continuous(limits=c(0,2))

extractPeriods(timeSteps[timeSteps$weeksince==(timeSteps$week-2),]%>%mutate(weeksince=week))

ggarrange(
  ggplot(timeSteps[timeSteps$weeksince==(timeSteps$week),])+
    geom_step(aes(x=weeksinceToDate(week),y=maxLk),color="red")+
    geom_hline(yintercept = 0.125,linetype=2)+
    theme_bw()+
    scale_y_continuous(limits=c(0,1.5))+
    labs(x="",y=expression(Delta*log(k)),title="Compared to this week"), 
  ggplot(timeSteps[timeSteps$weeksince==(timeSteps$week-1),])+
    geom_step(aes(x=weeksinceToDate(week),y=maxLk),color="red")+
    geom_hline(yintercept = 0.125,linetype=2)+
    theme_bw()+
    scale_y_continuous(limits=c(0,1.5))+
    labs(x="",y=expression(Delta*log(k)),title="Compared to one week before week"),
  ggplot(timeSteps[timeSteps$weeksince==(timeSteps$week-2),])+
    geom_step(aes(x=weeksinceToDate(week),y=maxLk),color="red")+
    geom_hline(yintercept = 0.125,linetype=2)+
    theme_bw()+
    scale_y_continuous(limits=c(0,1.5))+
    labs(x="",y=expression(Delta*log(k)),title="Compared to two weeks before week"),
  ggplot(timeSteps[timeSteps$weeksince==(timeSteps$week-3),])+
    geom_step(aes(x=weeksinceToDate(week),y=maxLk),color="red")+
    geom_hline(yintercept = 0.125,linetype=2)+
    theme_bw()+
    scale_y_continuous(limits=c(0,1.5))+
    labs(x="",y=expression(Delta*log(k)),title="Compared to three weeks before week")
  )
ggsave("LiveDetecion.svg")
ggsave("LiveDetecion.png")

## Resampling (Minimum number of isolates needed) #####

#Loop for resampling analysis
resample10<-resampledLK(DEsnps,10,verbose=T,repeats = 100)
resample25<-resampledLK(DEsnps,25,verbose=T,repeats=100)
resample50<-resampledLK(DEsnps,50,verbose=T,repeats=100)
resample100<-resampledLK(DEsnps,100,verbose=T,repeats=100)
resample250<-resampledLK(DEsnps,250,verbose=T,repeats=100)
resample500<-resampledLK(DEsnps,500,verbose=T,repeats=100)

#resample750<-resampledLK(DEsnps,750,verbose=T,repeats=10)
resample1000<-resampledLK(DEsnps,1000,verbose=T,repeats=10)

resample25[resample25$loop==1,]

timingSizes<-bind_rows(list(resample10,resample25,resample50,resample100,resample500,resample750,resample250,resample1000))%>%
  group_by(loop,sampleSize)%>%
  summarise(timing=timing[1]
            )

ggplot(timingSizes)+
  geom_point(aes(x=log(sampleSize),y=(timing)))

resample50%>%group_by(loop)%>%summarise(timing=timing[1])
resample100%>%group_by(loop)%>%summarise(timing=timing[1])
resample250%>%group_by(loop)%>%summarise(timing=timing[1])

estiBounds10<-tlLkQunatiles(resample10)
estiBounds25<-tlLkQunatiles(resample25)
estiBounds50<-tlLkQunatiles(resample50)
estiBounds100<-tlLkQunatiles(resample100)
estiBounds250<-tlLkQunatiles(resample250)
estiBounds500<-tlLkQunatiles(resample500)

ggplot(resample500)+
  geom_line(aes(x=weeksinceToDate(weeksince),y=maxLk, group=loop),color="#00000020")+
  geom_hline(yintercept = 0.125,linetype=2)+  
  theme_bw()
            
ggplot(estiBounds25)+
  geom_ribbon(aes(x=weeksinceToDate(weeksince),ymin=pmax(lowerRange,0),ymax=pmin(upperRange,3)),fill="#FF000050")+
  geom_ribbon(aes(x=weeksinceToDate(weeksince),ymin=pmax(lowerIQR,0),ymax=pmin(upperIQR,3)),fill="#FF000050")+
  geom_line(aes(x=weeksinceToDate(weeksince),y=median),color="red")+
  geom_hline(yintercept = 0.125,linetype=2)+
  scale_x_date(limit=c(as.Date("2021-01-01"),as.Date("2022-11-01")))+
  scale_y_continuous(limits=c(0,3))+
  theme_bw()
ggsave("Resampled25.svg")  


## time scaling of algorithm, based on resampling #####
allResampleData<-bind_rows(list(resample10,
     resample25,
     resample50,
     resample100,
     resample250,
     resample500))

resampleSpeed<-allResampleData%>%group_by(sampleSize,loop)%>%
  summarise(timing=mean(timing))

nrow(resample25[resample25$loop==1,])

ggplot(resampleSpeed)+
  geom_point(aes(x=sampleSize,y=timing),color="#00000020")+
  labs(x="Sample size (sequences per week)",y="Time (seconds)")+
  theme_bw()
ggsave("TimeScaling.svg")

#Add the countries (based on individual measurements):
countryTimings<-data.frame(name=c("US","UK","Denmark","Germany","South Africa"),
           time=c(128,82.7,38.8,60.4,10.5),
           samples=c(3769411,2774418,576296,797206,35770)/141
           )

ggplot(resampleSpeed)+
  geom_point(aes(x=sampleSize,y=timing),color="#00000020")+
  geom_point(data=countryTimings,aes(x=samples,y=time),shape=21,size=3,color="red",fill="yellow")+
  geom_point(data=countryTimings,aes(x=samples,y=time),shape=3,color="red")+
  labs(x="Sample size (sequences per week)",y="Time (seconds)")+
  scale_x_log10()+
  scale_y_log10()+
  theme_bw()
ggsave("TimeScaling-inclCountries.svg")


## Different lag periods ####
allProps<-calculateAllProportions(DEsnps)
logKdf<-lapply(1:51,function(l){
  print(paste0("lag ",l))
  calculateLogKs(allProps,lagVec=c(l))
})%>%bind_rows()
alldKs<-lapply(1:51,function(l){    
  logKdf%>%filter(lagT<=l)%>%
    group_by(weeksince)%>%
  summarise(meanLk=mean(diff(logK,lag=1)),
            maxLk=max(diff(logK,lag=1)),
            meanHp=mean(sumHp),
            meanHs=mean(sumHs)
  )%>%mutate(includeL=l,
             date=fConvWtoDate(weeksince)
             )%>%
  as.data.frame()}
)%>%bind_rows()

ggplot(alldKs%>%filter(includeL==7))+
  geom_line(aes(x=date,y=maxLk,group=includeL,color=includeL))+
  scale_x_date(limits=c(as.Date("2021-01-01"),as.Date("2022-11-01")))+
  scale_y_continuous(limits=c(-0.25,1))+
  geom_hline(yintercept = 0.125,linetype=2)+
  theme_bw()+theme(legend.position = "bottom")+
  labs(x="Date",y=expression(Delta*log(k)),color="Included timeframe (weeks previous)")

ggplot(alldKs%>%filter(includeL<=14))+
  geom_line(aes(x=date,y=maxLk,group=includeL,color=includeL))+
  scale_x_date(limits=c(as.Date("2021-01-01"),as.Date("2022-11-01")))+
  scale_y_continuous(limits=c(-0.25,1))+
  geom_hline(yintercept = 0.125,linetype=2)+
  theme_bw()+theme(legend.position = "bottom")+
  labs(x="Date",y=expression(Delta*log(k)),color="Included timeframe (weeks previous)")
ggsave("LagT-diffs.svg")

ggplot(alldKs%>%filter(includeL<=14))+
  geom_line(aes(x=date,y=meanLk,group=includeL,color=includeL))+
  scale_x_date(limits=c(as.Date("2021-01-01"),as.Date("2022-11-01")))+
  scale_y_continuous(limits=c(-0.25,1))+
  geom_hline(yintercept = 0.125,linetype=2)+
  theme_bw()+theme(legend.position = "bottom")+
  labs(x="Date",y=expression(Delta*log(k)),color="Included timeframe (weeks previous)")


## Other geographies ####

gc()
USsnps<-readRDS("USASNPs-blocked.rds")
max(USsnps$seqID)
USmeta<-readRDS("USAmetaData-blocked.rds")
system.time(
  dLk_USA<-estimDeltaLogK(USsnps,splitLogK = T)
)
max(dLk_USA$maxLk[!(dLk_USA$maxLk%in%c(Inf,-Inf,NA,NaN))])*5.5/7
fitnessPlot(dLk_USA[dLk_USA$weeksince>55,],startTF=(-8),endTF=(-2),serialInterval = 5,threshold=0.125,thresholdOnMax=T)
ggsave("UStimeline.svg")
saveRDS(dLk_USA,"US_dLk.rds")

rm(USsnps)

#### Denmark ####
DKsnps<-readRDS("DKSNPs-blocked.rds")
DKmeta<-readRDS("DKmetaData-blocked.rds")
max(DKsnps$seqID)
system.time(
  dLk_DK<-estimDeltaLogK(DKsnps)
)
max(dLk_DK$maxLk[!(dLk_DK$maxLk%in%c(Inf,-Inf,NA,NaN))])*5.5/7
fitnessPlot(dLk_DK[dLk_DK$weeksince>55,],startTF=(-8),endTF=(-2),serialInterval = 5,threshold=0.125,thresholdOnMax=T)
ggsave("DKtimeline.svg")
saveRDS(dLk_DK,"DK_dLk.rds")

rm(DKsnps)

#### United Kingdom ####
gc()
UKsnps<-readRDS("UKSNPs-blocked.rds")
UKmeta<-readRDS("UKmetaData-blocked.rds")
fullSeqMetaDataUK<-loadSeqMetaData(c("England","Wales","Scotland","Northern Ireland"))

max(UKsnps$seqID)
system.time(
  dLk_UK<-estimDeltaLogK(UKsnps,splitLogK = T)
)
fitnessPlot(dLk_UK[dLk_UK$weeksince>55,],startTF=(-8),endTF=(-2),serialInterval = 5,threshold=0.125,thresholdOnMax=T)
(dLk_UK$maxLk[!(dLk_UK$maxLk%in%c(Inf,-Inf,NA,NaN))])
ggsave("UKtimeline.svg")
saveRDS(dLk_UK,"UK_dLk.rds")

rm(UKsnps)
#### South Africa #####
RSAsnps<-readRDS("RSASNPs-blocked.rds")
fullSeqMetaDataRSA<-loadSeqMetaData("South Africa")
RSAmeta<-readRDS("RSAmetaData-blocked.rds")

max(RSAsnps$seqID)
system.time(
  dLk_RSA<-estimDeltaLogK(RSAsnps)
)
max(dLk_RSA$maxLk[!(dLk_RSA$maxLk%in%c(Inf,-Inf,NA,NaN))])*5.5/7

increasedPeriods<-extractPeriods(dLk_RSA[dLk_RSA$weeksince>55,],threshold=0.125)
fitnessPlot(dLk_RSA[dLk_RSA$weeksince>55,],startTF=(-8),endTF=(-2),serialInterval = 5,threshold=0.125,thresholdOnMax=T)#+theme(aspect.ratio = 1/2)
lineagePerWeekPlot(fullSeqMetaDataRSA,dLk_RSA,threshold=0.25)
ggsave("RSAtimeline.svg")

#### Baden-Württemberg #####
BWsnps<-readRDS("BW2SNPs-blocked.rds")
BWmeta<-readRDS("BW2metaData-blocked.rds")
system.time(
  dLk_BW<-estimDeltaLogK(BWsnps)
)
#fitnessBW<-estimDeltaLogK(BWsnps)
increasedPeriods<-extractPeriods(dLk_BW[dLk_BW$weeksince>55,],threshold=0.25)
fitnessPlot(dLk_BW[dLk_BW$weeksince>55,],
            startTF=(-8),endTF=(-2),serialInterval = 5.5,
            threshold=0.25,
            thresholdOnMax=T)
ggsave("BW-timeline.svg")

lineagePerWeekPlot(fullSeqMetaDataBW,dLk_BW,threshold=0.25)


# Export the used sequence labels for GISAID identifier ####


write.csv(DEmeta[,c("seqLabel","EPI")],"UsedSet-GISAID-DE.csv")
write.csv(DKmeta[,c("seqLabel","EPI")],"UsedSet-GISAID-DK.csv")
write.csv(RSAmeta[,c("seqLabel","EPI")],"UsedSet-GISAID-RSA.csv")
write.csv(UKmeta[,c("seqLabel","EPI")],"UsedSet-GISAID-UK.csv")
write.csv(USmeta[,c("seqLabel","EPI")],"UsedSet-GISAID-US.csv")
write.csv(BWmeta[,c("seqLabel","EPI")],"UsedSet-GISAID-BW.csv")

