getSeqSNPlocs<-function(seq,refChars,minPos=374,maxPos=50386){
  seqChars<-unlist(strsplit(as.character(seq), split = ""))
  locs<-(which(refChars!=seqChars))
  locs<-locs[locs>=minPos&locs<=maxPos]
  df<-data.frame(
    pos=locs,
    base=unlist(seqChars[locs]),
    origbase=unlist(refChars[locs])
  )
  #Note: insertions and deletions are included:
  #df<-df[df$base!="-",]
  #df<-df[df$origbase!="-",]
  
  #df<-df[df$base!="N",]
  #df<-df[df$origbase!="N",]
  df<-df[!(df$base %in% c("N","R","Y","K","M","S","W","B","D","H","V","-")),]  # exclude deletions
  df<-df[!(df$origbase %in% c("N","R","Y","K","M","S","W","B","D","H","V")),]
  if(nrow(df)>0)  rownames(df)<-(1:nrow(df))
  return(df)
}

#lapply(sBW[12:15],getSeqSNPlocs)
#Database consists of:
#metadata
# -tag
# -country
# -seqID (generated)
# -seqLabel
# -year
# -EPI
# -date
# -region
# -nSNPS
 
#seqdata
# - seqID
# - pos
# - base

logit<-function(p) log(p/(1-p))

buildSNPDF<-function(msaFilename,refFilename,maskingFilename="",firstNseq=-1){
  print("Opening reference data")
  sREF = readDNAStringSet(refFilename)[1]
  print("Opening sequence data file")
  sBW = readDNAStringSet(msaFilename)
  
  print("Converting reference sequence")
  refChars<-unlist(strsplit(as.character(sREF), split = ""))

  #### general front/end masking ####
  origReferencePos<-which(refChars!="-")
  #1-55
  #29804-29903
  firstExcl<-(origReferencePos[55])
  lastExcl<-(origReferencePos[29804])
  print(paste0("30th ref position in MSA is: ",firstExcl) )
  print(paste0("-30th ref position in MSA is: ",(lastExcl)) )
  print(as.character(refChars[(1:firstExcl)]))
  print(as.character(refChars[lastExcl:length(refChars)]))
  refChars[(1:firstExcl)]<-"N"
  refChars[lastExcl:length(refChars)]<-"N"
  #browser()
  
  ########################
  
  if(maskingFilename!=""){
    masking<-read.csv(maskingFilename,sep = "\t")
    origReferencePos<-which(refChars!="-")
    maskedAndCaution<-data.frame(
      refPos=masking$POS,
      msaPos=as.numeric(origReferencePos[masking$POS]),
      base=masking$REF,
      msaBase=as.character(refChars[origReferencePos[masking$POS]])
    )
    print("Excluding masked positions (By entering Ns)")
    refChars[maskedAndCaution$msaPos]<-"N"
  }
  if(firstNseq==-1){
    buildNSeqs<-length(sBW)
  } else {
    buildNSeqs<-firstNseq
  }
  print(paste0("Creating the dataframe, from ",buildNSeqs," records..."))
  print(paste0("Extracting the SNPs..."))
  #pb <- progress_bar$new(total = length(sBW)-104485)
  pb <- progress_bar$new(total = buildNSeqs,format = "  :current/:total [:bar] :percent eta: :eta")
  #for (i in 1:buildNSeqs) {
  snpDF<-lapply(1:buildNSeqs, function(i){
    thisSeqID<-i
    pb$tick()
    seqDF<-getSeqSNPlocs2(sBW[i],refChars)
    if(nrow(seqDF)>0) seqDF$seqID<-thisSeqID
    return(seqDF)
  })
  print("Extracting the metadata...")
  pb2 <- progress_bar$new(total = buildNSeqs,format = "  :current/:total [:bar] :percent eta: :eta")
  #for (i in 1:buildNSeqs) {
  metaDF<-lapply(1:buildNSeqs, function(i){
    pb2$tick()
    thisSeqID<-i
    thisTag<-as.character(names(sBW[i]))[[1]]
    alltags<-unlist(str_split(thisTag,fixed("|")))
    allSubTags<-unlist(str_split(alltags[1],fixed("/")))
    thisDate<-as.Date(alltags[3],format="%Y-%m-%d")

    metaDF<-data.frame(tag=thisTag,
                       country=allSubTags[2],
                       seqID=thisSeqID,
                       seqLabel=allSubTags[3],
                       year=as.numeric(allSubTags[4]),
                       week=paste0(isoyear(thisDate),"-", isoweek(thisDate)),
                       EPI=alltags[2],
                       date=thisDate,
                       region=alltags[4],
                       nSNPS=nrow(snpDF[[i]])
    )
  }) %>%bind_rows()
  print("Done")
  snpDF<-snpDF%>%bind_rows()
  
  metaDF$date<-as.Date(metaDF$date)
  metaDF$weeksince<-as.numeric(floor((metaDF$date-as.Date("2020-01-05"))/7))
  
  weekCounts<-
    metaDF%>%#filter(week=="2021-38")%>%
    mutate(weeksince=floor((as.numeric(date)-18266)/7))%>% #days difference to 5 Jan 2020 (A Sunday)
    select(seqID,weeksince)%>%
    group_by(weeksince)%>%
    mutate(nSeq=n())
  
  snpDF<-snpDF%>%left_join(weekCounts,by="seqID")
  
  return(list(
    SNPs=snpDF,
    meta=metaDF
  ))
}

buildSNPDF_blocks<-function(msaFilename,refFilename,blockSize=25000, maskingFilename="",firstNseq=-1,outputPrefix="TMP",cores=4){
  print("Opening reference data")
  sREF = readDNAStringSet(refFilename)[1]
  print("Opening sequence data file")
  #sBW = readDNAStringSet(msaFilename)
  
  print("Converting reference sequence")
  refChars<-unlist(strsplit(as.character(sREF), split = ""))
  #### general front/end masking ####
  origReferencePos<-which(refChars!="-")
  #1-55
  #29804-29903
  firstExcl<-(origReferencePos[55])
  lastExcl<-(origReferencePos[29804])
  print(paste0("30th ref position in MSA is: ",firstExcl) )
  print(paste0("-30th ref position in MSA is: ",(lastExcl)) )
  print(as.character(refChars[(1:firstExcl)]))
  print(as.character(refChars[lastExcl:length(refChars)]))
  refChars[(1:firstExcl)]<-"N"
  refChars[lastExcl:length(refChars)]<-"N"
  #browser()
  
  ########################
  
  if(maskingFilename!=""){
    masking<-read.csv(maskingFilename,sep = "\t")
    origReferencePos<-which(refChars!="-")
    maskedAndCaution<-data.frame(
      refPos=masking$POS,
      msaPos=as.numeric(origReferencePos[masking$POS]),
      base=masking$REF,
      msaBase=as.character(refChars[origReferencePos[masking$POS]])
    )
    print("Excluding masked positions (By entering Ns)")
    refChars[maskedAndCaution$msaPos]<-"N"
  }
  print("Estimatting number of records...")
  total_records <- as.integer(system2("wc",
                                        args = c("-l",
                                                 msaFilename,
                                                 " | awk '{print $1}'"),
                                        stdout = TRUE))
  
  total_records<-total_records/2
  print(paste0("Found ",total_records," records..."))
  print(paste0("That should result in ",ceiling(total_records/blockSize) ," blocks and files"))
  print("Splitting the file...")
  #split the file
  system(paste("split","-l", format(blockSize*2, scientific=FALSE),# 100000,",#blockSize*2,
                   msaFilename,
                   paste0(dirname(msaFilename),"/tempSplit-Seg")))
  
  if(firstNseq==-1){
    buildNSeqs<-total_records#length(sBW)
  } else {
    buildNSeqs<-firstNseq
  }
  #split into blocks

#  pb <- progress_bar$new(total = buildNSeqs,format = "  :current/:total [:bar] :percent eta: :eta")
 # pb2 <- progress_bar$new(total = buildNSeqs,format = "  :current/:total [:bar] :percent eta: :eta")
  print(paste0("Extracting the SNPs..."))
  nBlock=ceiling(buildNSeqs/blockSize)
  print(paste0("Divided in ",nBlock," blocks"))
  segmentFileList=list.files(dirname(msaFilename),"tempSplit-Seg")
  print(paste0("Found ",length(segmentFileList)," split files"))
  print(segmentFileList)
  if(length(segmentFileList)!=nBlock) stop("number of blocks and files don't match!")
  # loop 1:nBlock
  # each loop handles sequences 
   pbmclapply(1:nBlock, function(bl){
  #pbmclapply(1:2, function(bl){
  # lapply(1:nBlock, function(bl){
   # print(bl)
    beginPointBlock<-((blockSize*(bl-1))+1)
    endPointBlock<-min(blockSize*bl,buildNSeqs)
    
    #sBW = readDNAStringSet(msaFilename,skip=(beginPointBlock-1),nrec=blockSize)
    sBW = readDNAStringSet(paste0(dirname(msaFilename),"/",segmentFileList[bl]))#,skip=(beginPointBlock-1),nrec=blockSize)

    #blockSize
    print(paste0("Running block from ",beginPointBlock," to ",endPointBlock))
    
  #pb <- progress_bar$new(total = length(sBW)-104485)
  #for (i in 1:buildNSeqs) {
  snpDF<-lapply(beginPointBlock:endPointBlock, function(i){
    thisSeqID<-i
    thisSeq<-sBW[i-(beginPointBlock-1)]
  #  if(i%%100==0)pb$tick(100)
    seqDF<-getSeqSNPlocs2(thisSeq,refChars)
    if(nrow(seqDF)>0) seqDF$seqID<-thisSeqID
    return(seqDF)
  })

  print("Extracting the metadata...")
  
  #for (i in 1:buildNSeqs) {
  metaDF<-lapply(beginPointBlock:endPointBlock, function(i){
    #if(i%%100==0)pb2$tick(100)
    thisSeqID<-i
    thisSeq<-sBW[i-(beginPointBlock-1)]
    thisTag<-as.character(names(thisSeq))[[1]]
    alltags<-unlist(str_split(thisTag,fixed("|")))
    allSubTags<-unlist(str_split(alltags[1],fixed("/")))
    thisDate<-as.Date(alltags[3],format="%Y-%m-%d")
    
    data.frame(tag=thisTag,
                       country=allSubTags[2],
                       seqID=thisSeqID,
                       seqLabel=allSubTags[3],
                       year=as.numeric(allSubTags[4]),
                       week=paste0(isoyear(thisDate),"-", isoweek(thisDate)),
                       EPI=alltags[2],
                       date=thisDate,
                       region=alltags[4],
                       nSNPS=nrow(snpDF[[(i-beginPointBlock)+1]])
    )
  }) %>%bind_rows()
  snpDF<-snpDF%>%bind_rows()
  saveRDS(snpDF,paste0("SNPextract-block",bl,".rds"))
  saveRDS(metaDF,paste0("MetaDFextract-block",bl,".rds"))
  },mc.cores=cores
  )
  print("Done")

  #Combine the blocks:
  tempDump<-lapply(1:nBlock,
                   function(bl){
                     readRDS(paste0("SNPextract-block",bl,".rds"))
                   }
  ) %>% bind_rows()
  tempMetaDump<-lapply(1:nBlock,
                   function(bl){
                     readRDS(paste0("MetaDFextract-block",bl,".rds"))
                   }
  ) %>% bind_rows()
  
  #remove the block files
  lapply(1:nBlock,
         function(bl){
           file.remove(paste0("SNPextract-block",bl,".rds"))
         }
  )
  lapply(1:nBlock,
         function(bl){
           file.remove(paste0("MetaDFextract-block",bl,".rds"))
         }
  )
  print(file.remove(paste0(dirname(msaFilename),"/",segmentFileList)))
  #file.remove(paste0(dirname(msaFilename),"/",segmentFileList))
  tempMetaDump$date<-as.Date(tempMetaDump$date)
  tempMetaDump$weeksince<-as.numeric(floor((tempMetaDump$date-as.Date("2020-01-05"))/7))
  
  weekCounts<-
    tempMetaDump%>%#filter(week=="2021-38")%>%
    mutate(weeksince=floor((as.numeric(date)-18266)/7))%>% #days difference to 5 Jan 2020 (A Sunday)
    select(seqID,weeksince)%>%
    group_by(weeksince)%>%
    mutate(nSeq=n())
  
  tempDump<-tempDump%>%left_join(weekCounts,by="seqID")
  
  saveRDS(tempMetaDump,paste0(outputPrefix,"metaData-blocked.rds"))
  saveRDS(tempDump,paste0(outputPrefix,"SNPs-blocked.rds"))
  
  return(list(
    SNPs=tempDump,
    meta=tempMetaDump
  ))
}

getSeqSNPlocs2<-function(seq,refChars,minPos=374,maxPos=50386){
  #seqChars<-unlist(strsplit(as.character(seq), split = ""))
  seqChars<-str_extract_all(as.character(seq),".")[[1]]
  locs<-(1:length(refChars))[refChars!=seqChars]
#  locs<-(which(refChars!=seqChars))
  locs<-locs[locs>=minPos&locs<=maxPos]
  df<-data.frame(
    pos=locs,
    base=unlist(seqChars[locs]),
    origbase=unlist(refChars[locs])
  )
  #Note: insertions and deletions are included:
  #df<-df[df$base!="-",]
  #df<-df[df$origbase!="-",]
  
  #df<-df[df$base!="N",]
  #df<-df[df$origbase!="N",]
  df<-df[!(df$base %in% c("N","R","Y","K","M","S","W","B","D","H","V","-")),]  # exclude deletions
  df<-df[!(df$origbase %in% c("N","R","Y","K","M","S","W","B","D","H","V")),]
  if(nrow(df)>0)  rownames(df)<-(1:nrow(df))
  return(df)
}
