########################################################
### RJE_LOAD: Data loading R functions         ~~~~~ ###
### VERSION: 0.6.0                             ~~~~~ ###
### LAST EDIT: 20/08/24                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for loading different common data table types and manipulating data.

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version for seqrenamer.R
# v0.2.0 : Updated for different transformation table formats.
# v0.3.0 : Set shortname=TRUE for vector input.
# v0.4.0 : Added saveTable() function.
# v0.5.0 : Added functions for cross-referencing genomic features with regions.
# v0.6.0 : Added listToFasta() function for saving the sequence.
version = "v0.6.0"

####################################### ::: SUMMARY ::: ############################################
#i# loadTable(filename,delimit="ext") - Returns TD$headers vector of headers and TD$data tibble
#i# fastaToList(filename,seqlist=list(),shortname=TRUE) - Returns list of seqname=sequence
#i# seqVector(filename,intvec=TRUE,shortname=TRUE) - Returns list of seqname=vector
#i# filterTableSeq(D,seqnames,seqfield="SeqName",logid="#FILTER") - Returns filtered tibble/dataframe
#i# saveTable(D,suffix,desc) - Generic saving of data frame to TSV file.
#i# listToFasta(filename,seqlist,seqdesc=list(),append=FALSE) - Saves sequence from list to fasta file

####################################### ::: TODO ::: ############################################


################################# ::: LOAD FUNCTIONS ::: ############################################

### R code for loading data into tables.

#i# BUSCOv3 headings:
v3head = c("BuscoID","Status","Contig","Start","End","Score","Length")
#i# BUSCOv5 headings:
v5head = c("BuscoID","Status","Contig","Start","End","Strand","Score","Length","OrthoDBURL","Description")

#i# Universal data loading table for biological data types
loadTable <- function(filename,delimit="ext",comment.char = "#"){
  #<<# Returns TD$headers vector of headers and TD$data tibble
  TD = list()

  # Extract delimiter from file extension
  if(delimit == "ext"){
    ext <- strsplit(filename,".",fixed=TRUE)[[1]]
    ext <- ext[length(ext)]
    delimit <- "\t"
    if(tolower(ext) %in% c("csv")){
      delimit <- ","
    }
    if(tolower(ext) %in% c("txt")){
      delimit <- " "
    }
  }

  # Read file type  
  ftype <- "delim"
  indata <- readLines(filename)
  TD$headers <- indata[startsWith(indata,"#")]
  if(startsWith(indata[1],"# BUSCO")){
    ftype <- "busco"
  }
  if(sum(endsWith(filename,c("gff","gff3")))){
    ftype <- "gff"
  }
  #!# Add SAM file
  TD$type <- ftype
  
  # Read in data
  if(ftype == "gff"){
    gffdb = read.table(filename,sep="\t",fill=TRUE,header=FALSE,row.names = NULL,stringsAsFactors=FALSE)
    logWrite(paste('#GFF',nrow(gffdb),"GFF features loaded from",filename))
    colnames(gffdb) = c('SeqName', 'Source', 'FType', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes')
    gffdb$SeqName = as.character(gffdb$SeqName)
    gffdb$Start <- as.integer(gffdb$Start)
    gffdb$End <- as.integer(gffdb$End)
    TD$data <- gffdb
  }
  if(ftype == "busco"){
    #i# NOTE: URL and Description not always there.
    buscodb = read.table(filename,fill=TRUE,row.names = NULL,sep="\t",quote="",stringsAsFactors=FALSE)
    if(ncol(buscodb) > length(v3head)){
      logWrite(paste("#BUSCOV BUSCO v5 format"))
      colnames(buscodb) = v5head[1:ncol(buscodb)]
    }else{
      logWrite(paste("#BUSCOV BUSCO v3 format"))
      colnames(buscodb) = v3head
    }
    buscodb$Contig = as.character(buscodb$Contig)
    buscodb$Start <- as.integer(buscodb$Start)
    buscodb$End <- as.integer(buscodb$End)
    TD$data <- buscodb
  }
  if(ftype == "delim"){
    TD$data <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",stringsAsFactors=FALSE,comment.char=comment.char)
  }
  logWrite(paste('#LOAD',nrow(TD$data),TD$type,"data rows loaded from",filename))
  return(TD)
}



### ~ FASTA to sequence list (dictionary) ~~~~~~~~~~~~~~~~~~~~~~ ###
# > filename: fasta file
# > seqlist: list of name=sequence loaded from filename
# > shortname: whether to only use first word of name for seqlist name
fastaToList <- function(filename,seqlist=list(),shortname=TRUE){
  indata = readLines(filename)
  i = 1
  logWrite(paste0("Loading sequence data from ",filename,"..."))
  while(i < length(indata)){
    seqname = indata[i]
    if(startsWith(seqname,">")){
      cat(".", file = stderr())
      seqname = substr(seqname,2,nchar(seqname))
      if(shortname){
        seqname <- strsplit(seqname,split=" ",fixed=TRUE)[[1]][1]
      }
      seqlist[[seqname]] = indata[i+1]
      i = i + 2
    }else{
      seqlist[[names(seqlist)[length(seqlist)]]] <- paste0(seqlist[[names(seqlist)[length(seqlist)]]],seqname)
      i = i + 1
    }
  }
  cat("\n", file = stderr())
  logWrite(paste("#FASTA Sequence data for",length(seqlist),"sequences loaded from",filename))
  return(seqlist)
}


### ~ Sequence value list (dictionary) ~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load lengths of depths or kmers in SeqName=Vector list
#i# seqvec = seqVector(filename)
# > filename: numeric vector file
# > intvec: whether to convert to integer
seqVector <- function(filename,intvec=TRUE,shortname=TRUE){
  seqvec = list()
  indata = readLines(filename)
  i = 1
  logWrite(paste0("Loading numeric data from ",filename,"..."))
  cat("", file = stderr())
  while(i < length(indata)){
    cat(".", file = stderr())
    seqname = indata[i]
    seqname = substr(seqname,2,nchar(seqname))
    if(shortname){
      seqname <- strsplit(seqname,split=" ",fixed=TRUE)[[1]][1]
    }
    if(intvec){
      seqvec[[seqname]] = as.integer(strsplit(indata[i+1]," ")[[1]])
    }else{
      seqvec[[seqname]] = as.numeric(strsplit(indata[i+1]," ")[[1]])
    }
    i = i + 2
  }
  cat("\n", file = stderr())
  logWrite(paste("#SEQVEC Numeric data for",length(seqvec),"sequences loaded from",filename))
  return(seqvec)  
}

################################# ::: TABLE FILTERING ::: ############################################

### ~ Filter table to subset of sequences ~ ###
filterTableSeq <- function(D,seqnames,seqfield="SeqName",logid="#FILTER"){
  pren <- nrow(D)
  D <- D[D[[seqfield]] %in% seqnames,]
  if(nchar(logid) > 0){
    logWrite(paste0('#',logid," ",pren," -> ",nrow(buscodb)," entries following filtering to ",length(seqnames)," sequences."))
  }
  return(D)
}


################################# ::: SAVE FUNCTIONS ::: ############################################

### ~ Generic table saving function ~ ###
saveTable <- function(df,suffix,desc){
  if(nrow(df)>0){
    outfile = paste(settings$basefile,suffix,"tsv",sep=".",collapse=".")
    logWrite(paste("#SAVE",nrow(df),desc,"output to",outfile))
    write.table(df,outfile,sep="\t",quote=FALSE,row.names=FALSE)
  }
}

### ~ Sequence list (dictionary) to FASTA ~~~~~~~~~~~~~~~~~~~~~~ ###
# > filename: fasta file
# > seqlist: list of name=sequence loaded from filename
# > seqdesc: optional list of names to descriptions
listToFasta <- function(filename,seqlist,seqdesc=list(),append=FALSE){
  # Write output to filename
  for(seqname in sort(names(seqlist))){
    outname <- seqname
    if(seqname %in% names(seqdesc)){
      outname <- paste(seqname, seqdesc[[seqname]])
    }
    cat(paste0(">",outname,"\n"), file=filename,append=append)
    append=TRUE
    cat(paste0(seqlist[[seqname]],"\n"), file=filename, append=TRUE)
  }
}


################################# ::: FEATURE FUNCTIONS ::: ############################################

### ~~~~~~~~~~ Generic function for returning features spanning regions ~~~~~~~~~~~~~~~~~~~ ###
#i# > regD = Region tibble with SeqName, Start, End[, Strand]
#i# > ftD = Feature tibble with SeqName, Start, End[, Strand] -> Becomes ftStart, ftEnd, ftStrand
#i# > stranded=FALSE (TRUE/FALSE) : whether to restrict overlaps to matching strands. Strand "." will match either.
#i# > overlap=any (any/span/whole/exact) : whether to include any form of overlap, only ftD entries that span regD, ftD entries within a regD region, or exact matches only
#i# > flank=0 [int] : additional flank to add to each end of region to meet criteria. Negative flank will shrink region and thus require a buffer
#i# > maxsize=100 [float] : Maximum predicted intermediate object size (MB) to try quick/simple regFeatures.
#i# < spanD = Joined tibble where feature overlap conditions are met. 
#i# NOTE: Will contain all fields from regD and ftD. These should be unique else will be used for join in addition to SeqName.
#i# NOTE: The simple version uses a many-to-many join and will call regFeatures() if the predicted combo exceeds the maxsize setting (MB)
regFeaturesSimple <- function(regD,ftD,stranded=FALSE,overlap="any",flank=0){
  #i# Join by SeqName
  jD <- full_join(regD, ftD, relationship = "many-to-many")
  #i# Stranded
  if(stranded){
    jD <- jD %>% filter(! Strand=="+" & ftStrand == "-", ! Strand=="-" & ftStrand == "+")
  }
  #i# Flank adjust
  if(flank != 0){
    jD <- jD %>% mutate(Start=Start-!!flank, End=End+!!flank)
  }
  #i# Filter based on matches
  filtered <- FALSE
  if(overlap == "any"){
    filtered <- TRUE
    jD <- jD %>% filter(ftEnd>=Start,ftStart<=End)
  }
  if(overlap == "exact"){
    filtered <- TRUE
    jD <- jD %>% filter(ftEnd==Start,ftStart==End)
  }
  if(overlap == "span"){
    filtered <- TRUE
    jD <- jD %>% filter(ftEnd>=End,ftStart<=Start)
  }
  if(overlap == "whole"){
    filtered <- TRUE
    jD <- jD %>% filter(ftEnd<=End,ftStart>=Start)
  }
  if(! filtered){
    logWrite(paste("Did not recognise regFeatures overlap =",overlap))
    jD <- jD %>% mutate(Tmp=Start) %>% filter(Start > Tmp) %>% select(-Tmp)
  }  
  #i# Reverse flank adjust
  if(flank != 0){
    jD <- jD %>% mutate(Start=Start+!!flank, End=End-!!flank)
  }
  logWrite(paste(nrow(jD),"features with",overlap,"matches to regions."))
  return(jD)
}

#i# Full version will incorporate a slower option that does not include a full_join
regFeatures <- function(regD,ftD,stranded=FALSE,overlap="any",flank=0,maxsize=100){
  #i# Setup tibbles
  if(! "Strand" %in% rownames(regD)){ regD$Strand <- "." }
  if(! "Strand" %in% rownames(ftD)){ ftD$Strand <- "." }
  ftD <- ftD %>% rename(ftStart=Start, ftEnd=End, ftStrand=Strand)
  #i# Check size and try simple version
  jrows <- nrow(full_join(regD %>% select(SeqName), ftD %>% select(SeqName), relationship="many-to-many"))
  jsize <- ((object.size(regD) / nrow(regD)) + (object.size(ftD) / nrow(ftD))) * jrows / 1e6
  if(jsize <= maxsize){ return(regFeaturesSimple(regD,ftD,stranded,overlap,flank)) }
  else{
    logWrite(paste("Predicted", round(jsize,2), "MB Region-Feature tibble join"))
  }
  #!# Add more complex code. Do it by SeqName?
  seqnames <- unique(regD$SeqName)
  seqnames <- seqnames[seqnames %in% unique(ftD$SeqName)]
  jD <- tibble()
  for(i in 1:length(seqnames)){
    sD <- regFeaturesSimple(regD %>% filter(SeqName == !!seqnames[i]),ftD %>% filter(SeqName == !!seqnames[i]),stranded,overlap,flank)
    if(i < 2){
      jD <- sD
    }else{
      jD <- bind_rows(jD,sD)
    }
  }
  return(jD)
}