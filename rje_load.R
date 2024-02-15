########################################################
### RJE_LOAD: Data loading R functions         ~~~~~ ###
### VERSION: 0.2.0                             ~~~~~ ###
### LAST EDIT: 08/04/22                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for mapping/updating sequence identifiers

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version for seqrenamer.R
# v0.2.0 : Updated for different transformation table formats.
# v0.3.0 : Set shortname=TRUE for vector input.
version = "v0.3.0"

####################################### ::: SUMMARY ::: ############################################
#i# loadTable(filename,delimit="ext") - Returns TD$headers vector of headers and TD$data tibble
#i# fastaToList(filename,seqlist=list(),shortname=TRUE) - Returns list of seqname=sequence
#i# seqVector(filename,intvec=TRUE,shortname=TRUE) - Returns list of seqname=vector
#i# filterTableSeq(D,seqnames,seqfield="SeqName",logid="#FILTER") - Returns filtered tibble/dataframe

####################################### ::: TODO ::: ############################################

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


### ~ Filter table to subset of sequences ~ ###
filterTableSeq <- function(D,seqnames,seqfield="SeqName",logid="#FILTER"){
  pren <- nrow(D)
  D <- D[D[[seqfield]] %in% seqnames,]
  if(nchar(logid) > 0){
    logWrite(paste0('#',logid," ",pren," -> ",nrow(buscodb)," entries following filtering to ",length(seqnames)," sequences."))
  }
  return(D)
}


