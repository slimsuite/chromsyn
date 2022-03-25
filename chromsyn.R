########################################################
### ChromSyn Synteny Plot functions            ~~~~~ ###
### VERSION: 0.7.0                             ~~~~~ ###
### LAST EDIT: 23/03/22                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: richard.edwards@unsw.edu.au       ~~~~~ ###
########################################################

# This script is for homology-based chromosome synteny plots

####################################### ::: HISTORY ::: ############################################
# v0.0.0 : Initial version based on SynBad data exploration
# v0.1.0 : Added dev mode, expanded focus and added splitting of plots if too big.
# v0.2.0 : Added orphans=T/F to control whether scaffolds without BUSCO genes should be plotted. Added regdata=TSV input.
# v0.3.0 : Added excel file output. Added names to chromsyn plots.
# v0.4.0 : Added some additional controls for the chromosome labels, plus minlen=INT filter.
# v0.5.0 : Added optional reading of Tel5 and Tel3 files from sequences table. Swapped meaning of seqsort and seqorder. 
# v0.6.0 : Added optional TIDK parsing for additional telomere prediction.
# v0.7.0 : Add PNG output and optional gap and feature table parsing.
version = "v0.7.0"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript chromsyn.R [sequences=FOFN] [busco=FOFN] [focus=X] [seqorder=LIST] [order=LIST] [basefile=FILE] [minregion=INTbp] [align=X]
# : sequences=FOFN = File of PREFIX FILE with sequence names and lengths (name & length, or SeqName & SeqLen fields) [sequences.fofn]
# : busco=FOFN = File of PREFIX FILE with full BUSCO table results. Used to identify orthologous regions. [busco.fofn]
# : tidk=FOFN = Optional file of PREFIX FILE with TIDK search results. [tidk.fofn]
# : gaps=FOFN = Optional file of PREFIX FILE with TIDK search results. [gaps.fofn]
# : ft=FOFN = Optional file of PREFIX FILE with TIDK search results. [ft.fofn]
# : regdata=TSV = File of Genome, HitGenome, Seqname, Start, End, Strand, Hit, HitStart, HitEnd
# : focus=X = If given will orient all chromosomes to this assembly
# : orient=X = Mode for sequence orientation (none/focus/auto)
# : seqsort=none/focus/auto/FILE = Optional ordering strategy for other assemblies [auto]
# : seqorder=LIST = Optional ordering of the chromsomes for the focal assembly
# : order=LIST = File containing the Prefixes to include in vertical order. If missing will use sequences=FOFN.
# : basefile=FILE = Prefix for outputs [chromsyn]
# : plotdir=PATH = output path for graphics
# : minlen=INT = minimum length for a chromosome/scaffold to be included in synteny blocks/plots [0]
# : minregion=INT = minimum length for mapped regions to be included in plots [50000]
# : minbusco=INT = minimum number of BUSCO genes to be included in Syntenic block [1]
# : maxskip=0 = maximum number of BUSCO genes to skip and still be a syntenic block [0]
# : orphans=T/F = whether to include scaffolds that have no BUSCO genes [True]
# : tidkcutoff=INT = TIDK count cutoff for identifying a telomere [50]
# : align=X = alignment strategy for plotting chromosomes (left/right/centre/justify) [justify]
# : ygap=INT = vertical gap between chromosomes [4]
# : scale=X = units in basepairs for setting the x-axis scale [Mb]
# : textshift=NUM = offset for printing chromosome names [0.3]
# : ticks=INT = distance between tickmarks [5e7]
# : pdfwidth=NUM = PDF width [20]
# : pdfheight=NUM = over-ride for standard calculated PDF height [0]
# : pdfscale=NUM = over-ride for PDF output scale [1]
# : namesize=NUM = scaling factor for the Genome names in PDF plots [1]
# : labelsize=NUM = scaling factor for the chromosome names in PDF plots [1]
# : labels=T/F = whether to print chromosome name labels [TRUE]
# : opacity=NUM = Opacity of synteny regions (0-1) [0.3]
# : debug=T/F = whether to switch on additional debugging outputs [FALSE]
# : dev=T/F = whether to switch on dev mode during code updates [FALSE]

#i# Python code:
# slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
# rdir = '%slibraries/r/' % slimsuitepath
# os.popen('Rscript {0}chromsyn.R {1}'.format(rdir, optionstr)).readlines()

#i# Usage within R:
# Set an override vector of commandline arguments: this will replace argvec read from the commandline
# Use source() to run the script:
# source("$PATH/chromsyn.R")

####################################### ::: OUTPUTS ::: ############################################
#!# List the outputs here

####################################### ::: TODO ::: ############################################
#!# Test min number of BUSCO genes to compress into region? minbusco=INT
#!# seqorder=FILE with Genome and comma-separated sequences
#!# Save and load tables for pure plotting
#!# Add to BUSCOMP
#!# Test ability to skip missing BUSCO genes if needed
#!# Test settings for pdfwidth and pdfheight over-ride.
#!# Need to enable no loading of BUSCO data if regdata=TSV is given.
#!# Add output of synteny block table.
#!# Get bitmap output working. (Resolution is awful for unknown reasons so disabled)
# : pngscale=INT = PDF to PNG scaling [100]
# : pointsize=NUM = PNG text point size [24]
#!# Add loading of optional ftdb FOFN of features to plot with SeqName, Pos and optional Strand ("."), Color and Shape
#!# Add file checks and error messages!

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(busco="busco.fofn",sequences="sequences.fofn",order="",regdata="",
                minregion=50000,align="justify",ygap=4,seqorder="",orient="auto",seqsort="auto",
                #pngwidth=1200,pngheight=900,
                #pointsize=24,pngscale=100,
                plotdir="./",
                tidkcutoff=50,tidk="tidk.fofn",gaps="gaps.fofn",ft="ft.fofn",
                minbusco=1,maxskip=0,orphans=TRUE,minlen=0,opacity=0.3,
                pdfwidth=20,pdfheight=0,pdfscale=1,namesize=1,labelsize=1,labels=TRUE,
                scale = "Mb",textshift = 0.3,ticks=5e7,rscript=TRUE,
                basefile="chromsyn",focus="",debug=FALSE,dev=FALSE,
                outlog=stdout())

settings <- defaults
argvec = commandArgs(TRUE)
if("override" %in% ls()){
  argvec = override
  settings$rscript = FALSE
}
for(cmd in argvec){
  cmdv = strsplit(cmd,'=',TRUE)[[1]]
  if(length(cmdv) > 1){
    settings[[cmdv[1]]] = cmdv[2]
  }else{
    settings[[cmdv[1]]] = TRUE    
  }
}
for(cmd in c("pngwidth","pngheight","pointsize","minregion","ygap","minbusco","maxskip","minlen","pngscale","tidkcutoff")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
for(cmd in c("textshift","ticks","pdfwidth","pdfheight","pdfscale","namesize","labelsize","opacity")){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
if(settings$opacity < 0.1){
  settings$opacity <- 0.1
}
if(settings$opacity > 1){
  settings$opacity <- 1
}
for(cmd in c("order","seqorder")){
  if(length(settings[[cmd]]) > 0){
    settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
  }else{
    settings[[cmd]] <- defaults[[cmd]]
  }
}
for(cmd in c("debug","dev","orphans","labels")){
  settings[[cmd]] = as.logical(settings[[cmd]])
}

oldwarn <- getOption("warn")

if(settings$debug){
  writeLines(argvec)
}else{
  options(warn = -1)
}

### ~ logWrite function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite <- function(logstr){
  writeLines(paste0("[",date(),"] ",logstr),con=settings$outlog)
}
logWrite(paste("#RCODE ChromSyn.R:",version))
for(cmd in names(settings)[order(names(settings))]){
  logWrite(paste("CMD:",cmd,"=",paste0(settings[[cmd]],collapse=",")))
}
dir.create(settings$plotdir, showWarnings = FALSE)

### ~ Load packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! "tidyverse" %in% installed.packages()[,"Package"]){
  install.packages("tidyverse")
}
library(tidyverse)
library(RColorBrewer)
library(gtools)
settings$ggstatsplot = "ggstatsplot" %in% installed.packages()[,"Package"]
if(settings$ggstatsplot){
  library(ggstatsplot)
}
settings$writexl = "writexl" %in% installed.packages()[,"Package"]
if(settings$writexl){
  library(writexl)
}else{
  logWrite("#XLXS Install writexl package for compiled Excel file output.")
}

####################################### ::: FUNCTIONS ::: ############################################

##### ======================== Loading data functions ======================== #####

### ~ Load FOFN File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into FOFN table (genome, file)
#i# fofndb = fileTable(filename,delimit="\t")
fileTable <- function(filename,delimit=" ",genomes=c()){
  regdb <- read.table(filename,fill=TRUE,sep=delimit,header=FALSE,row.names = NULL,quote="\"",comment.char="") %>%
    select(1:2)
  colnames(regdb) <- c("Genome","Filename")
  logWrite(paste('#FOFN',nrow(regdb),"filenames loaded from",filename))
  if(length(genomes) > 0){
    regdb <- regdb[regdb$Genome %in% genomes,]
    logWrite(paste('#FOFN',nrow(regdb),"filenames after filtering to recognised genomes."))
  }
  return(regdb)
}

### ~ Load Sequences File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into Sequences table (name, length)
#i# seqdb = seqTable(filename,delimit="\t")
seqTable <- function(genome,filename,delimit="\t"){
  regdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="")
  regdb$Genome <- genome
  if("name" %in% colnames(regdb)){
    regdb <- regdb %>% rename(SeqName=name,SeqLen=length)
  }
  if(! "Tel5" %in% colnames(regdb)){
    regdb$Tel5 <- FALSE
    regdb$Tel3 <- FALSE
  }
  regdb$Tel5 <- as.logical(regdb$Tel5)
  regdb$Tel5[is.na(regdb$Tel5)] <- FALSE
  regdb$Tel3 <- as.logical(regdb$Tel3)
  regdb$Tel3[is.na(regdb$Tel3)] <- FALSE
  regdb <- regdb %>% select(Genome,SeqName,SeqLen,Tel5,Tel3)
  logWrite(paste('#SEQS',nrow(regdb),genome,"sequences loaded from",filename))
  regdb <- regdb[regdb$SeqLen >= settings$minlen,]
  logWrite(paste('#SEQS',nrow(regdb),genome,"sequences meet minlen cutoff of",settings$minlen,"bp"))
  return(regdb)
}

### ~ Load BUSCO File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load BUSCO table into tibble
#i# buscodb = buscoTable(filename)
#i# BUSCOv3 headings:
v3head = c("BuscoID","Status","Contig","Start","End","Score","Length")
#i# BUSCOv5 headings:
v5head = c("BuscoID","Status","Contig","Start","End","Strand","Score","Length","OrthoDBURL","Description")
#i# NOTE: URL and Description not always there.
buscoTable <- function(genome,filename,seqnames=c()){
  buscodb = read.table(filename,fill=TRUE,row.names = NULL,sep="\t",quote="")
  if(ncol(buscodb) > length(v3head)){
    logWrite(paste("#BUSCOV BUSCO v5 format"))
    colnames(buscodb) = v5head[1:ncol(buscodb)]
  }else{
    logWrite(paste("#BUSCOV BUSCO v3 format"))
    colnames(buscodb) = v3head
  }
  buscodb$Contig = as.character(buscodb$Contig)
  buscodb <- buscodb[buscodb$Status == "Complete",]
  logWrite(paste('#BUSCO',nrow(buscodb),genome,"Complete BUSCO genes loaded from",filename))
  buscodb$Genome <- genome
  buscodb$Start <- as.integer(buscodb$Start)
  buscodb$End <- as.integer(buscodb$End)
  buscodb <- buscodb %>% rename(SeqName=Contig) %>% select(Genome,SeqName,Start,End,Strand,BuscoID)
  #i# Filter to sequences if given
  if(length(seqnames) > 0){
    buscodb <- buscodb[buscodb$SeqName %in% seqnames,]
    logWrite(paste('#BUSCO',nrow(buscodb),genome,"Complete BUSCO genes following filtering to",length(seqnames),"sequences."))
  }
  return(buscodb)
}

### ~ Load Regions File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into Regions table (Genome, HitGenome, Seqname, Start, End, Strand, Hit, HitStart, HitEnd)
#i# regdb = regTable(filename,delimit="\t")
regTable <- function(filename,delimit="\t"){
  regdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="")
  if(! "Length" %in% colnames(regdb)){
    regdb <- mutate(regdb,Length=End-Start+1)
  }
  if(! "HitLength" %in% colnames(regdb)){
    regdb <- mutate(regdb,HitLength=HitEnd-HitStart+1)
  }
  regdb <- select(regdb,Genome, HitGenome, Seqname, Start, End, Strand, Hit, HitStart, HitEnd, Length, HitLength)
  #?# Check the fields
  logWrite(paste('#REGION',nrow(regdb),"linked regions loaded from",filename))
  return(regdb)
}

### ~ Load TIDK File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into TIDK table (Genome, Seqname, SeqWin, Tel5, Tel3, TelSeq)
#i# Headers: id,window,forward_repeat_number,reverse_repeat_number,telomeric_repeat
#i# tidkdb = tidkTable(genome,filename,delimit=",")
tidkTable <- function(genome,filename,delimit=","){
  tidkdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="")
  colnames(tidkdb) <- c("SeqName", "SeqWin", "Tel5", "Tel3", "TelSeq")
  tidkdb$Genome <- genome
  tidkdb <- tidkdb[tidkdb$Tel5 >= settings$tidkcutoff | tidkdb$Tel3 >= settings$tidkcutoff,] 
  tidk <- select(tidkdb, Genome, SeqName, SeqWin, Tel5, Tel3)
  #i# Strand and Pos conversion
  tel5 <- tidk[tidk$Tel5 >= settings$tidkcutoff,] %>% mutate(Pos=SeqWin, Strand="+") %>% select(Genome, SeqName, Pos, Strand)
  # if(nrow(tel5)){
  #   tel5$Strand <- "+"
  # }
  tel3 <- tidk[tidk$Tel3 >= settings$tidkcutoff,] %>% mutate(Pos=SeqWin, Strand="-") %>% select(Genome, SeqName, Pos, Strand)
  # if(nrow(tel3)){
  #   tel3$Strand <- "-"
  # }
  tidkdb <- bind_rows(tel5,tel3) 
  #?# Check the fields: Genome, SeqName, Pos, Strand
  logWrite(paste('#TIDK',nrow(tidkdb),"TIDK telomere windows loaded from",filename))
  return(tidkdb)
}

### ~ Load Feature/Gap File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into feature or gaps (Genome, Seqname, Pos, Strand, Col, Shape)
#i# ftdb = ftTable(genome,filename,delimit=",")
ftTable <- function(genome,filename,colour="darkgreen",shape=15){
  delimit <- "\t"
  if(endsWith(filename,".csv")){
    delimit <- ","
  }
  ftdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="")
  ftdb$Genome <- genome
  for(ftfield in c("SeqName","Pos","Strand","Col","Shape","Start","End")){
    if(tolower(ftfield) %in% colnames(ftdb)){
      ftdb[[ftfield]] <- ftdb[[tolower(ftfield)]]
    }
  }
  if(! "Pos" %in% colnames(ftdb)){
    ftdb <- mutate(ftdb,Pos=(End+Start)/2)
  }
  if(! "Strand" %in% colnames(ftdb)){
    ftdb$Strand <- "."  # Unstranded
  }
  if(! "Col" %in% colnames(ftdb)){
    ftdb$Col <- colour  # Defaults to steelblue
  }
  if(! "Shape" %in% colnames(ftdb)){
    ftdb$Shape <- shape  # Defaults to square (3 = + for gaps)
  }
  ftdb <- select(ftdb,Genome, SeqName, Pos, Strand, Col, Shape)
  logWrite(paste('#FT',nrow(ftdb),"features loaded from",filename))
  return(ftdb)
}

##### ======================== Sequence stats functions ======================== #####

### ~ Sequence length ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Call seqLen(seqdb,genome,seqname)
#i# Returns the length of sequence "seqname" in genome "genome"
seqLen <- function(seqdb,genome,seqname){
  seqlen <- seqdb %>% filter(Genome==genome,SeqName==seqname) %>% select(SeqLen)
  return(seqlen$SeqLen[1])
}

### ~ Sequence reversed? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Call seqRev(seqdb,genome,seqname)
#i# Returns whether sequence "seqname" in genome "genome" has been reversed
seqRev <- function(seqdb,genome,seqname){
  seqrev <- seqdb %>% filter(Genome==genome,SeqName==seqname) %>% select(Rev)
  return(seqrev$Rev[1])
}

################################## ::: PLOTTING FUNCTIONS ::: #######################################


################################## ::: RUN CODE ::: #######################################

##### ======================== Report key inputs ======================== #####

### ~ Key inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Sequences files = vital input
seqfofn = settings$sequences
logWrite(paste("Sequence FOFN File:",seqfofn))
if(! file.exists(seqfofn)){
  if(settings$rscript){
    quit("no",2)  
  }else{
    stop(paste("Cannot find sequence FOFN file:",seqfofn))
  }
}
#i# Region table (does not need BUSCO)
if(! file.exists(settings$regdata)){
  if(settings$regdata != ""){
    logWrite(paste("#REGDAT Cannot find region data file:",settings$regdata))
    settings$regdata <- ""
  }else{
    logWrite(paste("No region data file given."))
  }
}else{
  logWrite(paste("Region data file:",settings$regdata))
}
#i# BUSCO tables
buscofofn = settings$busco
logWrite(paste("BUSCO FOFN file:",buscofofn))
if(! file.exists(buscofofn) & settings$regdata == ""){
  if(settings$rscript){
    quit("no",2)  
  }else{
    stop(paste("Cannot find BUSCO FOFN file:",buscofofn))
  }
}
#i# Optional extra tables: TIDK, gaps, ft
filechecks <- list(tidk="TIDK", gaps="Assembly gap", ft="features table")
for(ftype in names(filechecks)){
  ffile <- settings[[ftype]]
  if(file.exists(ffile)){
    logWrite(paste(filechecks[[ftype]],"FOFN file:",ffile))
  }else{
    if(ffile != ""){
      logWrite(paste(filechecks[[ftype]],"FOFN file not found:",ffile))
    }
  }
}
logWrite('#RCODE Setup complete.')

##### ======================== Load data ======================== #####

### ~ Load sequence and BUSCO data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load file names
seqfiles <- fileTable(settings$sequences) %>% rename(SeqFile=Filename)
buscofiles <- fileTable(settings$busco,genomes=seqfiles$Genome) %>% rename(BUSCO=Filename)
if(length(settings$order)<1 | (length(settings$order)==1 & settings$order[1] == "")){
  settings$order <- seqfiles$Genome
}
logWrite(paste("Genomes (order=LIST):",paste(settings$order,collapse=", ")))
#i# Combine into table
gendb <- inner_join(seqfiles,buscofiles)
# - TIDK
if(file.exists(settings$tidk)){
  tidkfiles <- fileTable(settings$tidk,genomes=seqfiles$Genome) %>% rename(TIDK=Filename)
  gendb <- inner_join(gendb,tidkfiles)
}else{
  gendb$TIDK <- NA
}
# - gaps
if(file.exists(settings$gaps)){
  gapfiles <- fileTable(settings$gaps,genomes=seqfiles$Genome) %>% rename(gaps=Filename)
  gendb <- inner_join(gendb,gapfiles)
}else{
  gendb$gaps <- NA
}
# - features
if(file.exists(settings$ft)){
  ftfiles <- fileTable(settings$ft,genomes=seqfiles$Genome) %>% rename(features=Filename)
  gendb <- inner_join(gendb,ftfiles)
}else{
  gendb$features <- NA
}

#i# Update sequence order
settings$order[settings$order %in% gendb$Genome]
gendb <- gendb %>% filter(Genome %in% settings$order)
rownames(gendb) <- gendb$Genome
logWrite(paste("#GENOME ",nrow(gendb),"genomes:",paste(settings$order,collapse=", ")))
genomes <- settings$order
#i# Load actual sequence and busco data
seqdb <- data.frame()
busdb <- data.frame()
teldb <- data.frame()
gapdb <- data.frame()
ftdb <- data.frame()
for(genome in genomes){
  #i# Sequences
  filename <- gendb[genome,"SeqFile"]
  newseqdb <- seqTable(genome,filename)
  if(nrow(seqdb) > 0){
    seqdb <- bind_rows(seqdb,newseqdb)
  }else{
    seqdb <- newseqdb
  }
  #i# BUSCO
  filename <- gendb[genome,"BUSCO"]
  if(nrow(busdb) > 0){
    busdb <- bind_rows(busdb,buscoTable(genome,filename,newseqdb$SeqName))
  }else{
    busdb <- buscoTable(genome,filename,newseqdb$SeqName)
  }
  #i# TIDK
  filename <- gendb[genome,"TIDK"]
  if(! is.na(filename)){
    if(nrow(teldb) > 0){
      teldb <- bind_rows(teldb,tidkTable(genome,filename))
    }else{
      teldb <- tidkTable(genome,filename)
    }
  }
  #i# Gaps
  filename <- gendb[genome,"gaps"]
  if(! is.na(filename)){
    newdb <- ftTable(genome,filename,colour="red",shape=3)
    if(nrow(gapdb) > 0){
      gapdb <- bind_rows(gapdb,newdb)
    }else{
      gapdb <- newdb
    }
  }
  #i# Features
  filename <- gendb[genome,"features"]
  if(! is.na(filename)){
    newdb <- ftTable(genome,filename)
    if(nrow(ftdb) > 0){
      ftdb <- bind_rows(ftdb,newdb)
    }else{
      ftdb <- newdb
    }
  }
}
logWrite(paste("#SEQS ",nrow(seqdb),"sequences loaded in total."))
logWrite(paste("#BUSCO ",nrow(busdb),"BUSCO genes loaded in total."))
logWrite(paste("#TIDK ",nrow(teldb),"TIDK telomere windows loaded in total."))
busdb <- filter(busdb,SeqName %in% seqdb$SeqName)
logWrite(paste("#BUSCO ",nrow(busdb),"BUSCO genes for loaded sequences."))

### ~ Optional sequence filter based on BUSCO data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! settings$orphans){
  fulln <- nrow(seqdb)
  buscseq <- group_by(busdb,Genome,SeqName) %>% summarise(Genome=first(Genome),SeqName=first(SeqName),N=n())
  seqdb <- inner_join(seqdb,buscseq) %>% filter(N>0) %>% select(-N)
  logWrite(paste("#ORPHAN Removed orphan sequences w/o BUSCOs:",fulln,"->",nrow(seqdb),"sequences."))
}


##### ======================== Compile data ======================== #####
#i# Sort for the counter
busdb <- arrange(busdb,Genome,SeqName,Start,End)
#i# Join on BuscoID to make homologous regions
buscopy <- busdb %>% rename(HitGenome=Genome,Hit=SeqName,HitStart=Start,HitEnd=End,HitStrand=Strand)
regdb <- inner_join(busdb,buscopy)
#i# Filter and tidy the regdb table
regdb <- filter(regdb,Genome != HitGenome)
fwd <- regdb$Strand == regdb$HitStrand
regdb$Strand <- "-"
regdb$Strand[fwd] <- "+"
regdb <- select(regdb,Genome,HitGenome,SeqName,Start,End,Strand,Hit,HitStart,HitEnd,BuscoID) %>%
  arrange(Genome,HitGenome,SeqName,Start,End)

#?# Should the regdata table be added here? Currently: no. Will be better to have converters to make BUSCO-like tables of 
#?# regions to link.

#i# Add counter for each Genome (no need to restart numbering)
regdb$Counter <- 1:nrow(regdb)
regdb <- arrange(regdb,Genome,HitGenome,Hit,HitStart,HitEnd)
regdb$HitCounter <- 1:nrow(regdb)
regdb <- arrange(regdb,Genome,HitGenome,SeqName,Start,End)

#i# - Add synteny block ID
regdb$Block <- 1
#i# - compress by block
progv <- as.integer(1:20 * (nrow(regdb)/20))
logWrite('Generating synteny blocks...')
noblock <- settings$maxskip + 2
for(i in 2:nrow(regdb)){
  j <- i - 1
  #i# Updated this to allow a number of moved genes using settings$maxskip
  #block <- abs(regdb$Counter[i] - regdb$Counter[j]) == 1 & abs(regdb$HitCounter[i] - regdb$HitCounter[j]) == 1
  block <- abs(regdb$Counter[i] - regdb$Counter[j]) < noblock & abs(regdb$HitCounter[i] - regdb$HitCounter[j]) < noblock
  if(block){
    for(field in c("Genome","SeqName","HitGenome","Hit","Strand"))
      if(regdb[i,field] != regdb[j,field]){
        block <- FALSE
      }
  }
  if(block){
    regdb$Block[i] <- regdb$Block[j]
  }else{
    regdb$Block[i] <- regdb$Block[j] + 1
  }
  if(i %in% progv){
    cat(paste0("...",which(progv == i)*5,"%"), file = stderr())
  }
}
cat("\n", file = stderr())
logWrite(paste('#BLOCK Generated',max(regdb$Block),"synteny blocks."))
#i# - compress by block
backdb <- regdb
regdb <- backdb %>% group_by(Block) %>% 
  summarise(Genome=first(Genome),HitGenome=first(HitGenome),SeqName=first(SeqName),Start=min(Start),End=max(End),Strand=first(Strand),Hit=first(Hit),HitStart=min(HitStart),HitEnd=max(HitEnd),BuscoID=n()) %>% 
  select(-Block) %>%
  mutate(Length=End-Start+1,HitLength=HitEnd-HitStart+1) %>%
  filter(Length>=settings$minregion,HitLength>=settings$minregion)
logWrite(paste('#BLOCK Reduced to',nrow(regdb),"synteny blocks based on minregion=INT filtering."))
if(settings$minbusco > 1){
  regdb <- filter(regdb,BuscoID>=settings$minbusco)
  logWrite(paste('#BLOCK Reduced to',nrow(regdb),"synteny blocks based on minbusco=INT filtering."))
}


##### ======================== Add pre-generated regdb data ======================== #####
if(settings$regdata != ""){
  linkdb <- regTable(settings$regdata) %>%
    filter(Length>=settings$minregion,HitLength>=settings$minregion)
  logWrite(paste('#REGION Reduced to',nrow(linkdb),"synteny blocks based on minregion=INT filtering."))
  if(nrow(regdb) > 0){
    regdb <- bind_rows(regdb,linkdb)
    logWrite(paste('#BLOCK',nrow(regdb),"combined region+BUSCO synteny blocks."))
  }else{
    regdb <- linkdb
  }
}


##### ======================== Update Gap Colours ======================== #####
#!# To do:
# - Gaps within synteny blocks go blue
# - Gaps outside synteny blocks stay red


##### ======================== Calculate best matches between genomes ======================== #####

### ~ Identify best matching chromosome for each other genome ~ ###
bestdb <- regdb %>% mutate(iG=HitStart/1e9,jG=HitEnd/1e9) %>%
  mutate(HitPos=(HitLength)*(iG+jG)/2) %>% select(-iG,-jG) %>%
#bestdb <- bestdb %>% 
  group_by(Genome,HitGenome,SeqName,Hit,Strand) %>%
  summarise(Genome=first(Genome),HitGenome=first(HitGenome),SeqName=first(SeqName),Hit=first(Hit),Strand=first(Strand),Length=sum(Length),HitLength=sum(HitLength),BuscoID=sum(BuscoID),HitPos=sum(HitPos)) %>%
  arrange(Genome,HitGenome,SeqName,desc(Length))
bestdb$Fwd <- 0
bestdb$Bwd <- 0
fwd <- bestdb$Strand == "+"
bestdb$Fwd[fwd] <- bestdb$Length[fwd]
bestdb$Bwd[!fwd] <- bestdb$Length[!fwd]
bestdb <- group_by(bestdb,Genome,HitGenome,SeqName,Hit) %>%
  summarise(Genome=first(Genome),HitGenome=first(HitGenome),SeqName=first(SeqName),Hit=first(Hit),Length=sum(Length),BuscoID=sum(BuscoID),Fwd=sum(Fwd),Bwd=sum(Bwd),HitLength=sum(HitLength),HitPos=sum(HitPos)) %>%
  arrange(Genome,HitGenome,SeqName,desc(Length)) %>%
  mutate(HitPosMb=1000*HitPos/HitLength) %>% select(-HitLength,-HitPos)
#i# Make a copy as hitdb to be saved in Excel file
hitdb <- bestdb
#i# Can now assign each SeqName to best Hit per HitGenome
bestdb <- slice_max(bestdb,Length,n=1)

### ~ Establish focus ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# If focus not given, choose the one with maximal summed Fwd coverage
if(! settings$focus %in% seqdb$Genome){
  focdb <- group_by(bestdb,Genome) %>% summarise(Genome=first(Genome),Fwd=sum(Fwd)) %>%
    slice_max(Fwd,n=1)
  settings$focus <- focdb$Genome[1]
}
logWrite(paste('#FOCUS Focal genome for orientation:',settings$focus))

### ~ Clean up objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#rm(backdb)
rm(buscopy)

##### ======================== Orient data ======================== #####
#i# This next section will order orient the sequences from each genome relative to the focus.
#i# This is done based on the orient=X setting:
#i# - auto = Propagate out from focus genome and sort/orient off neighbour [default]
#i# - focus = Fix the order and direction of everything from the focus genome
#i# - none = Use the loaded order and direction

#i# The critical data/variables at this point are:
# bestdb = best Hit for (Genome,SeqName) in each HitGenome (used to orient/place SeqName based on Hit)
# - Genome HitGenome SeqName Hit Length BuscoID Fwd Bwd HitPosMb
# >> Reversing a SeqName will flip Fwd and Bwd; Reversing a Hit will invert HitPosMb during ordering
# busdb = Positions of the BUSCO genes
# - Genome SeqName Start End Strand BuscoID
# >> Not currently edited if reversed
# gendb = Genomes and files
# - Genome SeqFile BUSCO
# regdb = Compressed regions of overlap for plotting
# - Genome HitGenome SeqName Start End Strand Hit HitStart HitEnd BuscoID Length HitLength
# >> Reversing a SeqName will flip Strand and invert Start and End versus seqLen(Genome,SeqName)
# >> Reversing a Hit will flip Strand and invert HitStart and HitEnd versus seqLen(HitGenome,Hit)
# seqdb = Sequence names and lengths
# - Genome SeqName SeqLen
# >> Also has Rev (TRUE/FALSE) and text (Display text, with R if reversed) added

### ~ Reorient sequences if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
genfocus <- list()
gnum <- length(settings$order)
#i# Work out which sequences to flip (seqdb$Rev == TRUE) but do not actually edit any of the other tables
#i# - inversion of reversed sequences happens during plotting.
if(! settings$orient %in% c("auto","none","focus")){
  logWrite(paste0("Could not recognise orient=",settings$orient," (auto/focus/none): set to 'auto'."))
  settings$orient <- "auto"
}
if(settings$orient == "none"){
  seqdb$Rev <- FALSE
  for(i in 1:gnum){
    genfocus[[settings$order[i]]] <- settings$order[i]
  }
}
#i# focus mode will flip and order everything relative to the query
if(settings$orient == "focus"){
  revdb <- ungroup(bestdb) %>% filter(HitGenome==settings$focus,Bwd>Fwd) %>% select(Genome,SeqName) %>% 
    mutate(Rev=TRUE)
  logWrite(paste('#REVSEQ',nrow(revdb),'sequences need reverse orientation.'))
  seqdb <- left_join(seqdb,revdb)
  seqdb[is.na(seqdb$Rev),]$Rev <- FALSE
  for(i in 1:gnum){
    genfocus[[settings$order[i]]] <- settings$focus
  }
}
#i# Auto mode has more complicated pairwise assessment.
if(settings$orient == "auto"){
  seqdb$Rev <- FALSE
  #># Will flip and order relative to the previous sequence
  #># Will need to propagate out from focal sequence rather than from genome 1
  #i# First, establish which sequence is orienting off which
  #i# Focus orients off itself!
  genfocus[[settings$focus]] <- settings$focus
  #i# Cycle until all added
  while(length(genfocus) < length(settings$order)){
    for(i in 1:gnum){
      genome <- settings$order[i]
      if(! genome %in% names(genfocus)){
        if(i < gnum && settings$order[i+1] %in% names(genfocus)){
          genfocus[[genome]] <- settings$order[i+1]
        }
        if(i > 1 && settings$order[i-1] %in% names(genfocus)){
          genfocus[[genome]] <- settings$order[i-1]
        }
        #i# Now update seqdb as appropriate using bestdb if added to genfocus.
        if(genome %in% names(genfocus)){
          for(i in c(1:nrow(seqdb))[seqdb$Genome == genome]){
            rowbest <- filter(bestdb,Genome==seqdb$Genome[i],SeqName==seqdb$SeqName[i],HitGenome==genfocus[[seqdb$Genome[i]]])
            if(nrow(rowbest) > 0){
              if(seqRev(seqdb,rowbest$HitGenome[1],rowbest$Hit[1]) & rowbest$Fwd[1] >= rowbest$Bwd[1]){
                seqdb$Rev[i] = TRUE      
              }
              if(! seqRev(seqdb,rowbest$HitGenome[1],rowbest$Hit[1]) & rowbest$Fwd[1] < rowbest$Bwd[1]){
                seqdb$Rev[i] = TRUE      
              }
            }
          }
        }
      }
    }
  }
}
#i# Update the text for each name
seqdb$text <- seqdb$SeqName
for(i in 1:nrow(seqdb)){
  if(seqdb$Rev[i]){
    seqdb$text[i] <- paste0(seqdb$text[i],"R")
  }
}

### ~ Establish ordering of sequences based on focus ~~~~~~~~~~~~~~~~~~~ ###
#!# Will want to add settings$seqsort=FILE to over-ride individual genome ordering
if(! settings$seqsort %in% c("auto","none","focus")){
  if(file.exists(settings$seqsort)){
    logWrite(paste0("Not yet implemented seqorder=FILE (",settings$seqsort,"): set to 'auto'."))
    settings$seqsort <- "auto"
  }
  else{
    logWrite(paste0("Could not recognise/find seqorder=",settings$seqsort," (auto/focus/none/FILE): set to 'auto'."))
    settings$seqsort <- "auto"
  }
}
#i# Set up the seqorder list. This will have focus and anything loaded from settings$seqorder
seqorder <- list()
seqorder[[settings$focus]] <- c(settings$seqorder, seqdb[seqdb$Genome==settings$focus & ! seqdb$SeqName %in% settings$seqorder,]$SeqName)
#i# seqsort=none will just order by sequence name
if(settings$seqsort == "none"){
  #?# Or might want to have options for sorting by length
  for(genome in settings$order[!settings$order %in% names(seqorder)]){
    seqorder[[genome]] <- seqdb$SeqName[seqdb$Genome == genome]
    logWrite(paste('#ORDER',genome,'sequences:',paste(seqorder[[genome]],collapse=", ")))
  }
}
#i# seqsort=focus will map everything to the focus species
if(settings$seqsort == "focus"){
  #i# posdb is only used here, for ordering sequences.
  posdb <- ungroup(bestdb) %>% filter(HitGenome==settings$focus)
  posdb$Order <- match(posdb$Hit,seqorder[[settings$focus]])
  posdb <- arrange(posdb,Order,HitPosMb)
  for(genome in settings$order[!settings$order %in% names(seqorder)]){
    seqorder[[genome]] <- c(posdb[posdb$Genome==genome,]$SeqName, seqdb$SeqName[! seqdb$SeqName %in% posdb$SeqName & seqdb$Genome == genome])
    logWrite(paste('#ORDER',genome,'sequences:',paste(seqorder[[genome]],collapse=", ")))
  }
}
#i# seqsort=auto will use the orientation focus (genfocus)
if(settings$seqsort == "auto"){
  while(length(seqorder) < length(settings$order)){
    for(genome in settings$order[!settings$order %in% names(seqorder)]){
      myfocus <- genfocus[[genome]]
      if(! myfocus %in% names(seqorder)){
        next
      }
      posdb <- ungroup(bestdb) %>% filter(HitGenome==myfocus)
      posdb$Order <- match(posdb$Hit,seqorder[[myfocus]])
      #i# Invert posdb for reversed Hit sequences
      for(i in 1:nrow(posdb)){
        if(seqRev(seqdb,myfocus,posdb$Hit[i])){
          mblen <- seqLen(seqdb,myfocus,posdb$Hit[i]) / 1e6
          posdb$HitPosMb[i] <- mblen - posdb$HitPosMb[i]
        }
      }
      posdb <- arrange(posdb,Order,HitPosMb)
      seqorder[[genome]] <- c(posdb[posdb$Genome==genome,]$SeqName, seqdb$SeqName[! seqdb$SeqName %in% posdb$SeqName & seqdb$Genome == genome])
      logWrite(paste('#ORDER',genome,'sequences:',paste(seqorder[[genome]],collapse=", ")))
    }
  }
}
#i# seqorder now has the order of the chromosomes for each genome

#?# This is correct but seems to be plotting out of order sometimes.


##### ======================== Setup Plots ======================== #####
#i# settings$order sets the vertical order
#i# seqorder[[genome]] sets the horizontal order for each genome
#i# regdb sets the linkages

### ~ Establish max genome size, gaps and padding ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite("Establish max genome size, gaps and padding")
gendb <- group_by(seqdb,Genome) %>% summarise(Genome=first(Genome),GenLen=sum(SeqLen),SeqNum=n())
xgap <- 0.01 * max(gendb$GenLen)
maxlen <- max(gendb$GenLen + (gendb$SeqNum - 1) * xgap)
if(settings$align == "justify"){
  gendb$xgap <- (maxlen - gendb$GenLen) / (gendb$SeqNum - 1)
}else{
  gendb$xgap <- xgap
}
gendb$xshift <- 0
if(settings$align == "right"){
  gendb$xshift <- maxlen - (gendb$GenLen + (gendb$SeqNum - 1) * xgap)
}
if(settings$align %in% c("centre","center")){
  gendb$xshift <- maxlen - (gendb$GenLen + (gendb$SeqNum - 1) * xgap) / 2
}
gendb$yshift <- 0 - (match(gendb$Genome,settings$order) - 1) * settings$ygap
rownames(gendb) <- gendb$Genome

### ~ Establish xshift and yshift for each sequence ~~~~~~~~~~~~~~~~~~~ ###
seqdb$xshift <- gendb[match(seqdb$Genome,gendb$Genome),]$xshift
seqdb$yshift <- gendb[match(seqdb$Genome,gendb$Genome),]$yshift
for(genome in settings$order){
  xshift <- gendb[genome,]$xshift
  for(seqname in seqorder[[genome]]){
    seqdb[seqdb$Genome==genome & seqdb$SeqName==seqname,]$xshift <- xshift
    xshift <- xshift + gendb[genome,]$xgap + seqLen(seqdb,genome,seqname)
  }
}

### ~ TIDK Telomeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite("TIDK Telomeres...")
#i# Combine with seqdb to get xshift and yshift
#i# Update teldb to have the correct xshift and yshift information, reversing where needed
teldb <- right_join(teldb,seqdb)
teldb <- teldb[! is.na(teldb$Pos),]
if(nrow(teldb)){
  teldb$RevPos <- teldb$SeqLen - teldb$Pos + 1
  teldb$RevStrand <- '+'
  teldb[teldb$Strand == '+',]$RevStrand <- '-'
  #i# Reverse the strand and position where needed for plotting
  teldb[teldb$Rev, ]$Strand <- teldb[teldb$Rev, ]$RevStrand
  teldb[teldb$Rev, ]$Pos <- teldb[teldb$Rev, ]$RevPos
  teldb$yshift <- teldb$yshift + 0.25
  teldb[teldb$Strand == '-',]$yshift <- teldb[teldb$Strand == '-',]$yshift + 0.5
}

### ~ Gaps and Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite("Gaps and Features...")
#i# Combine with seqdb to get xshift and yshift
#i# Update ftdb to have the correct xshift and yshift information, reversing where needed
ftdb <- right_join(ftdb,seqdb)
ftdb <- bind_rows(ftdb,right_join(gapdb,seqdb))
ftdb <- ftdb[! is.na(ftdb$Pos),]
if(nrow(ftdb)){
  ftdb$RevPos <- ftdb$SeqLen - ftdb$Pos + 1
  ftdb$RevStrand <- '.'
  ftdb[ftdb$Strand == '+',]$RevStrand <- '-'
  ftdb[ftdb$Strand == '-',]$RevStrand <- '+'
  #i# Reverse the strand and position where needed for plotting
  ftdb[ftdb$Rev, ]$Strand <- ftdb[ftdb$Rev, ]$RevStrand
  ftdb[ftdb$Rev, ]$Pos <- ftdb[ftdb$Rev, ]$RevPos
  ftdb$yshift <- ftdb$yshift + 0.5
  ftdb[ftdb$Strand == '-',]$yshift <- ftdb[ftdb$Strand == '-',]$yshift - 0.4
  ftdb[ftdb$Strand == '+',]$yshift <- ftdb[ftdb$Strand == '+',]$yshift + 0.4
}


##### ======================== Save data to Excel ======================== #####
### ~ Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(settings$writexl){
  if(endsWith(settings$basefile,"chromsyn")){
    outfile = paste0(settings$basefile,".xlsx")
  }else{
    outfile = paste0(settings$basefile,".chromsyn.xlsx")
  }
  X <- list(Genomes=gendb,Sequences=seqdb,BUSCO=busdb,Regions=regdb,Hits=hitdb,Best=bestdb,Positions=posdb)
  write_xlsx(
    x = X,
    path = outfile
  )
  logWrite(paste("#XLSX Tabs:",paste(names(X),collapse=",")))
  logWrite(paste("#SAVE","All chromsyn data output to",outfile))
}


##### ======================== Generate Plots ======================== #####
#i# settings$order sets the vertical order
#i# seqorder[[genome]] sets the horizontal order for each genome
#i# regdb sets the linkages

### ~ Generate chromosome plot ~~~~~~~~~~~~~~~~~~~ ###
chromSynPlot <- function(gendb,seqdb,regdb,linkages=c()){
  if(length(linkages)<1){
    linkages <- 1:(length(settings$order)-1)
  }
  scaling <- list(bp=1, kb=1e3, Mb=1e6, Gb=1e9)
  rescale <- scaling[[settings$scale]]
  
  #i# First, plot the chromosomes
  pD <- seqdb %>% mutate(xmin=xshift/rescale,xmax=(xshift+SeqLen)/rescale,ymin=yshift,ymax=yshift+1) %>% 
    arrange(yshift,xshift) 
  pD$ymod[odd(1:nrow(pD))] <- 1+settings$textshift
  pD$ymod[even(1:nrow(pD))] <- -settings$textshift
  pD$texty <- pD$ymin + pD$ymod
  plt <- ggplot() +
    geom_rect(data=pD, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Genome), color="black")
  if(settings$labels){
    plt <- plt +  annotate("text",label=pD$text,x=pD$xmin,y=pD$texty,colour="black",size=3 * settings$labelsize,hjust=0)
  }
  #i# Add genome labels
  pG <- group_by(pD,Genome) %>% summarise(Genome=first(Genome),ymax=min(ymax),xmin=min(xmin))
  pG$xmin <- pG$xmin - (0.25 * settings$ticks / rescale)
  plt <- plt +
    annotate("text",label=pG$Genome,x=pG$xmin,y=pG$ymax,colour="black",size=4 * settings$namesize,hjust=1)
  #i# Add tick marks
  tickx <- c()
  ticky <- c()
  for(i in 1:nrow(seqdb)){
    ntick <- as.integer(seqdb$SeqLen[i]/settings$ticks)
    tickx <- c(tickx, 0:ntick * settings$ticks + seqdb$xshift[i])
    ticky <- c(ticky, rep(seqdb$yshift[i],ntick+1))
  }
  tickx <- tickx / rescale
  # Top ticks
  pD <- data.frame(x=tickx,y=ticky,yend=ticky-(settings$textshift/2))
  plt <- plt + geom_segment(data=pD,mapping=aes(x=x,xend=x,y=y,yend=yend),colour="black")
  # Bottom ticks
  pD <- data.frame(x=tickx,y=ticky+1,yend=ticky+1+(settings$textshift/2))
  plt <- plt + geom_segment(data=pD,mapping=aes(x=x,xend=x,y=y,yend=yend),colour="black")
  # Features
  #?# Could have a feature type field and map colour by type?
  if(nrow(ftdb)){
    pD <- ftdb %>% mutate(xpos=(xshift+Pos)/rescale)
    plt <- plt + geom_point(data=pD,mapping=aes(x=xpos,y=yshift),colour=pD$Col,shape=pD$Shape)
  }
  # Telomeres
  pD <- seqdb %>% mutate(xmin=xshift/rescale,xmax=(xshift+SeqLen)/rescale,y=yshift+0.5)
  pD <- pD[ (pD$Tel5 & ! pD$Rev) | (pD$Tel3 & pD$Rev), ]
  if(nrow(pD)){
    plt <- plt + geom_point(data=pD,mapping=aes(x=xmin,y=y),colour="black")
  }
  pD <- seqdb %>% mutate(xmin=xshift/rescale,xmax=(xshift+SeqLen)/rescale,y=yshift+0.5)
  pD <- pD[ (pD$Tel5 & pD$Rev) | (pD$Tel3 & ! pD$Rev), ]
  if(nrow(pD)){
    plt <- plt + geom_point(data=pD,mapping=aes(x=xmax,y=y),colour="black")
  }
  # TIDK internal windows
  if(nrow(teldb)){
    pD <- teldb %>% mutate(xpos=(xshift+Pos)/rescale)
    plt <- plt + geom_point(data=pD,mapping=aes(x=xpos,y=yshift),colour="blue")
  }

  cat("Generating plot", file = stderr())

  #i# Then, plot the linkages
  vnum <- length(linkages)
  for(v in linkages){
    #i# Setup the pair
    genomea <- settings$order[v]
    genomeb <- settings$order[v+1]
    ya <- gendb$yshift[gendb$Genome==genomea]
    yb <- gendb$yshift[gendb$Genome==genomeb] + 1
    plotdb <- regdb %>% filter(Genome==genomea,HitGenome==genomeb)
    #i# Process the links for the pair
    for(i in 1:nrow(plotdb)){
      # Setup
      xa1 <- plotdb$Start[i]
      xa2 <- plotdb$End[i]
      xb1 <- plotdb$HitStart[i]
      xb2 <- plotdb$HitEnd[i]
      # Check for reversals
      fwd <- plotdb$Strand[i] == "+"
      if(seqRev(seqdb,plotdb$Genome[i],plotdb$SeqName[i])){
        tmp <- seqLen(seqdb,plotdb$Genome[i],plotdb$SeqName[i]) - xa1 + 1
        xa1 <- seqLen(seqdb,plotdb$Genome[i],plotdb$SeqName[i]) - xa2 + 1
        xa2 <- tmp
        fwd <- ! fwd
      }
      if(seqRev(seqdb,plotdb$HitGenome[i],plotdb$Hit[i])){
        tmp <- seqLen(seqdb,plotdb$HitGenome[i],plotdb$Hit[i]) - xb1 + 1
        xb1 <- seqLen(seqdb,plotdb$HitGenome[i],plotdb$Hit[i]) - xb2 + 1
        xb2 <- tmp
        fwd <- ! fwd
      }
      #i# Add xshift and scale
      xa1 <- (xa1 + seqdb$xshift[seqdb$Genome==genomea & seqdb$SeqName==plotdb$SeqName[i]]) / rescale
      xa2 <- (xa2 + seqdb$xshift[seqdb$Genome==genomea & seqdb$SeqName==plotdb$SeqName[i]]) / rescale
      xb1 <- (xb1 + seqdb$xshift[seqdb$Genome==genomeb & seqdb$SeqName==plotdb$Hit[i]]) / rescale
      xb2 <- (xb2 + seqdb$xshift[seqdb$Genome==genomeb & seqdb$SeqName==plotdb$Hit[i]]) / rescale
      
      #i# Plot links
      if(fwd){
        pD <- data.frame(x=c(xa1,xa2,xb2,xb1), y=c(ya,ya,yb,yb))
        plt <- plt + 
          geom_polygon(data=pD,mapping=aes(x=x, y=y),fill="steelblue", color=NA, alpha=settings$opacity) 
      }else{
        pD <- data.frame(x=c(xa1,xa2,xb1,xb2), y=c(ya,ya,yb,yb))
        plt <- plt + 
          geom_polygon(data=pD,mapping=aes(x=x, y=y),fill="indianred", color=NA, alpha=settings$opacity) 
      }
      
    }
    cat(paste0("...(",v,") ",round(100*which(v==linkages)/vnum,1),"%"), file = stderr())
  }
  cat("\n", file = stderr())    

  #i# Add the theme
  plt <- plt + theme_bw() + 
    xlab(paste0("Position (",settings$scale,")")) + ylab("") +
    theme(#legend.position = "None", 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )

  return(plt)  
}

### ~ Output one or more plots based on C stack limit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
linknum <- (length(settings$order)-1)
plotsplits <- 0
plt <- NA
pdfwidth <- settings$pdfwidth
pdfheight <- 1 + 2*length(settings$order)
if(settings$pdfheight > 0.1){
  pdfheight <- settings$pdfheight
}

while(is.na(plt) & plotsplits < linknum){
  plotsplits <- (plotsplits + 1)
  #i# Try generating plots until fragmented enough
  tryCatch(
    expr = {
      cat(paste(plotsplits,"plot(s)..."), file = stderr())
      for(px in 1:plotsplits){
        plt <- NA
        if(plotsplits > 1){
          linkn <- round(linknum/plotsplits) + 1
          linkages <- 1:linkn + (linkn-1) * (px - 1)
          linkages <- linkages[linkages <= linknum]
          if(settings$debug){
            cat(linkages)
          }
          plt <- chromSynPlot(gendb,seqdb,regdb,linkages)
          plotbase <- paste0(settings$basefile,"-",px)
        }else{
          plt <- chromSynPlot(gendb,seqdb,regdb)
          plotbase <- settings$basefile
        }
        print(plt)
        plotfile <- paste0(plotbase,".pdf")
        ggsave(plotfile,plot=plt,device="pdf",path=settings$plotdir,width=pdfwidth,height=pdfheight,scale=settings$pdfscale)
        logWrite(paste0('#GGSAVE Saved output plot to ',settings$plotdir,plotfile))

        plotfile <- paste0(plotbase,".png")
        ggsave(plotfile,plot=plt,device="png",path=settings$plotdir,width=pdfwidth,height=pdfheight,scale=settings$pdfscale)
        logWrite(paste0('#GGSAVE Saved output plot to ',settings$plotdir,plotfile))
        
        # pngfile <- paste0(plotbase,".png")
        # pngwidth <- as.integer(pdfwidth * settings$pngscale)
        # pngheight <- as.integer(pdfheight * settings$pngscale)
        # #png(pngfile,width=pngwidth,height=pngheight,pointsize=settings$pointsize)
        # png(pngfile,width=pdfwidth,height=pdfheight,units="in",pointsize=settings$pointsize,res=settings$pngscale)
        # print(plt)
        # dev.off()
        # logWrite(paste0('#PNG Saved output plot to ',settings$plotdir,pngfile))
        # 
        # pngfile <- paste0(plotbase,".jpeg")
        # jpeg(pngfile,width=pdfwidth,height=pdfheight,units="in",pointsize=settings$pointsize,res=settings$pngscale)
        # print(plt)
        # dev.off()
        # logWrite(paste0('#PNG Saved output plot to ',settings$plotdir,pngfile))
      }
    },
    error = function(e){
      print(e)
      logWrite(paste("#ERROR",e,"=> splitting plot into",plotsplits+1))
      plt <- NA
    }
  )
  
}


##### ======================== Tidy and Finish ======================== #####

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE ChromSyn.R finished.")
#quit("no",0)

