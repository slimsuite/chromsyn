########################################################
### ChromSyn Synteny Plot functions            ~~~~~ ###
### VERSION: 1.6.1                             ~~~~~ ###
### LAST EDIT: 10/09/24                        ~~~~~ ###
### AUTHORS: Richard Edwards 2022              ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au           ~~~~~ ###
### CITE: Edwards et al. bioRxiv 2022.04.22.489119 ~ ###
###       doi: 10.1101/2022.04.22.489119       ~~~~~ ###
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
# v0.8.0 : Added ypad setting to offset the syntenic block plotting a little for clarity. Added fill for open symbols.
# v0.9.0 : Added chromfill setting to control the colouring of chromosomes by Genome, Type or Col.
# v0.9.1 : Added ftsize=NUM setting to control the size of telomere and feature points. Fixed rare seqorder bug and R default issues.
# v0.9.2 : Added scale=pc option to rescale all assemblies to be proportions of total length. Increased default ftsize.
# v0.9.3 : Fixed feature-plotting bug when chromosome fill set.
# v0.9.4 : Fixed regdata bugs. Moved orphan reduction to follow regdata incorporation.
# v0.9.5 : Changed default feature plot to a black square with white fill. Enabled missing subsets of files.
# v1.0.0 : Changed default plotting to remove legend. First release, corresponding to Myrtle Rust bioRxiv paper.
# v1.1.0 : Added duplicated=T/F option for plotting Duplicated BUSCO genes and added features to xlsx output. Added Duplicated Gap output.
# v1.1.1 : General bug fixes. Added minftlen=INT.
# v1.1.2 : Enforced loading sequence names as character.
# v1.1.3 : Fixed some bugs when there are missing data in places.
# v1.1.4 : Fixed bug with reverse chrom not plotting Duplicated genes.
# v1.1.5 : Fixed Duplicated BUSCO loading bug introduced by v1.1.4.
# v1.1.6 : Fixed gap-free plotting bug introduced by v1.1.5.
# v1.1.7 : Fixed bug with TSV input for TIDK. (Expected CSV.)
# v1.2.0 : Added reading and parsing of CN table for mapping features onto. Added restrict=LIST.
# v1.2.1 : Fixed some input reading bugs.
# v1.3.0 : Added regmirror=T/F function to allow unidirectional regdata and fixed order bug without BUSCO.
# v1.3.1 : Added BuscoIDList output to Regions sheet of Excel output. Fixed bug with orient=none.
# v1.3.2 : Fixed bug with scale=pc mode.
# v1.4.0 : Added plotting of repeats with a -2 shape as inverted triangles near centre. Fixed minor list argument parsing bug.
# v1.5.0 : Added more complex plotting of gaps based on SynBad ratings and qcmode=T/F. Added maxregions=INT setting and regdata collapse to dynamically sets the min region length.
# v1.6.0 : Moved linkage plotting to behind features.
# v1.6.1 : Fixed cndata bug when not all CN categories present.
version = "v1.6.1"

####################################### ::: USAGE ::: ############################################
# Example use:
# Rscript chromsyn.R [sequences=FOFN] [busco=FOFN] [focus=X] [seqorder=LIST] [order=LIST] [basefile=FILE] [minregion=INTbp] [align=X]
# : sequences=FOFN = File of PREFIX FILE with sequence names and lengths (name & length, or SeqName & SeqLen fields) [sequences.fofn]
# : busco=FOFN = File of PREFIX FILE with full BUSCO table results. Used to identify orthologous regions. [busco.fofn]
# : tidk=FOFN = Optional file of PREFIX FILE with TIDK search results. [tidk.fofn]
# : gaps=FOFN = Optional file of PREFIX FILE with gap positions. [gaps.fofn]
# : ft=FOFN = Optional file of PREFIX FILE with custom features. [ft.fofn]
# : regdata=TSV = File of Genome, HitGenome, SeqName, Start, End, Strand, Hit, HitStart, HitEnd
# : regmirror=T/F = Whether to mirror and collapse the regdata input [True]
# : cndata=TSV = File of Genome, SeqName, Start, End, CN (Other fields OK)
# : focus=X = If given will orient all chromosomes to this assembly
# : orient=X = Mode for sequence orientation (none/focus/auto)
# : seqsort=none/focus/auto/FILE = Optional ordering strategy for other assemblies [auto]
# : seqorder=LIST = Optional ordering of the chromsomes for the focal assembly
# : restrict=LIST = List of sequence names to restrict analysis to. Will match to any genome.
# : order=LIST = File containing the Prefixes to include in vertical order. If missing will use sequences=FOFN.
# : chromfill=X = Sequences table field to use for setting the colouring of chromosomes (e.g. Genome, SeqName, Type or Col) [Genome]
# : qcmode=T/F = QC mode that will hide gaps rated as "good" by SynBad [False]
# : synbad=FILE = File of custom shapes and colours for SynBad gap ratings []
# : runpath=PATH = Run ChromSyn in this path (will look for inputs etc. here) [./]
# : basefile=FILE = Prefix for outputs [chromsyn]
# : plotdir=PATH = output path for graphics
# : minlen=INT = minimum length for a chromosome/scaffold to be included in synteny blocks/plots [0]
# : minregion=INT = minimum length for mapped regions to be included in plots [50000]
# : maxregions=INT = maximum number of mapped regions to be included in plots [0]
# : minbusco=INT = minimum number of BUSCO genes to be included in Syntenic block [1]
# : minftlen=INT = minimum length of feature to be plotted [1]
# : maxskip=0 = maximum number of BUSCO genes to skip and still be a syntenic block [0]
# : orphans=T/F = whether to include scaffolds that have no BUSCO genes [True]
# : tidkcutoff=INT = TIDK count cutoff for identifying a telomere [50]
# : align=X = alignment strategy for plotting chromosomes (left/right/centre/justify) [justify]
# : ygap=INT = vertical gap between chromosomes [4]
# : ypad=NUM = proportion of ygap to extend synteny blocks before linking [0.1]
# : ybleed=NUM = proportion of chromosome to bleed synteny block plotting into [0]
# : scale=X = units in basepairs for setting the x-axis scale (bp/kb/Mb/Gb/pc) [Mb]
# : textshift=NUM = offset for printing chromosome names [0.3]
# : ticks=INT = distance between tickmarks [5e7]
# : pdfwidth=NUM = PDF width [20]
# : pdfheight=NUM = over-ride for standard calculated PDF height [0]
# : pdfscale=NUM = over-ride for PDF output scale [1]
# : ftsize=NUM setting to control the size of telomere and feature points. [1]
# : namesize=NUM = scaling factor for the Genome names in PDF plots [1]
# : labelsize=NUM = scaling factor for the chromosome names in PDF plots [1]
# : labels=T/F = whether to print chromosome name labels (True) or legend (False) [TRUE]
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
#!# Need to test no loading of BUSCO data if regdata=TSV is given.
#!# Add output of synteny block table.
#!# Add file checks and error messages!
#!# Add scale bar for tick marks.
#!# Add compilation of regdata synteny blocks like BUSCO blocks. Add a second regdata minlen cutoff (regfilter=X).
#!# Add different gap modes (+ or line)
#!# Add depth input for coloring of chromosome scaffolds.
#!# Add more in-flight checks for data integrity.
# : labels=X = Sequence table field to use as labels in plots [Label -> SeqName]
# : genlabels=FILE = Option file of Genome Label to use in plots [labels.txt]
#?# Update restrict to fill in unidrectional best hits if no restriction given for that Genome.
#!# Add bestlists ouput that will generate a text file with the best hits in other genomes for each sequence. (For future restrict)
#!# Add restrictmode=expand to include the hits above a % (of total) synteny to restrict? (Accounting for focus and order.)
#!# Add a helper script to compiles a features table from common sources, such as rRNA predictions.
#!# Add output of the BUSCO IDs within each synteny block .
#!# Add centromeres=FOFN input and special plotting. 

####################################### ::: SETUP ::: ############################################
### ~ Commandline arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
defaults = list(busco="busco.fofn",sequences="sequences.fofn",order="",
                regdata="",cndata="",restrict="",regmirror=TRUE,
                minregion=50000,maxregions=0,
                align="justify",ygap=4,ypad=0.1,ybleed=0.0,synbad="",
                seqorder="",orient="auto",seqsort="auto",chromfill="Genome",qcmode=FALSE,
                #pngwidth=1200,pngheight=900,
                #pointsize=24,pngscale=100,
                plotdir="./",duplicated=TRUE,minftlen=1,
                tidkcutoff=50,tidk="tidk.fofn",gaps="gaps.fofn",ft="ft.fofn",
                minbusco=1,maxskip=0,orphans=TRUE,minlen=0,opacity=0.3,ftsize=1,
                pdfwidth=20,pdfheight=0,pdfscale=1,namesize=1,labelsize=1,labels=TRUE,
                scale = "Mb",textshift = 0.3,ticks=5e7,rscript=TRUE,
                basefile="chromsyn",focus="",debug=FALSE,dev=FALSE,
                rdir="",runpath="",
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
    settings[[cmdv[1]]] = ""    
  }
}
#i# integer parameters
for(cmd in c("pngwidth","pngheight","pointsize","minregion","maxregions","ygap","minbusco","maxskip","minlen","pngscale","tidkcutoff","minftlen")){
  settings[[cmd]] = as.integer(settings[[cmd]])
}
#i# other numeric parameter
for(cmd in c("textshift","ticks","pdfwidth","pdfheight","pdfscale","namesize","labelsize","opacity","ypad","ybleed","ftsize")){
  settings[[cmd]] = as.numeric(settings[[cmd]])
}
#i# adjust parameters where needed
settings$ftsize <- settings$ftsize * 2
if(settings$opacity < 0.1){
  settings$opacity <- 0.1
}
if(settings$opacity > 1){
  settings$opacity <- 1
}
#i# list parameters
for(cmd in c("order","seqorder","restrict")){
  if(sum(grep(",",settings[[cmd]],fixed=TRUE)) > 0){
    settings[[cmd]] = strsplit(settings[[cmd]],',',TRUE)[[1]]
  }else{
    settings[[cmd]] <- settings[[cmd]]
  }
}
#i# logical parameters
for(cmd in c("debug","dev","orphans","labels","duplicated","regmirror","qcmode")){
  settings[[cmd]] = as.logical(settings[[cmd]])
}

#i# Set warnings based on debug
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
if(! settings$runpath == ""){
  setwd(settings$runpath)
}
logWrite(paste("#RCODE ChromSyn.R:",version))
logWrite(paste("#PATH Running from:",getwd()))
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

### ~ Load R libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
if(settings$rdir == ""){
  settings$rdir <- getScriptPath()
}
sfile <- paste0(settings$rdir,"/rje_load.R")
logWrite(sfile)
source(sfile)

####################################### ::: FUNCTIONS ::: ############################################

##### ======================== Loading data functions ======================== #####

### ~ Load FOFN File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into FOFN table (genome, file)
#i# fofndb = fileTable(filename,delimit="\t")
fileTable <- function(filename,delimit=" ",genomes=c()){
  regdb <- read.table(filename,fill=TRUE,sep=delimit,header=FALSE,row.names = NULL,quote="\"",comment.char="",stringsAsFactors=FALSE) %>%
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
#i# seqdb <- seqTable(filename,delimit="\t")
seqTable <- function(genome,filename,delimit="\t"){
  if(! file.exists(filename)){ return(data.frame()) }
  seqdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="",stringsAsFactors=FALSE)
  seqdb$Genome <- genome
  if("name" %in% colnames(seqdb)){
    seqdb <- seqdb %>% rename(SeqName=name,SeqLen=length)
  }
  seqdb$SeqName <- as.character(seqdb$SeqName)
  if(! "Tel5" %in% colnames(seqdb)){
    seqdb$Tel5 <- FALSE
    seqdb$Tel3 <- FALSE
  }
  seqdb$Tel5 <- as.logical(seqdb$Tel5)
  seqdb$Tel5[is.na(seqdb$Tel5)] <- FALSE
  seqdb$Tel3 <- as.logical(seqdb$Tel3)
  seqdb$Tel3[is.na(seqdb$Tel3)] <- FALSE
  #i# Check for settings$chromfill
  keepfields <- c("Genome","SeqName","SeqLen","Tel5","Tel3")
  if(settings$chromfill %in% colnames(seqdb) & ! settings$chromfill %in% keepfields){
    keepfields <- c(keepfields,settings$chromfill)
  }
  seqdb <- seqdb %>% select(all_of(keepfields))
  logWrite(paste('#SEQS',nrow(seqdb),genome,"sequences loaded from",filename))
  seqdb <- seqdb[seqdb$SeqLen >= settings$minlen,]
  logWrite(paste('#SEQS',nrow(seqdb),genome,"sequences meet minlen cutoff of",settings$minlen,"bp"))
  if(length(settings$restrict) > 1){
    seqdb <- seqdb[seqdb$SeqName %in% settings$restrict,]
    logWrite(paste('#SEQS',nrow(seqdb),genome,"sequences found in filter=LIST."))
  }
  return(seqdb)
}

### ~ Load BUSCO File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load BUSCO table into tibble
#i# buscodb = buscoTable(filename)
#i# BUSCOv3 headings:
v3head = c("BuscoID","Status","Contig","Start","End","Score","Length")
#i# BUSCOv5 headings:
v5head = c("BuscoID","Status","Contig","Start","End","Strand","Score","Length","OrthoDBURL","Description")
#i# NOTE: URL and Description not always there.
buscoTable <- function(genome,filename,seqnames=c(),btype="Complete"){
  if(! file.exists(filename)){ return(data.frame()) }
  buscodb = read.table(filename,fill=TRUE,row.names = NULL,sep="\t",quote="",stringsAsFactors=FALSE)
  if(ncol(buscodb) > length(v5head)){
    logWrite(paste("#BUSCOV Dropping",ncol(buscodb) - length(v5head),"input field(s), not found in v5 format"))
    buscodb <- buscodb[,1:length(v5head)]
  }
  if(ncol(buscodb) > length(v3head)){
    logWrite(paste("#BUSCOV BUSCO v5 format"))
    colnames(buscodb) = v5head[1:ncol(buscodb)]
  }else{
    logWrite(paste("#BUSCOV BUSCO v3 format"))
    colnames(buscodb) = v3head
  }
  buscodb$Contig = as.character(buscodb$Contig)
  buscodb <- buscodb[buscodb$Status == btype,]
  logWrite(paste('#BUSCO',nrow(buscodb),genome,btype,"BUSCO genes loaded from",filename))
  if(nrow(buscodb) > 0){
    buscodb$Genome <- genome
    buscodb$Start <- as.integer(buscodb$Start)
    buscodb$End <- as.integer(buscodb$End)
    buscodb <- buscodb[! is.na(buscodb$Start) & ! is.na(buscodb$End), ]
  }else{
    buscodb <- buscodb %>% mutate(Genome="")
  }
  buscodb <- buscodb %>% rename(SeqName=Contig) %>% select(Genome,SeqName,Start,End,Strand,BuscoID)
  #i# Filter to sequences if given
  if(length(seqnames) > 0){
    buscodb <- buscodb[buscodb$SeqName %in% seqnames,]
    logWrite(paste('#BUSCO',nrow(buscodb),genome,btype,"BUSCO genes following filtering to",length(seqnames),"sequences."))
  }
  return(buscodb)
}

### ~ Load Regions File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into Regions table (Genome, HitGenome, SeqName, Start, End, Strand, Hit, HitStart, HitEnd)
#i# regdb = regTable(filename,delimit="\t")
regTable <- function(filename,delimit="\t"){
  if(! file.exists(filename)){ return(data.frame()) }
  regdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="",stringsAsFactors=FALSE)
  if(! "Length" %in% colnames(regdb)){
    regdb <- mutate(regdb,Length=End-Start+1)
  }
  if(! "HitLength" %in% colnames(regdb)){
    regdb <- mutate(regdb,HitLength=HitEnd-HitStart+1)
  }
  regdb <- select(regdb,Genome, HitGenome, SeqName, Start, End, Strand, Hit, HitStart, HitEnd, Length, HitLength)
  regdb$SeqName <- as.character(regdb$SeqName)
  regdb$Hit <- as.character(regdb$Hit)
  #?# Check the fields
  logWrite(paste('#REGION',nrow(regdb),"linked regions loaded from",filename))
  #i# Mirror if appropriate
  if(settings$regmirror){
    mirrordb <- regdb %>% select(HitGenome, Genome, Hit, HitStart, HitEnd, Strand, SeqName, Start, End, HitLength, Length)
    colnames(mirrordb) <- colnames(regdb)
    regdb <- unique(bind_rows(regdb,mirrordb)) %>% filter(Genome != HitGenome)
    logWrite(paste('#REGION',nrow(regdb),"unique linked regions following region mirroring (regmirror=TRUE)"))
    # Collapse adjacent synteny regions
    regdb <- arrange(regdb, Genome,HitGenome,SeqName,Start,End)
    #i# Add counter for each Genome (no need to restart numbering)
    regdb$Counter <- 1:nrow(regdb)
    regdb <- arrange(regdb,Genome,HitGenome,Hit,HitStart,HitEnd)
    regdb$HitCounter <- 1:nrow(regdb)
    regdb <- arrange(regdb,Genome,HitGenome,SeqName,Start,End)
    #i# - Add synteny block ID
    regdb$Block <- 1
    #i# - compress by block
    progv <- as.integer(1:20 * (nrow(regdb)/20))
    logWrite(paste('Generating synteny blocks, skipping upto (maxskip =)',settings$maxskip,'regions...'))
    noblock <- settings$maxskip + 2
    if(settings$minregion < 0){ noblock <- 0 }
    for(i in 2:nrow(regdb)){
      j <- i - 1
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
    regdb <- regdb %>% group_by(Block) %>% 
      summarise(Genome=first(Genome),HitGenome=first(HitGenome),SeqName=first(SeqName),Start=min(Start),End=max(End),Strand=first(Strand),Hit=first(Hit),HitStart=min(HitStart),HitEnd=max(HitEnd),BuscoID=n(),BuscoIDList=toString(BuscoID)) %>% 
      select(-Block) %>%
      mutate(Length=End-Start+1,HitLength=HitEnd-HitStart+1)
  }
  return(regdb)
}

### ~ Load TIDK File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load delimited file into TIDK table (Genome, SeqName, SeqWin, Tel5, Tel3, TelSeq)
#i# Headers: id,window,forward_repeat_number,reverse_repeat_number,telomeric_repeat
#i# tidkdb = tidkTable(genome,filename,delimit=",")
tidkTable <- function(genome,filename,delimit=","){
  if(! file.exists(filename)){ return(data.frame()) }
  if(endsWith(filename,".tsv") | endsWith(filename,".tdt")){
    delimit <- "\t"
  }
  tidkdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="",stringsAsFactors=FALSE)
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
  if(nrow(tidkdb) > 0){
    tidkdb$SeqName <- as.character(tidkdb$SeqName)
  }
  #?# Check the fields: Genome, SeqName, Pos, Strand
  logWrite(paste('#TIDK',nrow(tidkdb),"TIDK telomere windows loaded from",filename))
  return(tidkdb)
}

### ~ Load Feature/Gap File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Set and load SynBad gap rating colours
#!# To be added
goodgap <- c("Syn","Aln","Span","Ins","Long","Div",'HiC','DupHiC')
gapcol <- rep("blue",length(goodgap))
gapshape <- rep(3,length(goodgap))
neutralgap <- c("Null","Term")
gapcol <- c(gapcol, rep("black",length(neutralgap)))
gapshape <- c(gapshape, rep(3,length(neutralgap)))
qcgap <- c("Dup","Inv","InvFix","InvDupFix","Brk","Tran","Frag","InvBrk")
gapcol <- c(gapcol, c("red","darkred","darkred","red","red","red","red","red"))
gapshape <- c(gapshape, c(3,4,13,13,3,3,3,3))
synbaddb <- tibble(SynBad=c(goodgap, neutralgap, qcgap),
                   Col=gapcol,Shape=gapshape)
if(settings$dev){
  #i# Going to trial adding a second row for the inverted triangles. Will both be plotted?
  qcadd <- tibble(SynBad=c("Inv","InvFix","InvDupFix","Brk","Tran"),
                  Col=c("darkred","darkred","darkred","red","red"),
                  Shape=c(-2,-2,-2,-2,-2))
  synbaddb <- bind_rows(synbaddb,qcadd)
}
if(file.exists(settings$synbad)){
  tmpdb <- loadTable(settings$synbad) %>% select(SynBad,Col,Shape)
  if(ncol(tmpdb) == 3){
    synbaddb <- tmpdb
    logWrite(paste("#SYNBAD Loading SynBad colour scheme from",settings$synbad))
  }
}

#i# Load delimited file into feature or gaps (Genome, Seqname, Pos, Strand, Col, Shape)
#i# ftdb = ftTable(genome,filename,delimit=",")
ftTable <- function(genome,filename,colour="white",shape=22){
  if(! file.exists(filename)){ return(data.frame()) }
  delimit <- "\t"
  if(endsWith(filename,".csv")){
    delimit <- ","
  }
  ftdb <- read.table(filename,fill=TRUE,sep=delimit,header=TRUE,row.names = NULL,quote="\"",comment.char="",stringsAsFactors=FALSE)
  if(nrow(ftdb) > 0){
    for(ftfield in c("SeqName","Pos","Strand","Col","Shape","Start","End")){
      if(tolower(ftfield) %in% colnames(ftdb)){
        ftdb[[ftfield]] <- ftdb[[tolower(ftfield)]]
      }
    }
    if(settings$minftlen > 1 && "Start" %in% colnames(ftdb) && "End" %in% colnames(ftdb)){
      prex <- nrow(ftdb)
      ftdb <- ftdb[(ftdb$End - ftdb$Start + 1) >= settings$minftlen,]
      logWrite(paste(prex,"features reduced to",nrow(ftdb),"on minftlen (length <",settings$minftlen,"bp)"))
    }
  }
  if(nrow(ftdb) > 0){
    ftdb$Genome <- genome
  }else{
    logWrite(paste('#FT',nrow(ftdb),"features loaded from",filename))
    return(data.frame())
  }
  
  if(! "Pos" %in% colnames(ftdb)){
    ftdb <- mutate(ftdb,Pos=(End+Start)/2)
  }
  if(! "Strand" %in% colnames(ftdb)){
    ftdb$Strand <- "."  # Unstranded
  }
  for(oddstrand in c("+/-","+-","-+")){
    if(oddstrand %in% ftdb$Strand){
      ftdb[ftdb$Strand == oddstrand,]$Strand <- "."
    }
  }
  # Special SynBad Gap colouration
  #i# When loading from the SynBad table, a colour of NA or shape of 0 will be filtered
  if(! "Col" %in% colnames(ftdb) & "SynBad" %in% colnames(ftdb) & nrow(synbaddb) > 0){
    # Special SynBad Gap colouration
    ftdb <- left_join(ftdb,synbaddb) %>% 
      mutate(Col=if_else(is.na(Col),"NA",Col)) %>%
      filter(Col != "NA", Shape != 0)
    # Add the inverted triangles
    strandme <- ftdb %>% filter(Shape == -2)
    if(nrow(strandme) > 0){
      ftdb <- bind_rows(ftdb %>% filter(Shape != -2),strandme %>% mutate(Strand="+"),strandme %>% mutate(Strand="-"))
    }
  }
  if(! "Col" %in% colnames(ftdb) & "SynBad" %in% colnames(ftdb)){
    ftdb <- ftdb %>%
      mutate(Col="darkred", Shape=3) %>%
    #i# Syn, Aln and Span are good (blue +)
    #i# Ins, Long, Div, also good (blue +)
      mutate(Col=if_else(SynBad %in% c("Syn","Aln","Span","Ins","Long","Div"),"blue",Col)) %>%
    #i# Null and Term are neutral (black +)
      mutate(Col=if_else(SynBad %in% c("Null","Term"),"black",Col)) %>%
    #i# Dupl is OK (red +)
      mutate(Col=if_else(SynBad %in% c("Dupl"),"red",Col)) %>%
    #i# Inv is bad (dark red X)
      mutate(Shape=if_else(SynBad %in% c("Inv"),4,Shape)) %>%
    #i# InvFix is special bad (dark red X in circle) [Or a blue X?]
      mutate(Shape=if_else(SynBad %in% c("InvFix"),13,Shape)) %>%
    #i# InvDupFix is special bad (red X in circle) [Or a red X?]
      mutate(Shape=if_else(SynBad %in% c("InvDupFix"),13,Shape)) %>%
      mutate(Col=if_else(SynBad %in% c("InvDupFix"),"red",Col))
    #i# All others are bad (dark red +)
    logWrite("#COL Added SynBad gap colours and shapes.")
  }
  #i# Special QC mode to drop Good rated gaps. Could also think of adding -2 for Fix rather than having squares?
  if(settings$qcmode & "SynBad" %in% colnames(ftdb)){
    dropn <- sum(ftdb$SynBad %in% c("Syn","Aln","Span","Ins","Long","Div"))
    ftdb <- ftdb %>% filter(! SynBad %in% c("Syn","Aln","Span","Ins","Long","Div"))
    logWrite(paste("#QC QC Mode: dropped",dropn,"`Good` rated SynBad gaps."))
  }

  if(! "Col" %in% colnames(ftdb)){
    ftdb$Col <- colour  # Defaults to white
  }
  if(! "Shape" %in% colnames(ftdb)){
    ftdb$Shape <- shape  # Defaults to square (3 = + for gaps)
  }
  ftdb$Fill <- ftdb$Col
  openrow <- ftdb$Shape %in% c(-1,-2,21:25)
  if(sum(openrow)){
    ftdb[openrow,]$Col <- "black"
  }
  if(! "Feature" %in% colnames(ftdb)){
    ftdb <- mutate(ftdb,Feature="Feature")
  }
  ftdb <- select(ftdb,Genome, SeqName, Pos, Strand, Feature, Col, Fill, Shape)
  if(nrow(ftdb) > 0){
    ftdb$SeqName <- as.character(ftdb$SeqName)
  }
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
  if(is.na(seqname)){
    logWrite(paste0('Warning: problem with missing seqname for ',genome))
    if(! file.exists(settings$busco)){
      logWrite(paste0('No BUSCO data. Check regdata and regmirror settings.'))
    }
    return(FALSE)
  }
  seqrev <- seqdb %>% filter(Genome==genome,SeqName==seqname) %>% select(Rev)
  if(length(seqrev$Rev) < 1 | is.na(seqrev$Rev[1])){
    logWrite(paste0('Warning: problem with SeqRev(',genome,",",seqname,")"))
    return(FALSE)
  }
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
  msg <- paste("Cannot find sequence FOFN file:",seqfofn)
  if(settings$rscript){
    logWrite(msg)
    quit("no",2)  
  }else{
    stop(msg)
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
#i# Optional extra tables: TIDK, gaps, ft
filechecks <- list(tidk="TIDK", gaps="Assembly gap", ft="features table", busco="BUSCO full")
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
#i# BUSCO tables
buscofofn = settings$busco
logWrite(paste("BUSCO FOFN file:",buscofofn))
if(! file.exists(buscofofn) & settings$regdata == ""){
  logWrite(paste("No BUSCO FOFN data file found."))
  if(settings$rscript){
    quit("no",2)  
  }else{
    stop(paste("Cannot find BUSCO FOFN file:",buscofofn))
  }
}
logWrite('#RCODE Setup complete.')

##### ======================== Load data ======================== #####

### ~ Load sequence and BUSCO data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Load file names
seqfiles <- fileTable(settings$sequences) %>% rename(SeqFile=Filename)
gendb <- seqfiles
if(file.exists(buscofofn)){
  buscofiles <- fileTable(settings$busco,genomes=seqfiles$Genome) %>% rename(BUSCO=Filename)
  #i# Combine into table
  gendb <- inner_join(seqfiles,buscofiles)
}
if(length(settings$order)<1 | (length(settings$order)==1 & settings$order[1] == "")){
  settings$order <- seqfiles$Genome
  logWrite(paste("Genomes (order=LIST):",paste(settings$order,collapse=", ")))
}
# - TIDK
if(file.exists(settings$tidk)){
  tidkfiles <- fileTable(settings$tidk,genomes=seqfiles$Genome) %>% rename(TIDK=Filename)
  gendb <- left_join(gendb,tidkfiles)
}else{
  gendb$TIDK <- NA
}
# - gaps
if(file.exists(settings$gaps)){
  gapfiles <- fileTable(settings$gaps,genomes=seqfiles$Genome) %>% rename(gaps=Filename)
  gendb <- left_join(gendb,gapfiles)
}else{
  gendb$gaps <- NA
}
# - features
if(file.exists(settings$ft)){
  ftfiles <- fileTable(settings$ft,genomes=seqfiles$Genome) %>% rename(features=Filename)
  gendb <- left_join(gendb,ftfiles)
}else{
  gendb$features <- NA
}

#i# Update sequence order
settings$order <- settings$order[settings$order %in% gendb$Genome]
gendb <- gendb %>% filter(Genome %in% settings$order)
rownames(gendb) <- gendb$Genome
logWrite(paste("#GENOME ",nrow(gendb),"genomes:",paste(settings$order,collapse=", ")))
if(settings$debug){ 
  logWrite(paste("#GENDB ",nrow(gendb),"genomes:",paste(gendb$Genome,collapse=", "))) 
  print(gendb)
}
genomes <- settings$order
#i# Load actual sequence and busco data
seqdb <- data.frame()
busdb <- data.frame()
dupdb <- data.frame()
teldb <- data.frame()
gapdb <- data.frame()
ftdb <- data.frame()
for(genome in genomes){
  logWrite(paste0(genome,"..."))
  #i# Sequences
  filename <- gendb[gendb$Genome == genome,]$SeqFile[1]
  if(settings$debug){ logWrite(paste0(filename,": ",file.exists(filename))) }
  newseqdb <- seqTable(genome,filename)
  if(nrow(seqdb) > 0){
    #i# Check for inconsistent presence of settings$chromfill field
    if(settings$chromfill %in% colnames(seqdb) & ! settings$chromfill %in% colnames(newseqdb)){
      newseqdb[[settings$chromfill]] <- "NA"
    }
    if(! settings$chromfill %in% colnames(seqdb) & settings$chromfill %in% colnames(newseqdb)){
      seqdb[[settings$chromfill]] <- "NA"
    }
    seqdb <- bind_rows(seqdb,newseqdb)
  }else{
    seqdb <- newseqdb
  }
  if("BUSCO" %in% colnames(gendb)){
    #i# BUSCO
    filename <- gendb[genome,"BUSCO"]
    if(nrow(busdb) > 0){
      newbusdb <- buscoTable(genome,filename,newseqdb$SeqName)
      if(nrow(newbusdb) > 0){
        busdb <- bind_rows(busdb,newbusdb)
      }
    }else{
      busdb <- buscoTable(genome,filename,newseqdb$SeqName)
    }
    #i# BUSCO Duplicated
    filename <- gendb[genome,"BUSCO"]
    if(nrow(dupdb) > 0){
      newbusdb <- buscoTable(genome,filename,newseqdb$SeqName,"Duplicated")
      if(nrow(newbusdb) > 0){
        dupdb <- bind_rows(dupdb,newbusdb)
      }
    }else{
      dupdb <- buscoTable(genome,filename,newseqdb$SeqName,"Duplicated")
    }
  }
  #i# TIDK
  filename <- gendb[genome,"TIDK"]
  if(! is.na(filename)){
    newdb <- tidkTable(genome,filename)
    if(nrow(teldb) > 0){
      if(nrow(newdb)){
        teldb <- bind_rows(teldb,newdb)
      }
    }else{
      teldb <- tidkTable(genome,filename)
    }
  }
  #i# Gaps
  filename <- gendb[genome,"gaps"]
  if(! is.na(filename)){
    newdb <- ftTable(genome,filename,colour="darkred",shape=3)
    if(nrow(gapdb) > 0){
      if(nrow(newdb)){
        gapdb <- bind_rows(gapdb,newdb)
      }
    }else{
      gapdb <- newdb
    }
  }
  #i# Features
  filename <- gendb[genome,"features"]
  if(! is.na(filename)){
    newdb <- ftTable(genome,filename)
    if(nrow(ftdb) > 0){
      if(nrow(newdb)){
        ftdb <- bind_rows(ftdb,newdb)
      }
    }else{
      ftdb <- newdb
    }
  }
}
#i# Duplicated table update
if(nrow(dupdb)){
  dupdb$Pos <- (dupdb$End + dupdb$Start) / 2
  dupdb <- select(dupdb, Genome, SeqName, Pos, Strand, BuscoID) %>% rename(Feature=BuscoID)
}
if(nrow(dupdb)){
  dupdb$Col <- "black"
  dupdb$Fill <- "yellow"
  dupdb$Shape <- -1
}

#i# Summary
logWrite(paste("#SEQS ",nrow(seqdb),"sequences loaded in total."))
logWrite(paste("#BUSCO ",nrow(busdb),"Complete BUSCO genes loaded in total."))
logWrite(paste("#DUPL ",nrow(dupdb),"Duplicated BUSCO genes loaded in total."))
logWrite(paste("#TIDK ",nrow(teldb),"TIDK telomere windows loaded in total."))
logWrite(paste("#FT ",nrow(ftdb),"features loaded in total."))
logWrite(paste("#GAPS ",nrow(gapdb),"assembly gaps loaded in total."))

#i# Filter BUSCO to loaded sequences
if(nrow(busdb)){
  busdb <- filter(busdb,SeqName %in% seqdb$SeqName)
}
logWrite(paste("#BUSCO ",nrow(busdb),"BUSCO genes for loaded sequences."))

#i# Check for appropriate fill field:
if(! settings$chromfill %in% colnames(seqdb)){
  newfill <- "Genome"
  if("Col" %in% colnames(seqdb)){
    newfill <- "Col"
  }
  logWrite(paste0("Cannot find chromfill field '",settings$chromfill,"' in sequence input tables: setting to chromfill=",newfill))
  settings$chromfill <- newfill
}
if(settings$chromfill == "Col" & "NA" %in% seqdb[[settings$chromfill]]){
  seqdb[seqdb[[settings$chromfill]] == "NA",][[settings$chromfill]] <- "white"
}


##### ======================== Compile data ======================== #####
regdb <- data.frame(Genome=c(),HitGenome=c(),SeqName=c(),Start=c(),End=c(),Strand=c(),Hit=c(),HitStart=c(),HitEnd=c(),BuscoID=c(),Length=c(),HitLength=c())
if(nrow(busdb)){
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
  if(settings$minregion < 0){ noblock <- 0 }
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
  regdb <- regdb %>% group_by(Block) %>% 
    summarise(Genome=first(Genome),HitGenome=first(HitGenome),SeqName=first(SeqName),Start=min(Start),End=max(End),Strand=first(Strand),Hit=first(Hit),HitStart=min(HitStart),HitEnd=max(HitEnd),BuscoID=n(),BuscoIDList=toString(BuscoID)) %>% 
    select(-Block) %>%
    mutate(Length=End-Start+1,HitLength=HitEnd-HitStart+1)
  if(settings$minbusco > 1){
    regdb <- filter(regdb,BuscoID>=settings$minbusco)
    logWrite(paste('#BLOCK Reduced to',nrow(regdb),"synteny blocks based on minbusco=INT filtering."))
  }
  setminregion <- FALSE
  regdb <- regdb %>% filter(Length>=settings$minregion,HitLength>=settings$minregion)
  while(settings$maxregions > 0 & nrow(regdb) > settings$maxregions){
    settings$minregion <- min(c(regdb$Length,regdb$HitLength))+1
    regdb <- regdb %>% filter(Length>=settings$minregion,HitLength>=settings$minregion)
    setminregion <- TRUE
  }
  if(setminregion){
    logWrite(paste('#MINREG minregion=INT increased to',settings$minregion,"bp to restrict to maxregions=INT",settings$maxregions,"synteny blocks."))
  }
  logWrite(paste('#BLOCK Reduced to',nrow(regdb),"synteny blocks based on minregion=INT filtering."))
}

##### ======================== Add pre-generated regdb data ======================== #####
if(settings$regdata != ""){
  linkdb <- regTable(settings$regdata) %>%
    filter(Length>=settings$minregion,HitLength>=settings$minregion)
  logWrite(paste('#REGION Reduced to',nrow(linkdb),"synteny blocks based on minregion=INT filtering."))
  if(nrow(linkdb)){
    linkdb <- inner_join(linkdb,seqdb %>% select(Genome,SeqName,SeqLen)) %>% select(-SeqLen)
    linkdb <- inner_join(linkdb,seqdb %>% select(Genome,SeqName,SeqLen) %>% rename(HitGenome=Genome,Hit=SeqName)) %>% select(-SeqLen)
    linkdb <- unique(linkdb) #?# Why is this needed?!
    logWrite(paste('#REGION Reduced to',nrow(linkdb),"synteny blocks after removing filtered sequences."))
    if(! "BuscoID" %in% colnames(linkdb)){
      linkdb <- linkdb %>% mutate(BuscoID=1,BuscoIDList="-")
    }
    linkdb <- linkdb %>% 
      select(Genome,HitGenome,SeqName,Start,End,Strand,Hit,HitStart,HitEnd,BuscoID,Length,HitLength,BuscoIDList) %>%
      arrange(Genome,HitGenome,SeqName,Start,End)
  }
  if(nrow(regdb) > 0){
    if(nrow(linkdb) > 0){
      regdb <- bind_rows(regdb,linkdb)
      logWrite(paste('#BLOCK',nrow(regdb),"combined region+BUSCO synteny blocks."))
    }else{
      logWrite(paste('#BLOCK',nrow(regdb),"BUSCO synteny blocks: no loaded region links to add."))
    }
  }else{
    regdb <- linkdb
  }
  setminregion <- FALSE
  regdb <- regdb %>% filter(Length>=settings$minregion,HitLength>=settings$minregion)
  while(settings$maxregions > 0 & nrow(regdb) > settings$maxregions){
    settings$minregion <- min(c(regdb$Length,regdb$HitLength))+1
    regdb <- regdb %>% filter(Length>=settings$minregion,HitLength>=settings$minregion)
    setminregion <- TRUE
  }
  if(setminregion){
    logWrite(paste('#MINREG minregion=INT increased to',settings$minregion,"bp to restrict to maxregions=INT",settings$maxregions,"synteny blocks."))
    logWrite(paste('#REGION Reduced to',nrow(regdb),"synteny blocks based on minregion=INT filtering."))
  }
}
#!# Currently need some linkages. Future versions should be able to plot just chromosomes.
if(! nrow(regdb) & ! settings$dev){
  msg <- "No BUSCO links or pre-generated syntenic regions loaded."
  if(settings$rscript){
    logWrite(msg)
    quit("no",2)  
  }else{
    stop(msg)
  }
}

### ~ Optional sequence filter based on BUSCO & Link data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! settings$orphans){
  fulln <- nrow(seqdb)
  buscseq <- group_by(regdb,Genome,SeqName) %>% summarise(Genome=first(Genome),SeqName=first(SeqName),N=n())
  seqdb <- inner_join(seqdb,buscseq) %>% filter(N>0) %>% select(-N)
  seqdb <- unique(seqdb) #? Getting weird duplication issues with inner_join. Why?
  logWrite(paste("#ORPHAN Removed orphan sequences w/o synteny blocks:",fulln,"->",nrow(seqdb),"sequences."))
}

### ~ Filter genomes without any synteny blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
seqgenomes <- unique(seqdb$Genome)
badgenomes <- gendb$Genome[! gendb$Genome %in% seqgenomes]
if(length(badgenomes) > 0){
  logWrite(paste('#NOSEQ',length(badgenomes),'genomes without sequences following filters: check minregion, minbusco, minlen and orphan settings'))
  seqdb <- seqdb[seqdb$Genome %in% seqgenomes,]
  gendb <- gendb[gendb$Genome %in% seqgenomes,]
  settings$order <- settings$order[settings$order %in% seqgenomes]
}
#gendb$Genome
#settings$order



##### ======================== Update Gap Colours ======================== #####
#!# To do:
# - Gaps within synteny blocks go blue or disappear
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
if("buscopy" %in% ls()){
  rm(buscopy)
}

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
revn <- sum(seqdb$Rev)
logWrite(paste('#REVSEQ',revn,"scaffolds flipped and R suffixes added."))

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
    #logWrite(paste('#ORDER',genome,'sequences:',paste(seqorder[[genome]],collapse=", ")))
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
    #logWrite(paste('#ORDER',genome,'sequences:',paste(seqorder[[genome]],collapse=", ")))
  }
}
#i# seqsort=auto will use the orientation focus (genfocus)
if(settings$seqsort == "auto"){
  while(length(seqorder) < length(settings$order)){
    for(genome in settings$order[!settings$order %in% names(seqorder)]){
      myfocus <- genfocus[[genome]]
      if(myfocus == genome){
        seqorder[[genome]] <- seqdb[seqdb$Genome==genome,]$SeqName
      }
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
      #logWrite(paste('#ORDER',genome,'sequences:',paste(seqorder[[genome]],collapse=", ")))
    }
  }
}
#i# seqorder now has the order of the chromosomes for each genome -> check integrity
for(genome in names(seqorder)){
  seqnames <- unique(seqdb[seqdb$Genome == genome,]$SeqName)
  badseq <- seqorder[[genome]][! seqorder[[genome]] %in% seqnames]
  if(length(badseq)>0){
    #?# Why is "" sometimes in the query genome order?! - Because settings$seqorder == ""
    if(length(badseq)>1 | badseq[1] != ""){ 
      logWrite(paste("Some ordered",genome,"sequences missing. Removed as orphans?:",paste(badseq,collapse=", ")))
    }
    seqorder[[genome]] <- seqorder[[genome]][seqorder[[genome]] %in% seqnames]
  }
  misseq <- seqnames[! seqnames %in% seqorder[[genome]]]
  if(length(misseq)>0){
    logWrite(paste("Missing",genome,"sequences added to order list:",paste(misseq,collapse=", ")))
    seqorder[[genome]] <- c(seqorder[[genome]],misseq)
  }
  logWrite(paste('#ORDER',genome,'sequences:',paste(seqorder[[genome]],collapse=", ")))
}
#i# Debug check
if(settings$debug){
  logWrite(paste("#SEQS",nrow(seqdb),"sequences in total."))
  logWrite(paste("#BUSCO",nrow(busdb),"BUSCO genes in total."))
  logWrite(paste("#DUPL",nrow(dupdb),"Duplicated BUSCO genes in total."))
  logWrite(paste("#TIDK",nrow(teldb),"TIDK telomere windows in total."))
  logWrite(paste("#FT",nrow(ftdb),"features in total."))
  logWrite(paste("#GAPS",nrow(gapdb),"assembly gaps in total."))
}


##### ======================== Setup Plots ======================== #####
#i# settings$order sets the vertical order
#i# seqorder[[genome]] sets the horizontal order for each genome
#i# regdb sets the linkages

### ~ Establish max genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite("Establish max genome size, gaps and padding")
gendb <- group_by(seqdb,Genome) %>% summarise(Genome=first(Genome),GenLen=sum(SeqLen),SeqNum=n())
#i# Check for sequences in order
badord <- settings$order[! settings$order %in% gendb$Genome]
if(length(badord)>0){
  logWrite(paste("Problem with ordered genomes lacking sequences:",paste(badord,collapse=",")))
  settings$order <- settings$order[settings$order %in% gendb$Genome]
}  
badgen <- gendb$Genome[! gendb$Genome %in% settings$order]
if(length(badgen)>0){
  logWrite(paste("Problem with sequences not matching ordered genomes for",paste(badgen,collapse=",")))
  gendb <- gendb[gendb$Genome %in% settings$order,]
}

### ~ Rescaling ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Rescaling if scale=pc - needs to change all positions and chromosome sizes
pcScale <- function(gendb,dbtable,dbtxt){
  if(! nrow(dbtable)){ return(dbtable) }
  pcdb <- left_join(dbtable,gendb %>% select(Genome,GenLen))
  fields <- c()
  for(field in c("Pos","Start","End","RevPos","SeqLen")){
    if(field %in% colnames(pcdb)){
      pcdb[[field]] <- 100 * pcdb[[field]] / pcdb$GenLen
      fields <- c(fields,field)
    }
  }
  logWrite(paste0("#SCALE Rescaled ",dbtxt," fields to scale=pc: ", paste(fields,collapse=", ")))
  return(pcdb)
}

if(settings$scale == "pc"){
  # Rescale seqdb, busdb, ftdb, gapdb, regdb, teldb
  seqdb <- pcScale(gendb,seqdb,"Sequences")
  busdb <- pcScale(gendb,busdb,"BUSCO")
  dupdb <- pcScale(gendb,dupdb,"Duplicated")
  ftdb <- pcScale(gendb,ftdb,"Features")
  gapdb <- pcScale(gendb,gapdb,"Gaps")
  regdb <- pcScale(gendb,regdb,"Synteny")
  teldb <- pcScale(gendb,teldb,"Telomeres")
  # Rescale hits in regdb
  pcdb <- left_join(regdb,gendb %>% select(Genome,GenLen) %>% rename(HitGenome=Genome,HitGenLen=GenLen))
  fields <- c()
  for(field in c("HitStart","HitEnd","HitLength")){
    if(field %in% colnames(pcdb)){
      pcdb[[field]] <- 100 * pcdb[[field]] / pcdb$HitGenLen
      fields <- c(fields,field)
    }
  }
  logWrite(paste0("#SCALE Rescaled Syntenty hit fields to scale=pc: ", paste(fields,collapse=", ")))
  regdb <- pcdb
  # No need to rescale HitPosMb within bestdb: not used again.
  # Rescale gendb
  gendb <- gendb %>% rename(RawGenLen=GenLen)
  gendb$GenLen <- 100
  # Rescale ticks
  if(settings$ticks > 100){
    settings$ticks <- 1
    logWrite('#TICKS Rescaled tick marks to 1%. Set ticks < 100 to avoid rescaling.')
  }
}

### ~ Establish gaps and padding ~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#i# Extablish xgaps
xgap <- 0.01 * max(gendb$GenLen)
maxlen <- max(gendb$GenLen + ((gendb$SeqNum - 1) * xgap))
if(settings$debug){ logWrite(paste("Maximum combined chromosome and gap length =",maxlen)) }
if(settings$align == "justify"){
  gendb$xgap <- (maxlen - gendb$GenLen) / (gendb$SeqNum - 1)
  if(settings$debug){ 
    print(summary( gendb$GenLen + ((gendb$SeqNum - 1) * gendb$xgap)  ))  
  }
}else{
  gendb$xgap <- xgap
}
if(settings$debug){ logWrite("Adding xshift") }
gendb$xshift <- 0
if(settings$align == "right"){
  gendb$xshift <- maxlen - (gendb$GenLen + (gendb$SeqNum - 1) * xgap)
}
if(settings$align %in% c("centre","center")){
  gendb$xshift <- maxlen - (gendb$GenLen + (gendb$SeqNum - 1) * xgap) / 2
}
gendb$yshift <- 0 - (match(gendb$Genome,settings$order) - 1) * settings$ygap
#X# rownames(gendb) <- gendb$Genome

### ~ Establish xshift and yshift for each sequence ~~~~~~~~~~~~~~~~~~~ ###
logWrite("Establish xshift and yshift for each sequence")
seqdb$xshift <- gendb[match(seqdb$Genome,gendb$Genome),]$xshift
seqdb$yshift <- gendb[match(seqdb$Genome,gendb$Genome),]$yshift
if(settings$debug){ logWrite("Adding xshift per genome") }
if(settings$debug){ paste(gendb$Genome) }
for(genome in settings$order){
  grow <- gendb$Genome == genome
  xshift <- gendb[grow,]$xshift
  if(settings$debug){ 
    logWrite(genome) 
    logWrite(paste(seqorder[[genome]],collapse=", "))
    logWrite(paste(unique(seqdb[seqdb$Genome == genome,]$SeqName),collapse=", "))
  }
  for(seqname in seqorder[[genome]]){
    seqdb[seqdb$Genome==genome & seqdb$SeqName==seqname,]$xshift <- xshift
    xshift <- xshift + gendb[grow,]$xgap + seqLen(seqdb,genome,seqname)
  }
}

### ~ TIDK Telomeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite("TIDK Telomeres...")
#i# Combine with seqdb to get xshift and yshift
#i# Update teldb to have the correct xshift and yshift information, reversing where needed
if(nrow(teldb)){
  teldb <- right_join(teldb,seqdb)
  teldb <- teldb[! is.na(teldb$Pos),]
}
if(nrow(teldb)){
  if(settings$scale == "pc"){
    teldb$RevPos <- teldb$SeqLen - teldb$Pos #+ (1 / teldb$GenLen)
  }else{
    teldb$RevPos <- teldb$SeqLen - teldb$Pos + 1
  }
  teldb$RevStrand <- '+'
  if('+' %in% teldb$Strand){
    teldb[teldb$Strand == '+',]$RevStrand <- '-'
  }
  #i# Reverse the strand and position where needed for plotting
  if(sum(teldb$Rev)>0){
    teldb[teldb$Rev, ]$Strand <- teldb[teldb$Rev, ]$RevStrand
    teldb[teldb$Rev, ]$Pos <- teldb[teldb$Rev, ]$RevPos
  }
  teldb$yshift <- teldb$yshift + 0.25
  if('-' %in% teldb$Strand){
    teldb[teldb$Strand == '-',]$yshift <- teldb[teldb$Strand == '-',]$yshift + 0.5
  }
}
teldb <- mutate(teldb,Size=settings$ftsize)
seqdb <- mutate(seqdb,Size=settings$ftsize)

### ~ Gaps and Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
logWrite("Gaps and Features...")
#i# Combine with seqdb to get xshift and yshift
#i# Update ftdb to have the correct xshift and yshift information, reversing where needed
gapdb <- mutate(gapdb,Size=1,Feature="Gap")
if(nrow(ftdb)){
  ftdb <- mutate(ftdb,Size=settings$ftsize)
  ftdb <- left_join(ftdb,seqdb %>% select(Genome,SeqName,SeqLen,Rev,xshift,yshift))
  if(nrow(gapdb)){
    ftdb <- bind_rows(ftdb,left_join(gapdb,seqdb %>% select(Genome,SeqName,SeqLen,Rev,xshift,yshift)))
  }
  ftdb <- ftdb[! is.na(ftdb$Pos),]
}else{
  if(nrow(gapdb)){
    ftdb <- left_join(gapdb,seqdb %>% select(Genome,SeqName,SeqLen,Rev,xshift,yshift))
  }
  ftdb <- ftdb[! is.na(ftdb$Pos),]
}
#i# Add Duplicated BUSCOs
if(nrow(ftdb)){
  ftdb <- ftdb %>% select(Genome,SeqName,Pos,Strand,Feature,Col,Fill,Shape,Size,SeqLen,Rev,xshift,yshift)
}
if(nrow(dupdb)){
  dupdb <- mutate(dupdb,Size=settings$ftsize) %>% select(Genome,SeqName,Pos,Strand,Feature,Col,Fill,Shape,Size)
}
if(nrow(ftdb)){
  if(nrow(dupdb)){
    ftdb <- bind_rows(ftdb,left_join(dupdb,seqdb %>% select(Genome,SeqName,SeqLen,Rev,xshift,yshift)))
  }
  ftdb <- ftdb[! is.na(ftdb$Pos),]
}else{
  if(nrow(dupdb)){
    ftdb <- left_join(dupdb,seqdb %>% select(Genome,SeqName,SeqLen,Rev,xshift,yshift))
  }
  ftdb <- ftdb[! is.na(ftdb$Pos),]
}
#i# Correct any dodgy data
if(nrow(ftdb)){
  ftdb$Shape <- as.integer(ftdb$Shape)
  if(sum(is.na(ftdb$Shape)) > 0){
    ftdb[is.na(ftdb$Shape),]$Shape <- 22
  }
}

#i# Reverse the features where required
logWrite("Reverse Features where required...")
if(settings$debug){
  print(head(ftdb))
}
if(nrow(ftdb)){
  if(settings$scale == "pc"){
    ftdb$RevPos <- ftdb$SeqLen - ftdb$Pos #+ (1 / ftdb$GenLen)
  }else{
    ftdb$RevPos <- ftdb$SeqLen - ftdb$Pos + 1
  }
  revstrand <- data.frame(Strand=c(".","+","-"),RevStrand=c(".","-","+"))
  ftdb <- left_join(ftdb,revstrand)
  #i# Reverse the strand and position where needed for plotting
  badft <- ftdb[is.na(ftdb$Rev),]
  if(settings$debug){
    print(summary(ftdb))
    logWrite(paste("Dropped:",paste(unique(badft$SeqName),collapse=", ")))
  }
  ftdb <- ftdb[! is.na(ftdb$Rev),]
  if(sum(ftdb$Rev) > 0){
    ftdb[ftdb$Rev, ]$Strand <- ftdb[ftdb$Rev, ]$RevStrand
    ftdb[ftdb$Rev, ]$Pos <- ftdb[ftdb$Rev, ]$RevPos
  }
  ftdb$yshift <- ftdb$yshift + 0.5
  #i# Update directional triangles for -1 shapes (Duplicated BUSCOs) and -2 shapes (Repeats or Breaks etc.)
  drows <- ftdb$Strand == '.' & ftdb$Shape == -2
  if(sum(drows) > 0){
    ftdb[drows,]$Strand <- "+"
    ftdb <- bind_rows(ftdb, ftdb[drows,] %>% mutate(Strand="-"))
  }  
  if("-" %in% ftdb$Strand){
    ftdb[ftdb$Strand == '-',]$yshift <- ftdb[ftdb$Strand == '-',]$yshift - 0.4
    # BUSCO Duplicates
    drows <- ftdb$Strand == '-' & ftdb$Shape == -1
    if(sum(drows) > 0){
      ftdb[drows,]$yshift <- ftdb[drows,]$yshift - 0.1
      ftdb[drows,]$Shape <- 25
    }
    # Special Markers
    drows <- ftdb$Strand == '-' & ftdb$Shape == -2
    if(sum(drows) > 0){
      ftdb[drows,]$yshift <- ftdb[drows,]$yshift + 0.3
      ftdb[drows,]$Shape <- 24
    }
  }
  if("+" %in% ftdb$Strand){
    ftdb[ftdb$Strand == '+',]$yshift <- ftdb[ftdb$Strand == '+',]$yshift + 0.4
    # BUSCO Duplicates
    drows <- ftdb$Strand == '+' & ftdb$Shape == -1
    if(sum(drows) > 0){
      ftdb[drows,]$yshift <- ftdb[drows,]$yshift + 0.1
      ftdb[drows,]$Shape <- 24
    }
    # Special Markers
    drows <- ftdb$Strand == '+' & ftdb$Shape == -2
    if(sum(drows) > 0){
      ftdb[drows,]$yshift <- ftdb[drows,]$yshift - 0.3
      ftdb[drows,]$Shape <- 25
    }
  }
}
#i# Debug check
if(settings$debug){
  logWrite(paste("#SEQS",nrow(seqdb),"sequences in total."))
  logWrite(paste("#BUSCO",nrow(busdb),"BUSCO genes in total."))
  logWrite(paste("#TIDK",nrow(teldb),"TIDK telomere windows in total."))
  logWrite(paste("#FT",nrow(ftdb),"features+gaps+dupl in total."))
  logWrite(paste("#GAPS",nrow(gapdb),"assembly gaps in total."))
}


##### ======================== Add Gap warnings ======================== #####
fdupdb <- data.frame()
if(nrow(ftdb) > 0){
  ftdb <- arrange(ftdb,Genome,SeqName,Pos)
  #!# Find gaps between synteny blocks
  #!# Find gaps between telomere pairs
  #!@ Find telomeres between synteny blocks
  ### ~ Gaps between Duplicated BUSCOs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
  fi <- 3:nrow(ftdb)
  dupgap <- ftdb$Feature[fi] %in% dupdb$Feature & ftdb$Feature[fi-1] == "Gap" & ftdb$Feature[fi] == ftdb$Feature[fi-2]
  logWrite(paste(sum(dupgap),"potential false duplication gap(s) (mis-scaffolding) identified"))
  if(sum(dupgap)>0){
    fdupdb <- ftdb[c(FALSE,dupgap,FALSE),]
  }
}

##### ======================== Add CN Data ======================== #####
#!# Currently won't work for pc scaling. What about RevComp? May need to move sooner! (Should be OK?)
#i# Setup CN colour palette
# grey, orange, lightblue, green, yellow, blue, red, pink
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# grey, yellow, lightblue, blue, pink, red
cbPalette <- cbPalette[c(1,5,3,6,8,7)]
regCol <- function(cn){
  if(cn == 0){
    return(cbPalette[1])
  }
  if(cn < 0.33){
    return(cbPalette[2])
  }
  if(cn < 0.75){
    return(cbPalette[3])
  }
  if(cn < 1.5){
    return(cbPalette[4])
  }
  if(cn < 2.5){
    return(cbPalette[5])
  }
  return(cbPalette[6])
}
#Create a custom color scale
cnColors <- cbPalette    # [c(1,5,3,6,8,7)] (Reordered above)
names(cnColors) <- c("0n", "<1n", "1n", "2n", "4n", "6+n")
colScale <- scale_colour_manual(name = "RegCN",values = cnColors)

#i# Add optional DepthSizer or DepthKopy CN parsing.
cndb <- data.frame()
# Parse fields SeqName/Contig, Start, End, CN.
if(file.exists(settings$cndata)){
  cndb <- loadTable(settings$cndata)
  cndb <- cndb$data
  if(! "SeqName" %in% colnames(cndb) & "Contig" %in% colnames(cndb)){
    cndb <- rename(cndb,SeqName=Contig)
  }
}
# Generate Pos field
if(nrow(cndb) > 0){
  if(! "Pos" %in% colnames(cndb)){
    cndb$Pos <- (cndb$End + cndb$Start) / 2
  }
  # Set colour by CN
  cndb$ColCN <- cbPalette[6]
  cndb <- cndb %>% 
    mutate(ColCN=if_else(cndb$CN < 2.5,cbPalette[5],ColCN)) %>%
    mutate(ColCN=if_else(cndb$CN < 1.5,cbPalette[4],ColCN)) %>%
    mutate(ColCN=if_else(cndb$CN < 0.75,cbPalette[3],ColCN)) %>%
    mutate(ColCN=if_else(cndb$CN < 0.33,cbPalette[2],ColCN)) %>%
    mutate(ColCN=if_else(cndb$CN == 0,cbPalette[1],ColCN))
  # Update the ftdb table by joining on Genome, SeqName and Pos.
  # Strand and Genome are optional
  if("Strand" %in% colnames(cndb) & "Genome" %in% colnames(cndb)){
    ftdb <- left_join(ftdb,cndb %>% select(Genome, SeqName, Pos, Strand, CN, ColCN))
  }
  if("Strand" %in% colnames(cndb) & ! "Genome" %in% colnames(cndb)){
    ftdb <- left_join(ftdb,cndb %>% select(SeqName, Pos, Strand, CN, ColCN))
  }
  if(! "Strand" %in% colnames(cndb) & "Genome" %in% colnames(cndb)){
    ftdb <- left_join(ftdb,cndb %>% select(Genome, SeqName, Pos, CN, ColCN))
  }
  if(! "Strand" %in% colnames(cndb) & ! "Genome" %in% colnames(cndb)){
    ftdb <- left_join(ftdb,cndb %>% select(SeqName, Pos, CN, ColCN))
  }
  # Update colour
  logWrite(paste0(length(is.na(ftdb$ColCN))," features are lacking Copy Number"))
  ftdb[! is.na(ftdb$ColCN),]$Fill <- ftdb[! is.na(ftdb$ColCN),]$ColCN
}


##### ======================== Save data to Excel ======================== #####
### ~ Save data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(settings$writexl){
  if(endsWith(settings$basefile,"chromsyn")){
    outfile = paste0(settings$basefile,".xlsx")
  }else{
    outfile = paste0(settings$basefile,".chromsyn.xlsx")
  }
  xlftdb <- ftdb
  if(nrow(ftdb)>0){
    xlftdb <- select(ftdb,Genome,SeqName,Pos,Strand,Feature,Col,Fill,Shape)
  }
  X <- list(Genomes=gendb,Sequences=seqdb,BUSCO=busdb,Regions=regdb,Hits=hitdb,Best=bestdb,Positions=posdb,Telomeres=teldb,Features=xlftdb)
  if(nrow(fdupdb)){
    X$DupGap <- fdupdb    
  }
  write_xlsx(
    x = X,
    path = outfile
  )
  logWrite(paste("#XLSX Tabs:",paste(names(X),collapse=",")))
  logWrite(paste("#SAVE","All chromsyn data output to",outfile))
}else{
  logWrite("#XLXS No Excel output: check writexl=T/F setting and writexl package installation if expected.")
}


##### ======================== Generate Plots ======================== #####
#i# settings$order sets the vertical order
#i# seqorder[[genome]] sets the horizontal order for each genome
#i# regdb sets the linkages

### ~ Generate chromosome plot ~~~~~~~~~~~~~~~~~~~ ###
chromSynPlot <- function(gendb,seqdb,regdb,linkages=c()){
  cat("Generating plot", file = stderr())
  if(length(linkages)<1){
    linkages <- 1:(length(settings$order)-1)
  }
  scaling <- list(bp=1, kb=1e3, Mb=1e6, Gb=1e9, pc=1)
  rescale <- scaling[[settings$scale]]
  revshift <- 1
  if(settings$scale == "pc"){ revshift <- 0 }
  
  #i# First, plot the chromosomes
  cat("...genomes", file = stderr())
  pD <- seqdb %>% mutate(xmin=xshift/rescale,xmax=(xshift+SeqLen)/rescale,ymin=yshift,ymax=yshift+1) %>% 
    arrange(yshift,xshift) 
  pD$ymod[odd(1:nrow(pD))] <- 1+settings$textshift
  pD$ymod[even(1:nrow(pD))] <- -settings$textshift
  pD$texty <- pD$ymin + pD$ymod
  pC <- pD
  plt <- ggplot()

  #i# Add chromosomes and set colours
  #geom_rect(data=pD, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Genome), color="black")
  if(settings$chromfill == "Col"){
    plt <- plt + geom_rect(data=pD, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=pD$Col, color="black")
  }else{
    plt <- plt + geom_rect(data=pD, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=!!sym(settings$chromfill)), color="black")
  }
  if(settings$labels){
    plt <- plt +  annotate("text",label=pD$text,x=pD$xmin,y=pD$texty,colour="black",size=3 * settings$labelsize,hjust=0)
  }

  #i# Then, plot the linkages
  cat("...linkages", file = stderr())
  ybleed <- max(0,min(1,settings$ybleed))
  vnum <- length(linkages)
  for(v in linkages){
    #i# Setup the pair
    genomea <- settings$order[v]
    genomeb <- settings$order[v+1]
    ya <- gendb$yshift[gendb$Genome==genomea] + ybleed
    yb <- gendb$yshift[gendb$Genome==genomeb] + 1 - ybleed
    ya2 <- ya + (settings$ypad * (yb - ya)) - ybleed
    yb2 <- yb - (settings$ypad * (yb - ya)) + ybleed
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
      # if(settings$debug){
      #   logWrite(paste(plotdb$Genome[i],plotdb$SeqName[i]))
      # }
      if(seqRev(seqdb,plotdb$Genome[i],plotdb$SeqName[i])){
        tmp <- seqLen(seqdb,plotdb$Genome[i],plotdb$SeqName[i]) - xa1 + revshift
        xa1 <- seqLen(seqdb,plotdb$Genome[i],plotdb$SeqName[i]) - xa2 + revshift
        xa2 <- tmp
        fwd <- ! fwd
      }
      # if(settings$debug){
      #   logWrite(paste(plotdb$HitGenome[i],plotdb$Hit[i]))
      # }
      if(seqRev(seqdb,plotdb$HitGenome[i],plotdb$Hit[i])){
        tmp <- seqLen(seqdb,plotdb$HitGenome[i],plotdb$Hit[i]) - xb1 + revshift
        xb1 <- seqLen(seqdb,plotdb$HitGenome[i],plotdb$Hit[i]) - xb2 + revshift
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
        if(settings$ypad){
          pD <- data.frame(x=c(xa1,xa1,xa2,xa2,xb2,xb2,xb1,xb1), y=c(ya2,ya,ya,ya2,yb2,yb,yb,yb2))
        }
        plt <- plt + 
          geom_polygon(data=pD,mapping=aes(x=x, y=y),fill="steelblue", color=NA, alpha=settings$opacity) 
      }else{
        pD <- data.frame(x=c(xa1,xa2,xb1,xb2), y=c(ya,ya,yb,yb))
        if(settings$ypad){
          pD <- data.frame(x=c(xa1,xa1,xa2,xa2,xb1,xb1,xb2,xb2), y=c(ya2,ya,ya,ya2,yb2,yb,yb,yb2))
        }
        plt <- plt + 
          geom_polygon(data=pD,mapping=aes(x=x, y=y),fill="indianred", color=NA, alpha=settings$opacity) 
      }
      
    }
    cat(paste0("...(",v,") ",round(100*which(v==linkages)/vnum,1),"%"), file = stderr())
  }

  #i# Add genome labels
  cat("...labels", file = stderr())    
  pG <- group_by(pC,Genome) %>% summarise(Genome=first(Genome),ymax=min(ymax),xmin=min(xmin))
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
  #!# Add plotting position numbers to tickmarks
  # Features
  cat("...features", file = stderr())    
  #?# Could have a feature type field and map colour by type?
  if(settings$debug){ logWrite(paste(nrow(ftdb),"features & gaps to plot")) }
  if(nrow(ftdb)){
    pD <- ftdb %>% mutate(xpos=(xshift+Pos)/rescale)
    plt <- plt + geom_point(data=pD,mapping=aes(x=xpos,y=yshift),colour=pD$Col,fill=pD$Fill,shape=pD$Shape,size=pD$Size)
    if("CN" %in% colnames(ftdb) & sum(! is.na(ftdb$CN)) > 0){
      plt <- plt + colScale
    }
  }
  # Telomeres
  cat("...telomeres", file = stderr())    
  pD <- seqdb %>% mutate(xmin=xshift/rescale,xmax=(xshift+SeqLen)/rescale,y=yshift+0.5)
  pD <- pD[ (pD$Tel5 & ! pD$Rev) | (pD$Tel3 & pD$Rev), ]
  if(settings$debug){ logWrite(paste(nrow(pD),"5' telomeres to plot")) }
  if(nrow(pD)){
    plt <- plt + geom_point(data=pD,mapping=aes(x=xmin,y=y),colour="black",size=pD$Size)
  }
  pD <- seqdb %>% mutate(xmin=xshift/rescale,xmax=(xshift+SeqLen)/rescale,y=yshift+0.5)
  pD <- pD[ (pD$Tel5 & pD$Rev) | (pD$Tel3 & ! pD$Rev), ]
  if(settings$debug){ logWrite(paste(nrow(pD),"3' telomeres to plot")) }
  if(nrow(pD)){
    plt <- plt + geom_point(data=pD,mapping=aes(x=xmax,y=y),colour="black",size=pD$Size)
  }
  # TIDK internal windows
  if(settings$debug){ logWrite(paste(nrow(teldb),"3' TIDK telomeres to plot")) }
  if(nrow(teldb)){
    pD <- teldb %>% mutate(xpos=(xshift+Pos)/rescale)
    plt <- plt + geom_point(data=pD,mapping=aes(x=xpos,y=yshift),colour="blue",size=pD$Size)
  }

  #i# Add the theme
  cat("...theme", file = stderr())    
  plt <- plt + theme_bw() + 
    xlab(paste0("Position (",settings$scale,")")) + ylab("") +
    theme(#legend.position = "None", 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
  if(settings$labels){  # & ! "CN" %in% colnames(ftdb)){ #!# Need to add CN legend
    plt <- plt + theme(legend.position = "None")
  }
  cat("\n", file = stderr())    
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

plotted <- FALSE
while(! plotted & plotsplits < linknum){
  plotsplits <- (plotsplits + 1)
  #i# Try generating plots until fragmented enough
  tryCatch(
    expr = {
      cat(paste(plotsplits,"plot(s)..."), file = stderr())
      for(px in 1:plotsplits){
        plt <- NA
        if(plotsplits > 1){
          linkn <- round(linknum/plotsplits) + 1
          linkages <- 1:linkn + max(0, linkn * (px - 1) -1 )
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
        ggsave(plotfile,plot=plt,device="pdf",path=settings$plotdir,width=pdfwidth,height=pdfheight,scale=settings$pdfscale,limitsize = FALSE)
        logWrite(paste0('#GGSAVE Saved output plot to ',settings$plotdir,plotfile))

        plotfile <- paste0(plotbase,".png")
        ggsave(plotfile,plot=plt,device="png",path=settings$plotdir,width=pdfwidth,height=pdfheight,scale=settings$pdfscale,limitsize = FALSE)
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
      plotted <- TRUE
    },
    error = function(e){
      print(e)
      logWrite(paste("#ERROR",e,"=> splitting plot into",plotsplits+1))
      plt <- NA
      plotted <- FALSE
    }
  )
  
}


##### ======================== Tidy and Finish ======================== #####
if(file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

### ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
options(warn = oldwarn)
logWrite("#RCODE ChromSyn.R finished.")
#quit("no",0)

