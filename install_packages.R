##################################################################################
### INSTALL_PACKAGES: Installs R packages needed for SLiMSuite R scripts ~~~~~ ###
### VERSION: 0.2.0                                                       ~~~~~ ###
### LAST EDIT: 09/05/23                                                  ~~~~~ ###
### AUTHORS: Richard Edwards 2023                                        ~~~~~ ###
### CONTACT: rich.edwards@uwa.edu.au / GitHub @slimsuite                 ~~~~~ ###
##################################################################################

# This script is for installing the R packages required by SLiMSuite R scripts.

####################################### ::: HISTORY ::: ############################################
# v0.1.0 : Initial version for ChromSyn and DepthSizer.
# v0.2.0 : Updated code to run either standalone or as called from another script with source().
myversion <- "v0.2.0"

####################################### ::: FUNCTIONS ::: ##########################################
### ~ logWrite function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
if(! "logWrite" %in% ls()){
  logWrite <- function(logstr){
    writeLines(paste0("[",date(),"] ",logstr),con=stdout())
  }
}
logWrite(paste("#RCODE install_packages.R:",myversion))
logWrite(paste("#PATH Running from:",getwd()))

####################################### ::: INSTALL ::: ############################################
#i# Check and install packages
if(! "packages" %in% ls()){
  packages <- c("tidyverse","RColorBrewer","gtools","ggstatsplot","ggrepel","ggridges","writexl","doParallel","foreach","DT")
}
if("settings" %in% ls() && "packages" %in% names(settings)){
  packages <- settings$packages
}
for(rpack in packages){
  if(! rpack %in% installed.packages()[,"Package"]){
    install.packages(rpack, repos = "https://cloud.r-project.org/")
  }
}
#i# Report installation status
for(rpack in packages){
  logWrite(paste('#RPACK Package',rpack,'installed:',rpack %in% installed.packages()[,"Package"]))
}

logWrite("#RCODE install_packages.R finished.")


