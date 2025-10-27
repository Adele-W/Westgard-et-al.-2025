#### Clean data and identify laboratory-grown chambers ####

# Copyright: Adele Westgaard 2025
# Citation: 


##################################################
#Program to:
#Read trace element data csv (export from iolite) in ppm
#Remove initial contamination 
#Find laser cut-through point 
# Save clean data as individual csv files
################################################

# Load packages
library("mixtools")
library("zoo")
library(tidyverse)

#Set working directory, to folder containing laser ablation data
setwd("")

#Prepare lists to categorise chambers
NotThrough <- list("")
CutThrough <- list("")
Na_count <- list("")
Contam <- list("")
Cleaned <- list()

#Molar mass for converting ppm to molar
MBa138 <- 137.905
MBa135 <- 134.906
MBa137 <- 136.906
MCa43 <- 42.959
MMg24 <- 24.305
MNa23 <- 22.989769
MAl27 <- 26.9815
MMn55 <- 54.9380 


#### Access data file ####
#Loop: access laser ablation data for each chamber (individual csv files).
#Add total number of laser profiles below (numeric): 
chambers = 
  for (i in 1: chambers) {
    
  #Load new laser profile, edit to fit the name of files/Laser ID. 
  ForamName <- paste0("Foram-", i, ".csv")
    
  #Skips file that doesn't exist
  path = file.path(Fname)
  if(!file.exists(path)) {
     next
  }
    
  #read file into data frame
  Foram <- read.csv(ForamName, skip = 1)
    
  #Remove any chamber with significant contamination. 
  #Adjust cut-off value to user needs
  #Other element ratios can be added as required. e.g., Mn/Ca or Al/Ca
  Na <- mean(Foram[,5])/MNa23
  Mg <- mean(Foram[,6])/MMg24
  Mn <- mean(Foram[,10])/MMn55
  Al <- mean(Foram[,7])/MAl27
  Ca <- (mean(Foram[,8])/MCa43)/1000
  NaCa_avg_mol <- Na/Ca 
  MgCa_avg_mol <- Mg/Ca
  MnCa_avg_mol <- Mn/Ca
  AlCa_avg_mol <- Al/Ca
    
  if (NaCa_avg_mol > 20 | MgCa_avg_mol > 10) { #contaminated
    next }
  
  #Remove #'s below to use Mn/Ca or Al/Ca 
  # if (NaCa_avg_mol > X | MgCa_avg_mol > X |
  # MnCa_avg_mol > X | MnCa_avg_mol > X) { #replace X with appropriate values
    #next }
    
    
  #### Prepare data ####
  #Load common trace elements to be used, 
  #edit column number based on file structure
  abl_time <- Foram[,2] #ablation time
  Ca_ppm <- Foram[,8]
  Mg_ppm <- Foram[,6]
  Na_ppm <- Foram[,5]
  Al_ppm <- Foram[,7]
  Mn_ppm <- Foram[,10]
    
  #Convert to molar and calculate ratios
  Ca_mol <- list() 
  Mg_mmol <- list()
  Al_mmol <- list()
  Mn_mmol <- list()
    
  BaCa_M <- list()
  MgCa_M <- list()
  AlCa_M <- list()
  MnCa_M <- list()
    
  for (p in 1:length(Ca_ppm)){
    ##elements
    #Calcium to mol
    mol <- (Foram[p, 8]/MCa43)/1000
    Ca_mol <- append(Ca_mol, mol)
    
    #Magnesium to mmol
    Mg_M <- (Foram[p, 6]/MMg24)
    Mg_mmol <- append(Mg_mmol, Mg_M)
      
    #Aluminium to mmol
    Al_M <- (Foram[p, 7]/MAl27)
    Al_mmol <- append(Al_mmol, Al_M)
      
    #Manganese to mmol
    Mn_M <- (Foram[p, 10]/MMn55)
    Mn_mmol <- append(Mn_mmol, Mn_M)
      
      
    ## Ratios
    #Mg/Ca mmol/mol
    Mgratio <- (Mg_M/mol)
    MgCa_M <- append(MgCa_M, Mgratio)
      
    #Al/Ca mmol/mol
    Alratio <- (Al_M/mol)
    AlCa_M <- append(AlCa_M, Alratio)
      
    #Mn/Ca mmol/mol
    Mnratio <- (Mn_M/mol)
    MnCa_M <- append(MnCa_M, Mnratio)
  }
    
    
  #Calculate rolling average for Mg/Ca mmol/mol
  MgCa_mol <- as.numeric(Mg_mmol)/as.numeric(Ca_mol) 
  MgCa_roll_mol <- rollmean(MgCa_mol, 
                            5, na.pad = TRUE, align = "center")
  
    
  #create new data frame to add clean data to
  Foram_clean <- data.frame()
      
  #Remove initial contamination, 
  #REMOVING first 3 sec.
  Foram_trim <- filter(Foram, Foram[ , 2]>3)
    
  #then based on high element ratios in the first 10 sec
  Ten_sec <- which(Foram[,2] >= 10)[1]
  short_NaCa <- rollmean((Foram[,5]/Foram[,8])[1:Ten_sec], 
                         5, na.pad = TRUE, align = "center")
  short_NaCa <- na.omit(short_NaCa)
  if (any(short_NaCa > 0.01, na.rm = TRUE)) {
    if (any(short_NaCa < 0.01)) {
      index_NaCa <- (which(short_NaCa < 0.01, )[1])+3
      Foram_trim <- Foram_trim[index_NaCa:(length(Foram[,5]/Foram[,8])), ]
    } else {
      #All NaCa are >0.01, the first 10 sec of this foram is contaminated 
      #      and shouldn't be used
      #The data will be cleaned nevertheless
      Contam <- append(Contam, ForamName) #document contaminated shell
    }
  }
      
   
  #Check if laser cut through and look for contamination 
  #Through using increase in NaCa, i.e., analysing mounting tape/background
    NaCa <- Foram_trim[,5]/Foram_trim[,8] #Na/Ca ppm of cleaned file
    ma_NaCa <- rollmean(NaCa, 5, na.pad = TRUE, align = "center")
    ma_NaCa_NA <- na.omit(ma_NaCa)
    inc_NaCa <- diff(ma_NaCa_NA)
      
    #this will cut off a little too late or not at all
    #Here Mn/Ca and Al/Ca can be added based on the limit (X/Y) the user sets.
    #Remove # from lines 179-182 below to use other elements & adapt as needed. 
    #Remove lines you 184-185 if you use lines 179-182
    
    #if (any(ma_NaCa_NA > 0.01, na.rm = TRUE) |
      #any(MgCa_roll_mol > 8, na.rm = TRUE) |
      #any(AlCa_M > X, na.rm = TRUE) |
      #any(MnCa_M > Y, na.rm = TRUE)) {
      
    if (any(ma_NaCa_NA > 0.01, na.rm = TRUE) | 
        any(MgCa_roll_mol > 8, na.rm = TRUE)) {
      
        CutThrough <- append(CutThrough, ForamName) #document
      
      if (any(ma_NaCa_NA > 0.01, na.rm = TRUE)) {
        index_inc_NaCa <- which(ma_NaCa_NA > 0.01)[1] + 1
        Foram_trim <- Foram_trim[1:index_inc_NaCa, ]
        
        Na_count <- append(Na_count, ForamName) #document
        } 
        
      sec_15 <- which(abl_time > 15, abl_time)[1]
      rollMg_15 <- na.omit(MgCa_roll_mol[sec_15:length(MgCa_roll_mol)])
      
      if (any(rollMg_15 > 8, na.rm = TRUE)) {
        index_inc_MgCa <- which(rollMg_15 > 8, rollMg_15)[1] 
        Foram_trim <- Foram_trim[1:index_inc_MgCa, ]
        }
      
      #Remove "#" below if you want to remove contamination based on Al or Mn
      #Adjust cut-off value (X/Y) as needed
      #if (any(AlCa_M > X, na.rm = TRUE)) { 
        #index_inc_AlCa <- which(AlCa_M > X, AlCa_M)[1] 
        #Foram_trim <- Foram_trim[1:index_inc_AlCa, ]
      #}
      
      
      #if (any(MnCa_M > Y, na.rm = TRUE)) {
        #index_inc_MnCa <- which(MnCa_M > Y, MnCa_M)[1] 
        #Foram_trim <- Foram_trim[1:index_inc_MnCa, ]
      #}
      
      Foram_clean <- Foram_trim
        
    } else {
       Foram_clean <- Foram_trim #Leave data as is.
       NotThrough <- append(NotThrough, ForamName) #document
       }
      
     Cleaned <- append(Cleaned, ForamName) #document
      
     #### create csv file with clean data  ####
     # NB! Need to match with chamber ID for IDCrust.
     csv_file <- paste0("Foram_clean-", i, ".csv")
     write.csv(Foram_clean, file = csv_file, row.names = FALSE)
      
} #end loop

#### Keep Record ####
#Add record-lists to summary of all chambers from folder. 
Foram_Sum_list <- list(NotThrough, CutThrough, Na_count, Cleaned)

Foram_Sum <- data.frame(lapply(Foram_Sum_list, function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(Foram_Sum_list))
  return(x)
}))

#Edit to include only relevant categories. 
colnames(Foram_Sum) <- 
  c("Not through", "Cut through", "Na Cut-through", "Cleaned")

# Create csv file for records. edit file name
csv_summary <- "File_Name.csv"
write.csv(Foram_Sum, file = csv_summary, row.names = FALSE)

#### END ####
