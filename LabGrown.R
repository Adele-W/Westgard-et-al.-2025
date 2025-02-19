#### Clean data and identify laboratory-grown chambers ####

# Copyright: Adele Westgard 2025
# Citation: 


##################################################
#Program to:
  # Clean laser ablation foram data
  # Check if specimens are lab grown
  # Make record of which specimens are lab-grown and not
  # Save clean data as individual csv files
################################################

# Load packages
library("mixtools")
library("zoo")
library(tidyverse)

#Set working directory, to folder containing laser ablation data
setwd("")

#Prepare lists to categorise chambers
DidNOtGrow <- list("")
NotThrough <- list("")
CutThrough <- list("")
Na_count <- list("")
Contam <- list("")
PartLab <- list("")

#Molar mass for converting ppm to molar
MBa138 <- 137.905
MBa135 <- 134.906
MBa137 <- 136.906
MCa43 <- 42.959
MMg24 <- 24.305
MNa23 <- 22.989769


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
  
  #Read file into data frame
  Foram <- read.csv(ForamName, skip = 1)
  
#### Prepare data ####
  #Load common trace elements to be used, 
      #edit column number based on file structure
  abl_time <- Foram[,2] #ablation time
  Ca_ppm <- Foram[,8]
  Ba_ppm <- Foram[,12] #135Ba
  Mg_ppm <- Foram[,6]
  Na_ppm <- Foram[,5]
  
  #Convert to molar and calculate ratios
  Ca_mol <- list() 
  Ba_mmol <- list()
  Mg_mmol <- list()
  BaCa_M <- list()
  MgCa_M <- list()

  for (p in 1:length(Ca_ppm)){
    #Calcium to mol
    mol <- (Foram[p, 8]/MCa43)/1000
    Ca_mol <- append(Ca_mol, mol)
    
    #Barium 135 to mmol
    mmol <- (Foram[p, 12]/MBa135)
    Ba_mmol <- append(Ba_mmol, mmol)
    
    #Magnesium to mmol
    Mg_M <- (Foram[p, 6]/MMg24)
    Mg_mmol <- append(Mg_mmol, Mg_M)
    
    #135Ba/Ca mmol/mol
    Baratio <- (mmol/mol)
    BaCa_M <- append(BaCa_M, Baratio)
    
    #Mg/Ca mmol/mol
    Mgratio <- (Mg_M/mol)
    MgCa_M <- append(MgCa_M, Mgratio)
  }
  
  #Calculate rolling average for Mg/Ca mmol/mol
  MgCa_mol <- as.numeric(Mg_mmol)/as.numeric(Ca_mol) 
  MgCa_roll_mol <- rollmean(MgCa_mol, 
                            5, na.pad = TRUE, align = "center")
  
  
#### Is the specimen labgrown? ####
  #Find mean 135Ba/Ca of initial ~10 seconds
  Ten_sec <- which(abl_time >= 10, abl_time)[1] #find index of first 10 sec
  short_BaCa <- BaCa_M[5:Ten_sec]          #extract first ~10 sec of 135Ba/Ca
  avgBaCa_10sec <- mean(as.numeric(short_BaCa)) #Find mean of first 10 sec
  
  #If first 10 sec mean 135Ba/Ca < 0.008 mmol/mol: not laboratory-grown. 
  if  (avgBaCa_10sec < 0.008) {
    #make record AND skip to next profile
    DidNOtGrow <- append(DidNOtGrow, ForamName)
    next #next chamber
    
  } else { #The chamber is at least partially laboratory-grown
    
    #create new data frame to add clean data to
    Foram_clean <- data.frame()
    
#### Remove initial contamination: ####
    #Remove initial 3 seconds 
    Foram_trim <- filter(Foram, Foram[ , 2]>3)
    
    #Then remove based on excessively high Na/Ca in the first 10 sec
    short_NaCa <- rollmean((Na_ppm/Ca_ppm)[1:Ten_sec], 
                           5, na.pad = TRUE, align = "center")
    short_NaCa <- na.omit(short_NaCa)
    
    #If any Na/Ca in first 10 sec are higher than 0.01, there is contamination
    if (any(short_NaCa > 0.01, na.rm = TRUE)) {
      #removes the high values initially
      if (any(short_NaCa < 0.01)) {
        index_NaCa <- (which(short_NaCa < 0.01, )[1])+3
        Foram_trim <- Foram_trim[index_NaCa:(length(Na_ppm/Ca_ppm)), ]
      } else {
        #All NaCa are >0.01, the first 10 sec of this chamber is contaminated 
        #Record as contaminated
        Contam <- append(Contam, ForamName)
      }
    }
    
#### Which part is laboratory-grown? ####
    
    #Find rolling mean of 135Ba/Ca mmol/mol
    BaCa_M_num <- as.numeric(BaCa_M) 
    #rolling mean with 5 values in each bracket
    ma_BaCa <- rollmean(BaCa_M_num, 5, na.pad = TRUE, align = "center")
    ma_BaCa_NA <- na.omit(ma_BaCa)
    
    #rolling mean with 8 values in each bracket  
    means <- rollmean(BaCa_M_num, 8, na.pad = TRUE, align = "center")
    means <- na.omit(means)  
    
    #Labgrown transition
    rate_change <- diff(ma_BaCa)
    index_max_Rate <- (which.max(rate_change) + 3)
    
    #Through using increase in NaCa, i.e., analysing mounting tape
    NaCa <- Foram_trim[,5]/Foram_trim[,8] #Na/Ca ppm of cleaned file
    ma_NaCa <- rollmean(NaCa, 5, na.pad = TRUE, align = "center")
    ma_NaCa_NA <- na.omit(ma_NaCa)
    inc_NaCa <- diff(ma_NaCa_NA)
    
    #All lab-grown and not analysing something else than intended chamber?
    if (any(means < 0.008, na.rm = TRUE) |        
        any(ma_NaCa_NA > 0.01, na.rm = TRUE) | 
        any(MgCa_roll_mol > 8, na.rm = TRUE)) {
      
      #Laser cut through the shell wall:
      if (any(ma_NaCa_NA > 0.01, na.rm = TRUE)) {
        index_inc_NaCa <- which(ma_NaCa_NA > 0.01)[1] + 1
        Foram_trim <- Foram_trim[1:index_inc_NaCa, ]
        
        #Keep record
        Na_count <- append(Na_count, ForamName)
        CutThrough <- append(CutThrough, ForamName)
      }

      #rolling mean of  Mg/Ca after 10 sec 
      rollMg_10 <- na.omit(MgCa_roll_mol[Ten_sec:length(MgCa_roll_mol)])
      
      #Laser cut through the shell wall:
      if (any(rollMg_10 > 8, na.rm = TRUE)) {
        roll_10 <- which(rollMg_10 > 8, rollMg_10)[1]
        temp_ind <- which(MgCa_roll_mol == rollMg_10[roll_10], MgCa_roll_mol)
        
        if (length(temp_ind) > 1) {
          index_inc_MgCa <- which(temp_ind > Ten_sec)
        } else {
          index_inc_MgCa <- temp_ind}
        
        if (index_inc_MgCa <= length(Foram_trim[,2])) {
          Foram_trim <- Foram_trim[1:index_inc_MgCa, ] #remove cut-through
          CutThrough <- append(CutThrough, ForamName) } #keep record
      }
      
      #Check if remaining data after cut-thorugh is all lab-grown
      len_foram <- length(Foram_trim[,2])
      BaCa_M_trim <- BaCa_M_num[1:len_foram]
      means_trim <- rollmean(BaCa_M_trim, 8, na.pad = TRUE, align = "center")
      
      #Remove data from not-laboratory-grown component
      if (any(means_trim < 0.008, na.rm = TRUE)) {
        index_labgrown <- (which(means_trim < 0.008, means_trim)[1] + 4)
        Foram_trim <- Foram_trim[1:index_labgrown, ] 
        
        #Keep record, if true, this is partially laboratory-grown
        if (!(ForamName %in% CutThrough) & !(ForamName %in% Na_count)) {
          PartLab <- append(PartLab, ForamName)}
      }
      
      #If the remaining data has less than 5 data points, do not keep.
      if (length(Foram_trim[,2]) < 5) { #Too few datapoints - skip
        Contam <- append(Contam, ForamName)
        next
        
      #More than 5 data points: keep
      } else { #Keep data
        Foram_clean <- Foram_trim }
     
    #Laser did not cut through the shell wall AND it is lab-grown 
    } else {
      #Leave data as is, only initial contamination removed
      Foram_clean <- Foram_trim 
      NotThrough <- append(NotThrough, ForamName) #keep record
    }
    
#### create csv file with clean data  ####
    # NB! Need to match with chamber ID for IDCrust.
    csv_file <- paste0("Foram_clean-", i, ".csv")
    write.csv(Foram_clean, file = csv_file, row.names = FALSE)
    
  } #End laboratory- grown statement
} #end loop

#### Keep Reocord ####
#Add record-lists to summary of all chambers from folder. 
Foram_Sum_list <- list(DidNOtGrow, 
                       NotThrough, CutThrough, Na_count, Contam, PartLab)

Foram_Sum <- data.frame(lapply(Foram_Sum_list, function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(Foram_Sum_list))
  return(x)
}))

#Edit to include only relevant categories. 
colnames(Foram_Sum) <- 
  c("Not labgrown", "Not through", 
    "Cut through", "Na Cut-through", "Contaminated", "Part labgrown")

# Create csv file for records. edit file name
csv_summary <- "File_Name.csv"
write.csv(Foram_Sum, file = csv_summary, row.names = FALSE)

#### END ####
