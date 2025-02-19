#### Identify and distinguish crust and lamellar calcite ####

# Copyright: Adele Westg??rd 2025
# Citation: 

#N.B.! This overview file is based on laboratory-grown N.pachyderma.
  # see script LabGrown_Fossil for fossil shells 

##################################################
# Program to: 
  # Read pre-cleaned data
  # Distinguish crust and lamellar calcite
  # Separate the two components
  # Convert trace element data from ppm to mmol/mol
  # Determine element/Ca for each chamber: crust, lamellar calcite, shell wall
  # Report mean values in new data frame, categorised by chamber component
##################################################

library(tidyverse)
library("readxl") 
library(strucchange)

#Molar masses
MCa43 <- 42.959
MMg24 <- 24.305
MLi <- 6.941
MB <- 10.811
MNa <- 22.989769
MAl <- 26.982
MMn <- 54.938044
MSr <- 87.62
MBa135 <- 134.906
MBa137 <- 136.906
MBa138 <- 137.905


#### Access data ####

#Set working directory, to folder containing laser ablation data
setwd("")

#N.B! An overview file was created in excel which lists which analysis session 
#   data was obtained from, filenames of cleaned data, and the water conditions
#   laboratory-grown foraminifera grew in with individual chamber ID's.

#Read overview file: 
# This is used to access data and creates the base for saving new data. 
# Foram ID, chamber ID (matches ID from LabGrown), growth conditions
overview <- read_excel("All foraminifera info sheet.xlsx") #edit name

#Find total length of data sheet
len <- as.numeric(count(overview[1])) 
#Convert relevant columns of data into numeric values
col_to_num <- 11:40
overview[, col_to_num] <- lapply(overview[, col_to_num], as.numeric)


#Loop: uses overview as base for accessing data files.
# Each analytical session has its own folder of data.
for (i in 2:len) { #from row below titles to end of sheet
  
  #Analytical session 1: 
  if (overview[i, 1] == "Session 1" ){
    
    #Set working directory
    setwd("Session 1")
    Fname <- paste0(overview[i, 2], ".csv")
    
    #Skip files that doesn't exits
    path = file.path(Fname)
    if(!file.exists(path)) {
      next
    }
    
    #Skip contaminated chambers, obtained from script 1 XYZ 
    if (Fname == 'Foram_clean-X.csv') {next}
    
    #Read csv of cleaned chamber laser profile
    Foram <- read.csv(Fname)
    
    #Access the not cleaned version of the same chamber
    not_clean <- paste0(overview[i, 2], ".csv")
    not_clean <- gsub("_clean", "", not_clean)
    
    #Read file containing record of cleaning from script 1 XYZ
    slide <- read.csv("Cleaning_summary.csv")
  } # END session 1 if statement
  
  #Analytical session 2: !!! Make more of these as needed.  
  if (overview[i, 1] == "Session2" ){
    
    #Set working directory
    setwd("Session 2")
    Fname <- paste0(overview[i, 2], ".csv")
    
    #Skip files that doesn't exits
    path = file.path(Fname)
    if(!file.exists(path)) {
      next
    }
    
    #Skip contaminated chambers, obtained from script 1 XYZ 
    if (Fname == 'Foram_clean-X.csv') {next}
    
    #Read csv of cleaned chamber laser profile
    Foram <- read.csv(Fname)
    
    #Access the not cleaned version of the same chamber
    not_clean <- paste0(overview[i, 2], ".csv")
    not_clean <- gsub("_clean", "", not_clean)
    
    #Read file containing record of cleaning from script 1 XYZ
    slide <- read.csv("Cleaning_summary.csv")
  } # END session 2 if statement
  
  
#### Prepare data ####
  #Prepare elements, edit column numbers as needed. 
  abl_time <- na.omit(Foram[,2])
  Mg_ppm <- Foram[,6]
  Ca43_ppm <- Foram[,8]
  
  Li_ppm <- Foram[,3]
  B_ppm <- Foram[,4]
  Na_ppm <- Foram[,5]
  Al_ppm <- Foram[,7]
  Mn_ppm <- Foram[,10]
  Sr_ppm <- Foram[,11]
  Ba135_ppm <- Foram[,12]
  Ba137_ppm <- Foram[,13]
  Ba138_ppm <- Foram[,14]  
  
  #Rolling mean Mg/Ca
  MgCa <- na.omit(Mg_ppm / Ca43_ppm)
  means <- rollmean(MgCa, 5, na.pad = TRUE, align = "center")
  
  #Index first 5 seconds of profile, if longer than 5 seconds. 
  if (abl_time[1] < 5.5) {
    sec_5 <- which(abl_time < 5.5, abl_time)
    sec_5 <- sec_5[length(sec_5)]
  } else {
    sec_5 <- 1
  } 
  
  #remove initial 5 seconds  
  shell_Mg <- Mg_ppm[sec_5:length(Mg_ppm)]
  shell_Ca <- Ca43_ppm[sec_5:length(Ca43_ppm)]
  shell_Ba <- Ba135_ppm[sec_5:length(Ba135_ppm)]
  
  #Convert to molarity
  Ca_mol <- c()
  for(c in 1:(length(shell_Ca))){
    Ca <- (shell_Ca[c]/MCa43)/1000
    Ca_mol <- c(Ca_mol, Ca)
  }
  
  Ba_mol <- c()
  for(b in 1:(length(shell_Ba))){
    Ba <- shell_Ba[b]/MBa135
    Ba_mol <- c(Ba_mol, Ba)
  }
  #Calculate 135Ba/Ca and its rolling mean
  BaCa_mol <- na.omit(Ba_mol)/na.omit(Ca_mol)
  BaCa_mol_roll <- na.omit(rollmean(BaCa_mol, 5, na.pad = TRUE, align = "center"))
  
  
  #Confirm that all not labgrown data is removed
  if (any(BaCa_mol_roll < 0.005)){
    not_labg <- which(BaCa_mol_roll < 0.005) 
    not_labg <- not_labg[1]+2 
    
    #Remove any not labgrown data if profile is longer than 5 seconds
    if (abl_time[not_labg] > abl_time[sec_5]) { 
      
      abl_time <- abl_time[sec_5:not_labg]
      shell_MgCa <- MgCa[sec_5:not_labg]  #call shell for now... = total
      
      shell_Mg <- Mg_ppm[sec_5:not_labg]
      shell_Ca <- Ca43_ppm[sec_5:not_labg]
      shell_Na <- Na_ppm[sec_5:not_labg]
      
      avg_Mg_ppm <- mean(na.omit(shell_Mg))
      avg_Ca_ppm <- mean(na.omit(shell_Ca))
      
      avg_Li_ppm <- mean(na.omit(Li_ppm[sec_5:not_labg]))
      avg_B_ppm <- mean(na.omit(B_ppm[sec_5:not_labg]))
      avg_Na_ppm <- mean(na.omit(Na_ppm[sec_5:not_labg]))
      avg_Al_ppm <- mean(na.omit(Al_ppm[sec_5:not_labg]))
      avg_Mn_ppm <- mean(na.omit(Mn_ppm[sec_5:not_labg]))
      avg_Sr_ppm <- mean(na.omit(Sr_ppm[sec_5:not_labg]))
      avg_Ba135_ppm <- mean(na.omit(Ba135_ppm[sec_5:not_labg]))
      avg_Ba137_ppm <- mean(na.omit(Ba137_ppm[sec_5:not_labg]))
      avg_Ba138_ppm <- mean(na.omit(Ba138_ppm[sec_5:not_labg]))
    } else {
      #profile < 5 sec and non-labg part is within this: do not include.
      next   
    }
    
  } else {  
    #All is laboratory-grown, remove initial 5 seconds of all (contamination)
    
    abl_time <- abl_time[sec_5:length(abl_time)]
    shell_MgCa <- MgCa[sec_5:length(MgCa)]  #call shell for now... = total
    
    shell_Mg <- Mg_ppm[sec_5:length(Mg_ppm)]
    shell_Ca <- Ca43_ppm[sec_5:length(Ca43_ppm)]
    shell_Na <- Na_ppm[sec_5:length(Na_ppm)]
    
    avg_Mg_ppm <- mean(na.omit(shell_Mg))
    avg_Ca_ppm <- mean(na.omit(shell_Ca))
    
    avg_Li_ppm <- mean(na.omit(Li_ppm[sec_5:length(Li_ppm)]))
    avg_B_ppm <- mean(na.omit(B_ppm[sec_5:length(B_ppm)]))
    avg_Na_ppm <- mean(na.omit(Na_ppm[sec_5:length(Na_ppm)]))
    avg_Al_ppm <- mean(na.omit(Al_ppm[sec_5:length(Al_ppm)]))
    avg_Mn_ppm <- mean(na.omit(Mn_ppm[sec_5:length(Mn_ppm)]))
    avg_Sr_ppm <- mean(na.omit(Sr_ppm[sec_5:length(Sr_ppm)]))
    avg_Ba135_ppm <- mean(na.omit(Ba135_ppm[sec_5:length(Ba135_ppm)]))
    avg_Ba137_ppm <- mean(na.omit(Ba137_ppm[sec_5:length(Ba137_ppm)]))
    avg_Ba138_ppm <- mean(na.omit(Ba138_ppm[sec_5:length(Ba138_ppm)]))
  }
  
 #If profile has less than 20 data points (~5 seconds of data)
  if (length(abl_time) < 20) {next}
  
  #Convert to mol and calculate ratios after shortening data
  Ca_mol <- c()
  for(c in 1:(length(shell_Ca))){
    Ca <- (shell_Ca[c]/MCa43)/1000
    Ca_mol <- c(Ca_mol, Ca)
  }
  
  Mg_mol <- c()
  for(M in 1:(length(shell_Mg))){
    Mg <- shell_Mg[M]/MMg24
    Mg_mol <- c(Mg_mol, Mg)
  }
  
  Na_mol <- c()
  for(N in 1:(length(na.omit(shell_Na)))){
    Na <- shell_Na[N]/MNa
    Na_mol <- c(Na_mol, Na)
  }
  
  NaCa_mol <- na.omit(Na_mol)/na.omit(Ca_mol)
  MgCa_mol <- na.omit(Mg_mol)/na.omit(Ca_mol)
  
  avg_Ca43_mol <- (avg_Ca_ppm/MCa43)/1000
  avg_Mg_mol <- (avg_Mg_ppm/MMg24)
  
  Li_mmol <- (avg_Li_ppm/MLi)
  B_mmol <- (avg_B_ppm/MB)
  Na_mmol <- (avg_Na_ppm/MNa)
  Al_mmol <- (avg_Al_ppm/MAl)
  Mn_mmol <- (avg_Mn_ppm/MMn)
  Sr_mmol <- (avg_Sr_ppm/MSr)
  Ba135_mmol <- (avg_Ba135_ppm/MBa135)
  Ba137_mmol <- (avg_Ba137_ppm/MBa137)
  Ba138_mmol <- (avg_Ba138_ppm/MBa138)
  
  totMgCa_mol <- avg_Mg_mol/avg_Ca43_mol
  NaCa <- Na_mmol/avg_Ca43_mol
  BCa <- B_mmol/avg_Ca43_mol
  SrCa <- Sr_mmol/avg_Ca43_mol
  Ba135Ca <- Ba135_mmol/avg_Ca43_mol
  Ba137Ca <- Ba137_mmol/avg_Ca43_mol
  Ba138Ca <- Ba138_mmol/avg_Ca43_mol
  MnCa <- Mn_mmol/avg_Ca43_mol
  AlCa <- Al_mmol/avg_Ca43_mol
  LiCa <- Li_mmol/avg_Ca43_mol
  
  #If mean Na/Ca > 20 mmol/ mol = contaminated
  if (NaCa > 20) {
    next 
  }
  
  
#### Distinguishing crust and lamellar calcite ####
  #find length of ablation in sec
  abl_max <- max(abl_time)
  
  # If profile is longer than 15 sec, likely has crust
  if (abl_max > 15){
    
    #If laser did not cut through the shell, it likely has crust only. 
    if (between(abl_max, 39, 41) | between(abl_max, 59, 61) ) { #Not through
      #"decision" = tag to run flow of script
      decision <- "crust_only" 
    
    #Maybe have both crust and lamellar calcite
    } else {
      decision <- "both" 
      
      # Decide if crusted using "breakpoints"
      breaks <- breakpoints(MgCa_mol ~ 1, h = 10)
      
      #Use line below to find breakpoints based on Na/Ca
      #breaks <- breakpoints(NaCa_mol ~ 1, h = 10)
      
      # Are there breakpoints? 
      if (is.na(breaks$breakpoints[1]) == FALSE) {
        
        # Find index for 15 second of ablation time 
        sec_15 <- which(abl_time >= 15, abl_time)
        sec_15 <- sec_15[1]
        
        #Are there breakpoints after 15 seconds of ablation? 
        if (breaks$breakpoints[length(breaks$breakpoints)] > sec_15) {
          
          # Assign crust start to breakpoint after 15 sec. 
          break1 <- which(breaks$breakpoints > sec_15, breaks$breakpoints)
          
          # if >=2 breakpoints after 15 seconds, use the second.
          if (length(break1) >= 2){
            break1 <- as.numeric(breaks$breakpoints[break1[2]])
          } else {
            break1 <- as.numeric(breaks$breakpoints[break1[1]]) }
          
          #Assign index for crust/lamellar transition location. 
          crust_loc <- break1 
          
        } else  { #No breakpoint after 20 seconds.
          decision <- "crust_only" }
       
      #No breakpints were idenitfied, longer than 15 seconds: crust only
      } else { 
        decision <- "crust_only" } 
      
      # shorter than 15 sec: lamellar calcite only. 
    }} else { 
      decision <- "lam_only" }
  
  
#### use "decision" to assign data to each shell component in dataframe ####
  
  #Crust only
  if (decision == "crust_only") {
    crust_MgCa_mol <- totMgCa_mol
    crust_NaCa <- NaCa
    crust_BCa <- BCa
    crust_SrCa <- SrCa
    crust_Ba135Ca <- Ba135Ca
    crust_Ba137Ca <- Ba137Ca
    crust_Ba138Ca <- Ba138Ca
    crust_MnCa <- MnCa
    crust_AlCa <- AlCa
    crust_LiCa <- LiCa
    
    overview[i, 12] <- crust_MgCa_mol #crust
    overview[i, 15] <- crust_NaCa
    overview[i, 18] <- crust_BCa
    overview[i, 21] <- crust_SrCa
    overview[i, 24] <- crust_Ba135Ca
    overview[i, 27] <- crust_Ba137Ca
    overview[i, 30] <- crust_Ba138Ca
    overview[i, 33] <- crust_MnCa
    overview[i, 36] <- crust_AlCa
    overview[i, 39] <- crust_LiCa
    
  #Lamellar calcite only   
  } else if (decision == "lam_only") { 
    lam_MgCa_mol <- totMgCa_mol
    lam_NaCa <- NaCa
    lam_BCa <- BCa
    lam_SrCa <- SrCa
    lam_Ba135Ca <- Ba135Ca
    lam_Ba137Ca <- Ba137Ca
    lam_Ba138Ca <- Ba138Ca
    lam_MnCa <- MnCa
    lam_AlCa <- AlCa
    lam_LiCa <- LiCa
    
    overview[i, 13] <- lam_MgCa_mol #lamellar
    overview[i, 16] <- lam_NaCa
    overview[i, 19] <- lam_BCa
    overview[i, 22] <- lam_SrCa
    overview[i, 25] <- lam_Ba135Ca
    overview[i, 28] <- lam_Ba137Ca
    overview[i, 31] <- lam_Ba138Ca
    overview[i, 34] <- lam_MnCa
    overview[i, 37] <- lam_AlCa
    overview[i, 40] <- lam_LiCa 
    
  #Has both crust and lamellar calcite  
  } else { 
    
    #Crust component
    crust <- MgCa[sec_5:(crust_loc)]
    crust_Mg <- Mg_ppm[sec_5:(crust_loc)]
    crust_Ca <- Ca43_ppm[sec_5:(crust_loc)]
    
    crust_Mg_ppm <- mean(na.omit(crust_Mg))
    crust_Ca_ppm <- mean(na.omit( crust_Ca))
    
    #convert to mmol/mol
    crust_Ca43_mol <- (crust_Ca_ppm/MCa43)/1000
    crust_Mg_mol <- (crust_Mg_ppm/MMg24)
    
    crust_MgCa_mol <- crust_Mg_mol/crust_Ca43_mol
    
    crust_Li_ppm <- mean(na.omit(Li_ppm[sec_5:crust_loc]))
    crust_B_ppm <- mean(na.omit(B_ppm[sec_5:crust_loc]))
    crust_Na_ppm <- mean(na.omit(Na_ppm[sec_5:crust_loc]))
    crust_Al_ppm <- mean(na.omit(Al_ppm[sec_5:crust_loc]))
    crust_Mn_ppm <- mean(na.omit(Mn_ppm[sec_5:crust_loc]))
    crust_Sr_ppm <- mean(na.omit(Sr_ppm[sec_5:crust_loc]))
    crust_Ba135_ppm <- mean(na.omit(Ba135_ppm[sec_5:crust_loc]))
    crust_Ba137_ppm <- mean(na.omit(Ba137_ppm[sec_5:crust_loc]))
    crust_Ba138_ppm <- mean(na.omit(Ba138_ppm[sec_5:crust_loc]))
    
    crust_Li_mmol <- (crust_Li_ppm/MLi)
    crust_B_mmol <- (crust_B_ppm/MB)
    crust_Na_mmol <- (crust_Na_ppm/MNa)
    crust_Al_mmol <- (crust_Al_ppm/MAl)
    crust_Mn_mmol <- (crust_Mn_ppm/MMn)
    crust_Sr_mmol <- (crust_Sr_ppm/MSr)
    crust_Ba135_mmol <- (crust_Ba135_ppm/MBa135)
    crust_Ba137_mmol <- (crust_Ba137_ppm/MBa137)
    crust_Ba138_mmol <- (crust_Ba138_ppm/MBa138)
    
    crust_NaCa <- crust_Na_mmol/crust_Ca43_mol
    crust_BCa <- crust_B_mmol/crust_Ca43_mol
    crust_SrCa <- crust_Sr_mmol/crust_Ca43_mol
    crust_Ba135Ca <- crust_Ba135_mmol/crust_Ca43_mol
    crust_Ba137Ca <- crust_Ba137_mmol/crust_Ca43_mol
    crust_Ba138Ca <- crust_Ba138_mmol/crust_Ca43_mol
    crust_MnCa <- crust_Mn_mmol/crust_Ca43_mol
    crust_AlCa <- crust_Al_mmol/crust_Ca43_mol
    crust_LiCa <- crust_Li_mmol/crust_Ca43_mol
    
    #Lamellar component
    lamellar <- MgCa[crust_loc: length(MgCa)]
    lamellar_Mg <- Mg_ppm[crust_loc: length(Mg_ppm)]
    lamellar_Ca <- Ca43_ppm[crust_loc: length(Ca43_ppm)]
    
    lam_Mg_ppm <- mean(na.omit(lamellar_Mg))
    lam_Ca_ppm <- mean(na.omit(lamellar_Ca))
    
    #convert to mmol/mol
    lam_Ca43_mol <- (lam_Ca_ppm/MCa43)/1000
    lam_Mg_mol <- (lam_Mg_ppm/MMg24)
    
    lam_MgCa_mol <- lam_Mg_mol/lam_Ca43_mol
    
    lam_Li_ppm <- mean(na.omit(Li_ppm[crust_loc: length(Li_ppm)]))
    lam_B_ppm <- mean(na.omit(B_ppm[crust_loc: length(B_ppm)]))
    lam_Na_ppm <- mean(na.omit(Na_ppm[crust_loc: length(Na_ppm)]))
    lam_Al_ppm <- mean(na.omit(Al_ppm[crust_loc: length(Al_ppm)]))
    lam_Mn_ppm <- mean(na.omit(Mn_ppm[crust_loc: length(Mn_ppm)]))
    lam_Sr_ppm <- mean(na.omit(Sr_ppm[crust_loc: length(Sr_ppm)]))
    lam_Ba135_ppm <- mean(na.omit(Ba135_ppm[crust_loc: length(Ba135_ppm)]))
    lam_Ba137_ppm <- mean(na.omit(Ba137_ppm[crust_loc: length(Ba137_ppm)]))
    lam_Ba138_ppm <- mean(na.omit(Ba138_ppm[crust_loc: length(Ba138_ppm)]))
    
    lam_Li_mmol <- (lam_Li_ppm/MLi)
    lam_B_mmol <- (lam_B_ppm/MB)
    lam_Na_mmol <- (lam_Na_ppm/MNa)
    lam_Al_mmol <- (lam_Al_ppm/MAl)
    lam_Mn_mmol <- (lam_Mn_ppm/MMn)
    lam_Sr_mmol <- (lam_Sr_ppm/MSr)
    lam_Ba135_mmol <- (lam_Ba135_ppm/MBa135)
    lam_Ba137_mmol <- (lam_Ba137_ppm/MBa137)
    lam_Ba138_mmol <- (lam_Ba138_ppm/MBa138)
    
    lam_NaCa <- lam_Na_mmol/lam_Ca43_mol
    lam_BCa <- lam_B_mmol/lam_Ca43_mol
    lam_SrCa <- lam_Sr_mmol/lam_Ca43_mol
    lam_Ba135Ca <- lam_Ba135_mmol/lam_Ca43_mol
    lam_Ba137Ca <- lam_Ba137_mmol/lam_Ca43_mol
    lam_Ba138Ca <- lam_Ba138_mmol/lam_Ca43_mol
    lam_MnCa <- lam_Mn_mmol/lam_Ca43_mol
    lam_AlCa <- lam_Al_mmol/lam_Ca43_mol
    lam_LiCa <- lam_Li_mmol/lam_Ca43_mol
    
    #add crust and lamellar value to overview here
    overview[i, 12] <- crust_MgCa_mol #crust
    overview[i, 13] <- lam_MgCa_mol #lamellar
    overview[i, 15] <- crust_NaCa
    overview[i, 16] <- lam_NaCa
    overview[i, 18] <- crust_BCa
    overview[i, 19] <- lam_BCa
    overview[i, 21] <- crust_SrCa
    overview[i, 22] <- lam_SrCa
    overview[i, 24] <- crust_Ba135Ca
    overview[i, 25] <- lam_Ba135Ca
    overview[i, 27] <- crust_Ba137Ca
    overview[i, 28] <- lam_Ba137Ca
    overview[i, 30] <- crust_Ba138Ca
    overview[i, 31] <- lam_Ba138Ca
    overview[i, 33] <- crust_MnCa
    overview[i, 34] <- lam_MnCa
    overview[i, 36] <- crust_AlCa
    overview[i, 37] <- lam_AlCa
    overview[i, 39] <- crust_LiCa
    overview[i, 40] <- lam_LiCa
    
    #Total for the shell wall: combined crust and lamellar calcite
    overview[i, 11] <- totMgCa_mol  # total
    overview[i, 14] <- NaCa
    overview[i, 17] <- BCa
    overview[i, 20] <- SrCa
    overview[i, 23] <- Ba135Ca
    overview[i, 26] <- Ba137Ca
    overview[i, 29] <- Ba138Ca
    overview[i, 32] <- MnCa
    overview[i, 35] <- AlCa
    overview[i, 38] <- LiCa
  }
  
} #End Loop

#### End loop distinguishing crust and lamellar calcite ####

#Save data table as csv file
setwd("") # set working directory for saving data frame of all shells
write.csv(overview, file = "Foram_Overview.csv", row.names = FALSE)

#### END ####
