#libload
packages <- c("tidyverse", "zoo", "imputeTS","stringr","magick","ggplot2")

install.packages(setdiff(packages, rownames(installed.packages()))) 

library(tidyverse)
library(zoo)
library(imputeTS)
library(stringr)
library(magick)
library(grid)
library(ggplot2)
library(fs)

#filesystem_create



folder <- "data"

if (file.exists(folder)) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create(folder)
  
}
folder2 <- "Master FMD Export"

if (file.exists(folder2)) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create(folder2)
  
}
folder4 <- "Master QC Plots"
if (file.exists(folder4)) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create(folder4)
  
}


#dataload labchart files downsample by 333 to match DIA, blockheaderON, Time selected as true
dia_data <- read.delim("data/FMD_DIA.txt", header=F,sep=",", skip=49)
dia_data <- rename(dia_data, "diameter"= "V2", "index"="V1", "participant_id"="V3","condition_id"="V5","reader"="V10")
lc_data <- read.delim("data/FMD_LC.txt", header=F, sep="\t",skip=9)
lc_bl_data <- read.delim("data/BL_LC.txt", header=F, sep="\t",skip=9)
visco_data <- read.csv("data/VISC.csv", header=F,skip=45)
dia_data <- subset(dia_data, select = c("diameter","index","participant_id","condition_id","reader" ))
#file_id
participant.id<- as.character(dia_data[1,3])
participant.id <- str_replace_all(string=participant.id, pattern=" ", repl="")
condition.id<- as.character(dia_data[1,4])
condition.id <- str_replace_all(string=condition.id, pattern=" ", repl="")
reader.id<- as.character(dia_data[1,5])
reader.id <- str_replace_all(string=reader.id, pattern=" ", repl="")
file.id<- paste(participant.id,condition.id)
file.id <- str_replace_all(string=file.id, pattern=" ", repl="_")

fmd_data <- cbind(head(dia_data,5000),head(lc_data,5000))
fmd_data <- rename(fmd_data,"time"="V1","flow_vel"="V2","fing_pres"="V3")

bl_data <- read.delim("data/BL_DIA.txt", header=F,sep=",",skip=49)
bl_data <- rename(bl_data, "diameter"= "V2", "index"="V1")
bl_data <- subset(bl_data, select = c("diameter","index" ))
bl_data <- cbind(head(lc_bl_data,900),head(bl_data,900))
bl_data <- rename(bl_data, "time"="V1", "flow_vel"="V2", "fing_pres"="V3")



visco_data <-cbind(visco_data$V8,visco_data$V4)
visco_data <- as.data.frame(visco_data)
visco_data <- rename(visco_data, "visc_sr"="V1", "visc_visc"="V2")


#output folders
folder3 <- "results"

if (file.exists(file.path(getwd(),"Analyzed",participant.id,file.id,folder3), recursive = TRUE)) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create(file.path(getwd(),"Analyzed",participant.id,file.id,"results"), recursive = T)
  dir.create(file.path(getwd(),"Analyzed",participant.id,file.id,"plots"), recursive = T)
  dir.create(file.path(getwd(),"Analyzed",participant.id,file.id,"data"), recursive = T)
}


#smooth
smo_index <- seq(1, 5000, by = 1) 
smo_dia <- rollmean(fmd_data$diameter, k = 90, fill = NA)
smo_Qvel <- rollmean(fmd_data$flow_vel, k=90, fill= NA)
smo_fing_pres <- rollmean(fmd_data$fing_pres, k=90, fill= NA)
smo_fmd <- cbind(smo_dia,smo_fing_pres,smo_Qvel,smo_index)
smo_fmd <- as.data.frame(smo_fmd)
fmd_clean <- as.data.frame(smo_fmd)
fmd_clean$smo_dia <- na_seadec(fmd_clean$smo_dia)
fmd_clean <- rename(fmd_clean,"diameter"="smo_dia","flow_vel"="smo_Qvel","fing_pres"="smo_fing_pres","index"="smo_index")





#variablecreate

visc_225 <- visco_data[3,2]

#fmd
fmd_clean$shear_rate <- 8*(fmd_clean$flow_vel/fmd_clean$diameter)
fmd_clean$bulk_flow <- fmd_clean$flow_vel*(pi*((fmd_clean$diameter/20)^2))*60
fmd_clean$fvc <- fmd_clean$bulk_flow/fmd_clean$fing_pres
fmd_clean$shear_stress <- fmd_clean$shear_rate*visc_225
 
#bl
bl_data$shear_rate <- 8*(bl_data$flow_vel/bl_data$diameter)
bl_data$bulk_flow <- bl_data$flow_vel*(pi*((bl_data$diameter/20)^2))*60
bl_data$fvc <- bl_data$bulk_flow/bl_data$fing_pres
bl_data$shear_stress <- bl_data$shear_rate*visc_225

#outcomes
bl_diameter <- mean(bl_data$diameter,na.rm=TRUE)
bl_fvc <- mean(bl_data$fvc,na.rm=TRUE)
bl_flow <- mean(bl_data$bulk_flow,na.rm=TRUE)
bl_sr <- mean(bl_data$shear_rate,na.rm=TRUE)
bl_ss <- mean(bl_data$shear_stress,na.rm=TRUE)

peak_diameter <- max(fmd_clean$diameter,na.rm=TRUE)
peak_fvc <- max(fmd_clean$fvc,na.rm=TRUE)
peak_flow <- max(fmd_clean$bulk_flow,na.rm=TRUE)
peak_sr <- max(fmd_clean$shear_rate,na.rm=TRUE)
peak_ss <- max(fmd_clean$shear_stress,na.rm=TRUE)

delta_fvc <- peak_fvc-bl_fvc
delta_flow <- peak_flow-bl_flow

fmd_mm <- peak_diameter-bl_diameter
fmd_per <- fmd_mm/bl_diameter*100

#RESULTS
results <- cbind.data.frame(bl_diameter, bl_fvc, bl_flow, bl_sr, bl_ss,peak_diameter, peak_fvc, peak_flow, peak_sr,peak_ss, fmd_mm,fmd_per,visc_225,participant.id,condition.id,reader.id)

#plotsFMD
dir.create("plots")
png("plots/DIA_BL.png")
plot(bl_data$index,bl_data$diameter,xlab=paste(file.id,"index"),ylab="diameterBL",col = "#2E9FDF",pch=16)
dev.off()
png("plots/SR_BL.png")
plot(bl_data$index,bl_data$shear_rate,xlab=paste(file.id,"index"),ylab="shear rateBL",col = "#2E9FDF",pch=16)
dev.off()
png("plots/FV_BL.png")
plot(bl_data$index,bl_data$flow_vel,xlab=paste(file.id,"index"),ylab="flow velocityBL",col = "#2E9FDF",pch=16)
dev.off()
png("plots/SS_BL.png")
plot(bl_data$index,bl_data$shear_stress,xlab=paste(file.id,"index"),ylab="shear sressPO",col = "#2E9FDF",pch=16)
dev.off()

png("plots/DIA_PO.png")
plot(fmd_clean$index,fmd_clean$diameter,xlab=paste(file.id,"index"),ylab="diameterPO",col = "#FF5733",pch=16)
dev.off()
png("plots/SR_PO.png")
plot(fmd_clean$index,fmd_clean$shear_rate,xlab=paste(file.id,"index"),ylab="shear ratePO",col = "#FF5733",pch=16)
dev.off()
png("plots/FV_PO.png")
plot(fmd_clean$index,fmd_clean$flow_vel,xlab=paste(file.id,"index"),ylab="flow velocityPO",col = "#FF5733",pch=16)
dev.off()
png("plots/SS_PO.png")
plot(fmd_clean$index,fmd_clean$shear_stress,xlab=paste(file.id,"index"),ylab="shear stressPO",col = "#FF5733",pch=16)
dev.off()
#combine and report
# Load the six images into R using the png package
img1 <- image_read("plots/DIA_BL.png")
img2 <- image_read("plots/SR_BL.png")
img3 <- image_read("plots/FV_BL.png")
img4 <- image_read("plots/SS_BL.png")
img5 <- image_read("plots/DIA_PO.png")
img6 <- image_read("plots/SR_PO.png")
img7 <- image_read("plots/FV_PO.png")
img8 <- image_read("plots/SS_PO.png")

# Combine the images into a single image
combined_img1 <- image_append(c(img1, img2, img3, img4), stack = FALSE)
combined_img2 <- image_append(c(img5, img6, img7, img8), stack = FALSE)
combined_img <- image_append(c(combined_img1,combined_img2), stack=TRUE)
# Save the combined image as a PNG file
image_write(combined_img, "QCPlots.png")
image_write(combined_img, file.path(getwd(),"Analyzed",participant.id,file.id,"plots","QCPlots.png"))
image_write(combined_img, file.path(getwd(),folder4,paste0(file.id,".png")))

#raw_data_output
if (file.exists(file.path(getwd(),"Analyzed",participant.id,file.id,"data","BL_DIA.txt")) )
{
  
  cat("The folder already exists")
  
} else {
data_files <- list.files("data")
file_copy(file.path(getwd(),"data",data_files),file.path(getwd(),"Analyzed",participant.id,file.id,"data"),data_files)
}

#results_print csv
write.csv(results,file.path(getwd(),"Analyzed",participant.id,file.id,"results","FMD_results.csv"), row.names=FALSE)
write.csv(bl_data,file.path(getwd(),"Analyzed",participant.id,file.id,"results","BL_data_clean.csv"), row.names=FALSE)
write.csv(fmd_clean,file.path(getwd(),"Analyzed",participant.id,file.id,"results","PO_data_clean.csv"), row.names=FALSE)



if (file.exists("Master FMD Export/MASTER_FMD_DATA.csv") ){
  
  write.table( results,  
               file="Master FMD Export/MASTER_FMD_DATA.csv", 
               append = T, 
               sep=',', 
               row.names=F, 
               col.names=F )
  
  
} else {
  write.table( results,  
               file="Master FMD Export/MASTER_FMD_DATA.csv", 
               append = T, 
               sep=',', 
               row.names=F, 
               col.names=T )
}


