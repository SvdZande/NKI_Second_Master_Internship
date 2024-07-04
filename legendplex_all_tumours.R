#analysis section for legendplex over all tumours

source("~/Documents/R scripts/Functions.R")
library(readxl)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(forstringr)


#--------
#read in normalised data
lu023_TLS014 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/LU023_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu023_TLS014) <- lu023_TLS014[,1]
lu023_TLS014[,1] <- NULL
rownames(lu023_TLS014) <- gsub("Unstim", "unstim", rownames(lu023_TLS014))

lu088 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/LU088_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu088) <- lu088[,1]
lu088[,1] <- NULL
rownames(lu088) <- gsub("Unstim", "unstim", rownames(lu088))

lu091 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/LU091_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu091) <- lu091[,1]
lu091[,1] <- NULL
rownames(lu091) <- gsub("Unstim", "unstim", rownames(lu091))

lu035_TLS013 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS013_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu035_TLS013) <- lu035_TLS013[,1]
lu035_TLS013[,1] <- NULL

lu035_TLS015 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS015_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu035_TLS015) <- lu035_TLS015[,1]
lu035_TLS015[,1] <- NULL

lu066_TLS015 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU066_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu066_TLS015) <- lu066_TLS015[,1]
lu066_TLS015[,1] <- NULL

lu052_TLS010 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS010_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu052_TLS010) <- lu052_TLS010[,1]
lu052_TLS010[,1] <- NULL
rownames(lu052_TLS010) <- gsub("hrIFNy", "hrIFNy_high", rownames(lu052_TLS010))

lu052_TLS014 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS014_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu052_TLS014) <- lu052_TLS014[,1]
lu052_TLS014[,1] <- NULL


lu023_TLS013 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LU023_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu023_TLS013) <- lu023_TLS013[,1]
lu023_TLS013[,1] <- NULL

lu086_2 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LU086-2_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu086_2) <- lu086_2[,1]
lu086_2[,1] <- NULL

lu083_TLS015 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LU083_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu083_TLS015) <- lu083_TLS015[,1]
lu083_TLS015[,1] <- NULL

lu086_1 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU086-1_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu086_1) <- lu086_1[,1]
lu086_1[,1] <- NULL

lu066_TLS012 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU066_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu066_TLS012) <- lu066_TLS012[,1]
lu066_TLS012[,1] <- NULL

lu083_TLS012 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU083_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu083_TLS012) <- lu083_TLS012[,1]
lu083_TLS012[,1] <- NULL

lu080 <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU080_normalised_windsor.xlsx", sheet = 3)) #corrected values
rownames(lu080) <- lu080[,1]
lu080[,1] <- NULL

#read in TLS information
TLS_status <- read.csv("~/Documents/Analyses/TLS/fragment_table.csv")
TLS_status$X <- NULL
TLS_status$heatmap <- paste0(TLS_status$tumour, " ", TLS_status$condition, " fragment " ,TLS_status$fragment)
TLS_status$heatmap <- gsub("aIFNyR1\\+aPD1", "aPD1+aIFNyR1", TLS_status$heatmap)
TLS_status$heatmap[TLS_status$experiment == "TLS013" & TLS_status$tumour == "LU035"] <- paste0("TLS013 ", TLS_status$heatmap[TLS_status$experiment == "TLS013" & TLS_status$tumour == "LU035"])
TLS_status$heatmap[TLS_status$experiment == "TLS015" & TLS_status$tumour == "LU035"] <- paste0("TLS015 ", TLS_status$heatmap[TLS_status$experiment == "TLS015" & TLS_status$tumour == "LU035"])
TLS_status$heatmap[TLS_status$experiment == "TLS010" & TLS_status$tumour == "LU052"] <- paste0("TLS010 ", TLS_status$heatmap[TLS_status$experiment == "TLS010" & TLS_status$tumour == "LU052"])
TLS_status$heatmap[TLS_status$experiment == "TLS014" & TLS_status$tumour == "LU052"] <- paste0("TLS014 ", TLS_status$heatmap[TLS_status$experiment == "TLS014" & TLS_status$tumour == "LU052"])

#denote if LU086 is replicate 1 or 2
TLS_status$heatmap[TLS_status$experiment == "TLS011" & grepl("LU086_1",TLS_status$Name)] <- gsub("LU086", "LU086-1", TLS_status$heatmap[TLS_status$experiment == "TLS011" & grepl("LU086_1",TLS_status$Name)])
TLS_status$heatmap[TLS_status$experiment == "TLS011" & grepl("LU086_2",TLS_status$Name)] <- gsub("LU086", "LU086-2", TLS_status$heatmap[TLS_status$experiment == "TLS011" & grepl("LU086_2",TLS_status$Name)])

#remove tumours that were not assessed (at this point only TLS009 LU086)
TLS_status <- TLS_status[TLS_status$experiment != "TLS009",]

#remove the RN440-15 and RN440-16 conditions
TLS_status <- TLS_status[!(TLS_status$condition %in% c("RN440-15", "RN440-16")),]

#Get TLS segments out
TLS_023_013 <- TLS_status[TLS_status$tumour == "LU023"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS013","heatmap"]
int_023_013 <- TLS_status[TLS_status$tumour == "LU023"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS013","heatmap"]
tumour_023_013 <- TLS_status[TLS_status$tumour == "LU023"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS013","heatmap"]

TLS_023_014 <- TLS_status[TLS_status$tumour == "LU023"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS014","heatmap"]
int_023_014 <- TLS_status[TLS_status$tumour == "LU023"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS014","heatmap"]
tumour_023_014 <- TLS_status[TLS_status$tumour == "LU023"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS014","heatmap"]

TLS_035_013 <- TLS_status[TLS_status$tumour == "LU035"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS013","heatmap"]
int_035_013 <- TLS_status[TLS_status$tumour == "LU035"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS013","heatmap"]
tumour_035_013 <- TLS_status[TLS_status$tumour == "LU035"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS013","heatmap"]

TLS_035_015 <- TLS_status[TLS_status$tumour == "LU035"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS015","heatmap"]
int_035_015 <- TLS_status[TLS_status$tumour == "LU035"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS015","heatmap"]
tumour_035_015 <- TLS_status[TLS_status$tumour == "LU035"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS015","heatmap"]

TLS_052_010 <- TLS_status[TLS_status$tumour == "LU052"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS010","heatmap"]
int_052_010 <- TLS_status[TLS_status$tumour == "LU052"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS010","heatmap"]
tumour_052_010 <- TLS_status[TLS_status$tumour == "LU052"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS010","heatmap"]

TLS_052_014 <- TLS_status[TLS_status$tumour == "LU052"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS014","heatmap"]
int_052_014 <- TLS_status[TLS_status$tumour == "LU052"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS014","heatmap"]
tumour_052_014 <- TLS_status[TLS_status$tumour == "LU052"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS014","heatmap"]

TLS_066_015 <- TLS_status[TLS_status$tumour == "LU066"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS015","heatmap"]
int_066_015 <- TLS_status[TLS_status$tumour == "LU066"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS015","heatmap"]
tumour_066_015 <- TLS_status[TLS_status$tumour == "LU066"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS015","heatmap"]

TLS_066_012 <- TLS_status[TLS_status$tumour == "LU066"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS012","heatmap"]
int_066_012 <- TLS_status[TLS_status$tumour == "LU066"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS012","heatmap"]
tumour_066_012 <- TLS_status[TLS_status$tumour == "LU066"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS012","heatmap"]

TLS_080 <- TLS_status[TLS_status$tumour == "LU080" & TLS_status$segment == "TLS" ,"heatmap"]
int_080 <- TLS_status[TLS_status$tumour == "LU080"& TLS_status$segment == "intermediate","heatmap"]
tumour_080 <- TLS_status[TLS_status$tumour == "LU080"& TLS_status$segment == "tumour_bed","heatmap"]

TLS_083_012 <- TLS_status[TLS_status$tumour == "LU083"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS012","heatmap"]
int_083_012 <- TLS_status[TLS_status$tumour == "LU083"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS012","heatmap"]
tumour_083_012 <- TLS_status[TLS_status$tumour == "LU083"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS012","heatmap"]

TLS_083_015 <- TLS_status[TLS_status$tumour == "LU083"& TLS_status$segment == "TLS" & TLS_status$experiment == "TLS015","heatmap"]
int_083_015 <- TLS_status[TLS_status$tumour == "LU083"& TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS015","heatmap"]
tumour_083_015 <- TLS_status[TLS_status$tumour == "LU083"& TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS015","heatmap"]

TLS_086_1 <- TLS_status[grepl("LU086_1", TLS_status$Name) & TLS_status$segment == "TLS" & TLS_status$experiment == "TLS011","heatmap"]
int_086_1 <- TLS_status[grepl("LU086_1", TLS_status$Name) & TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS011","heatmap"]
tumour_086_1 <- TLS_status[grepl("LU086_1", TLS_status$Name) & TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS011","heatmap"]

TLS_086_2 <- TLS_status[grepl("LU086_2", TLS_status$Name) & TLS_status$segment == "TLS" & TLS_status$experiment == "TLS011","heatmap"]
int_086_2 <- TLS_status[grepl("LU086_2", TLS_status$Name) & TLS_status$segment == "intermediate" & TLS_status$experiment == "TLS011","heatmap"]
tumour_086_2 <- TLS_status[grepl("LU086_2", TLS_status$Name) & TLS_status$segment == "tumour_bed" & TLS_status$experiment == "TLS011","heatmap"]

TLS_088 <- TLS_status[TLS_status$tumour == "LU088" & TLS_status$segment == "TLS" ,"heatmap"]
int_088 <- TLS_status[TLS_status$tumour == "LU088"& TLS_status$segment == "intermediate","heatmap"]
tumour_088 <- TLS_status[TLS_status$tumour == "LU088"& TLS_status$segment == "tumour_bed","heatmap"]

TLS_091 <- TLS_status[TLS_status$tumour == "LU091" & TLS_status$segment == "TLS" ,"heatmap"]
int_091 <- TLS_status[TLS_status$tumour == "LU091"& TLS_status$segment == "intermediate","heatmap"]
tumour_091 <- TLS_status[TLS_status$tumour == "LU091"& TLS_status$segment == "tumour_bed","heatmap"]

#---------------------------
#Run QCs (from source file Functions.R)
file_list <- c("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/LU023.xlsx", "~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/LU088.xlsx", 
"~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/LU091.xlsx", "~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS013.xlsx", 
"~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS015.xlsx", "~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU066.xlsx", 
"~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS010.xlsx", "~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS014.xlsx",
"~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LU023_TLS013.xlsx","~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LU083_TLS015.xlsx",
"~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LU086-2_TLS011.xlsx","~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU066.xlsx",
"~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU080.xlsx","~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU083.xlsx",
"~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LU086-1.xlsx") 

MFI_file_list <- c("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/MFI LP063.xlsx","~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/MFI LP063.xlsx",
                   "~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023/MFI LP063.xlsx","~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/MFI LP064.xlsx",
                   "~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/MFI LP064.xlsx","~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/MFI LP065.xlsx",
                   "~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/MFI LP065.xlsx","~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/MFI LP065.xlsx",
                   "~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LP066_MFI.xlsx","~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LP066_MFI.xlsx",
                   "~/Documents/Experiments/TLS/LP066-LU023-LU083-LU086/LP066_MFI.xlsx","~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LP060_MFI.xlsx",
                   "~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LP060_MFI.xlsx","~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LP060_MFI.xlsx",
                   "~/Documents/Experiments/TLS/LP060-RW_PK_01042024/LP060_MFI.xlsx")

tumours <- c("LU023", "LU088", "LU091", "LU035", "LU035", "LU066", "LU052", "LU052", "LU023", "LU083", "LU086-2", "LU066", "LU080", "LU083", "LU086-1")
experiments <- c("TLS014", "TLS016", "TLS016", "TLS013", "TLS015", "TLS015", "TLS010", "TLS014" ,"TLS013","TLS015","TLS011","TLS012","TLS011","TLS012","TLS011")
legendplex <- c("LP063","LP063","LP063","LP064","LP064","LP065","LP065","LP065","LP066","LP066","LP066","LP060","LP060","LP060","LP060")

conc_plots <- list()
MFI_plots <- list()
stat_list <- list()

for (i in seq_along(file_list)){
  pdf(paste0("~/Documents/Analyses/TLS/QC plots/QC_", tumours[[i]], "_",experiments[[i]],".pdf"), width = 15)
  print((paste0(tumours[[i]], " ", experiments[[i]], " report")))
  
  cat("Based on concentration:","\n")
  conc_plots[[length(conc_plots)+1]] <- legendplex_qc_concentration(file_list[[i]])
  names(conc_plots)[[length(conc_plots)]] <- paste0(tumours[[i]], "_", experiments[[i]])
  print(conc_plots[[i]])
  
  cat("Based on MFI:","\n")
  MFI_plots[[length(MFI_plots)+1]] <- legendplex_qc_MFI(MFI_file = MFI_file_list[[i]], experiment = experiments[[i]], tumour_name = tumours[[i]])
  names(MFI_plots)[[length(MFI_plots)]] <- paste0(tumours[[i]], "_", experiments[[i]])
  print(MFI_plots[[i]])
  
  dev.off()
  
  stat_list[[length(stat_list)+1]] <- legendplex_fragment_control(file_list[[i]], experiment = experiments[[i]])
  names(stat_list)[[length(stat_list)]] <- paste0(tumours[[i]], "_", experiments[[i]])
}

#----------------------------

#Get separate Legendplex scores per segment
TLS_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% TLS_023_013,]
TLS_023_013_data <- colMeans(TLS_023_013_data) #average score per cytokine across all conditions
int_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% int_023_013,]
int_023_013_data <- colMeans(int_023_013_data)
tumour_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% tumour_023_013,]
tumour_023_013_data <- colMeans(tumour_023_013_data)

TLS_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% TLS_023_014,]
TLS_023_014_data <- colMeans(TLS_023_014_data) #average score per cytokine across all conditions
int_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% int_023_014,]
int_023_014_data <- colMeans(int_023_014_data)
tumour_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% tumour_023_014,]
tumour_023_014_data <- colMeans(tumour_023_014_data)

TLS_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% TLS_035_013,]
TLS_035_013_data <- colMeans(TLS_035_013_data) #average score per cytokine across all conditions
int_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% int_035_013,]
int_035_013_data <- colMeans(int_035_013_data)
tumour_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% tumour_035_013,]
tumour_035_013_data <- colMeans(tumour_035_013_data)

TLS_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% TLS_035_015,]
TLS_035_015_data <- colMeans(TLS_035_015_data) #average score per cytokine across all conditions
int_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% int_035_015,]
int_035_015_data <- colMeans(int_035_015_data)
tumour_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% tumour_035_015,]
tumour_035_015_data <- colMeans(tumour_035_015_data)

TLS_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% TLS_052_010,]
TLS_052_010_data <- colMeans(TLS_052_010_data) #average score per cytokine across all conditions
int_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% int_052_010,]
int_052_010_data <- colMeans(int_052_010_data)
tumour_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% tumour_052_010,]
tumour_052_010_data <- colMeans(tumour_052_010_data)

TLS_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% TLS_052_014,]
TLS_052_014_data <- colMeans(TLS_052_014_data) #average score per cytokine across all conditions
int_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% int_052_014,]
int_052_014_data <- colMeans(int_052_014_data)
tumour_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% tumour_052_014,]
tumour_052_014_data <- colMeans(tumour_052_014_data)

#This tumour has no intermediate fragments. Therefore, this condition is skipped.
TLS_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% TLS_066_012,]
TLS_066_012_data <- colMeans(TLS_066_012_data) #average score per cytokine across all conditions
tumour_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% tumour_066_012,]
tumour_066_012_data <- colMeans(tumour_066_012_data)

TLS_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% TLS_066_015,]
TLS_066_015_data <- colMeans(TLS_066_015_data) #average score per cytokine across all conditions
int_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% int_066_015,]
int_066_015_data <- colMeans(int_066_015_data)
tumour_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% tumour_066_015,]
tumour_066_015_data <- colMeans(tumour_066_015_data)

#LU080 has no TLS so we ignore this condition
int_080_data <- lu080[rownames(lu080) %in% int_080,]
int_080_data <- colMeans(int_080_data)
tumour_080_data <- lu080[rownames(lu080) %in% tumour_080,]
tumour_080_data <- colMeans(tumour_080_data)

TLS_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% TLS_083_012,]
TLS_083_012_data <- colMeans(TLS_083_012_data) #average score per cytokine across all conditions
int_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% int_083_012,]
int_083_012_data <- colMeans(int_083_012_data)
tumour_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% tumour_083_012,]
tumour_083_012_data <- colMeans(tumour_083_012_data)

TLS_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% TLS_083_015,]
TLS_083_015_data <- colMeans(TLS_083_015_data) #average score per cytokine across all conditions
int_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% int_083_015,]
int_083_015_data <- colMeans(int_083_015_data)
tumour_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% tumour_083_015,]
tumour_083_015_data <- colMeans(tumour_083_015_data)

#the only TLS of LU086-1 were in RN440-15/RN440-16, which will not be considered. Therefore, we ignore this condition.
int_086_1_data <- lu086_1[rownames(lu086_1) %in% int_086_1,]
int_086_1_data <- colMeans(int_086_1_data)
tumour_086_1_data <- lu086_1[rownames(lu086_1) %in% tumour_086_1,]
tumour_086_1_data <- colMeans(tumour_086_1_data)

#the only TLS of LU086-2 was in RN440-15, which was not taken along in the Legendplex. Therefore, we ignore this condition.
int_086_2_data <- lu086_2[rownames(lu086_2) %in% int_086_2,]
int_086_2_data <- colMeans(int_086_2_data)
tumour_086_2_data <- lu086_2[rownames(lu086_2) %in% tumour_086_2,]
tumour_086_2_data <- colMeans(tumour_086_2_data)

TLS_088_data <- lu088[rownames(lu088) %in% TLS_088,]
TLS_088_data <- colMeans(TLS_088_data) #average score per cytokine across all conditions
int_088_data <- lu088[rownames(lu088) %in% int_088,]
int_088_data <- colMeans(int_088_data)
tumour_088_data <- lu088[rownames(lu088) %in% tumour_088,]
tumour_088_data <- colMeans(tumour_088_data)

TLS_091_data <- lu091[rownames(lu091) %in% TLS_091,]
TLS_091_data <- colMeans(TLS_091_data) #average score per cytokine across all conditions
int_091_data <- lu091[rownames(lu091) %in% int_091,]
int_091_data <- colMeans(int_091_data)
tumour_091_data <- lu091[rownames(lu091) %in% tumour_091,]
tumour_091_data <- colMeans(tumour_091_data)

#Average across all tumours
TLS_scores <- rbind(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                    TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data) #missing LU080, LU086-1, LU086-2
TLS_scores <- colMeans(TLS_scores)

int_scores <- rbind(int_023_013_data,int_023_014_data, int_035_013_data, int_035_015_data, int_052_010_data, int_052_014_data, int_066_015_data,
                    int_080_data, int_083_012_data, int_083_015_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data) #missing LU066 (TLS012)
int_scores <- colMeans(int_scores)

tumour_scores <- rbind(tumour_023_013_data,tumour_023_014_data, tumour_035_013_data, tumour_035_015_data, tumour_052_010_data, tumour_052_014_data, tumour_066_012_data,tumour_066_015_data,
                       tumour_080_data, tumour_083_012_data, tumour_083_015_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data)
tumour_scores <- colMeans(tumour_scores)

heatmap_scores <- rbind(TLS_scores, int_scores, tumour_scores)
rownames(heatmap_scores) <- c("TLS", "intermediate", "tumour bed")

#Generate a heatmap of z-scores per cytokine
library(pheatmap)
  #z-score is (value - average)/standard deviation. (is this the delta or is the average per condition?)
  colSD <- apply(heatmap_scores, 2, sd)
  colMean <- colMeans(heatmap_scores)
  z_scores <- sweep(heatmap_scores, 2, colMean, "-") #subtract means per column
  z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
  z_scores <- t(z_scores)
  
  rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"
  rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y"
  
  #heatmap
  heatmap <- suppressWarnings(pheatmap(z_scores))  
  heatmap
  
  
  save_pheatmap(heatmap, "~/Documents/Analyses/TLS/final_all_fragments/segment_normalised_windsor.pdf")
  
  
 #------------------------- 
  #split per condition
  TLS_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% TLS_023_013,]
  int_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% int_023_013,]
  tumour_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% tumour_023_013,]

  TLS_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% TLS_023_014,]
  int_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% int_023_014,]
  tumour_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% tumour_023_014,]

  TLS_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% TLS_035_013,]
  int_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% int_035_013,]
  tumour_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% tumour_035_013,]

  TLS_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% TLS_035_015,]
  int_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% int_035_015,]
  tumour_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% tumour_035_015,]

  TLS_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% TLS_052_010,]
  int_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% int_052_010,]
  tumour_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% tumour_052_010,]
  #LU052_TLS010_tumour_bed has extremely high CXCL9/CXCL10, what does the heatmap look like if I take the IFN_high fragments out?
  TLS_052_010_data <- TLS_052_010_data[!grepl("hrIFNy_high", rownames(TLS_052_010_data)),]
  int_052_010_data <- int_052_010_data[!grepl("hrIFNy_high", rownames(int_052_010_data)),]
  tumour_052_010_data <- tumour_052_010_data[!grepl("hrIFNy_high", rownames(tumour_052_010_data)),]
  
  TLS_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% TLS_052_014,]
  int_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% int_052_014,]
  tumour_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% tumour_052_014,]

  #This tumour has no intermediate fragments. Therefore, this condition is skipped.
  TLS_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% TLS_066_012,]
  tumour_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% tumour_066_012,]

  TLS_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% TLS_066_015,]
  int_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% int_066_015,]
  tumour_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% tumour_066_015,]

  #LU080 has no TLS so we ignore this condition
  int_080_data <- lu080[rownames(lu080) %in% int_080,]
  tumour_080_data <- lu080[rownames(lu080) %in% tumour_080,]

  TLS_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% TLS_083_012,]
  int_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% int_083_012,]
  tumour_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% tumour_083_012,]

  TLS_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% TLS_083_015,]
  int_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% int_083_015,]
  tumour_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% tumour_083_015,]

  #the only TLS of LU086-1 were in RN440-15/RN440-16, which will not be considered. Therefore, we ignore this condition.
  int_086_1_data <- lu086_1[rownames(lu086_1) %in% int_086_1,]
  tumour_086_1_data <- lu086_1[rownames(lu086_1) %in% tumour_086_1,]

  #the only TLS of LU086-2 was in RN440-15, which was not taken along in the Legendplex. Therefore, we ignore this condition.
  int_086_2_data <- lu086_2[rownames(lu086_2) %in% int_086_2,]
  tumour_086_2_data <- lu086_2[rownames(lu086_2) %in% tumour_086_2,]

  TLS_088_data <- lu088[rownames(lu088) %in% TLS_088,]
  int_088_data <- lu088[rownames(lu088) %in% int_088,]
  tumour_088_data <- lu088[rownames(lu088) %in% tumour_088,]

  TLS_091_data <- lu091[rownames(lu091) %in% TLS_091,]
  int_091_data <- lu091[rownames(lu091) %in% int_091,]
  tumour_091_data <- lu091[rownames(lu091) %in% tumour_091,]

  #First for TLS
  TLS_list <- list(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                   TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data) #missing LU080, LU086-1, LU086-2
  names(TLS_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012","LU066_TLS015",
                       "LU083_TLS012", "LU083_TLS015", "LU088", "LU091")  
  
  unstim_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Calculate condition averages per tumour
  for(i in seq_along(TLS_list)){
    data <- TLS_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
    unstim_average_TLS <- rbind(unstim_average_TLS, unstim)
    rownames(unstim_average_TLS)[[nrow(unstim_average_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_average_TLS <- rbind(aPD1_average_TLS, aPD1)
      rownames(aPD1_average_TLS)[[nrow(aPD1_average_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_average_TLS <- rbind(aPD1_aIFNR_average_TLS, aPD1_aIFNR)
      rownames(aPD1_aIFNR_average_TLS)[[nrow(aPD1_aIFNR_average_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_average_TLS <- rbind(hrIFNy_average_TLS, hrIFNy)
      rownames(hrIFNy_average_TLS)[[nrow(hrIFNy_average_TLS)]] <- names(TLS_list)[[i]]
    }

  }
  
  #Average across all tumours
  colnames(unstim_average_TLS) <- colnames(TLS_023_014_data)
  unstim_average_TLS <- colMeans(unstim_average_TLS)
  
  colnames(aPD1_average_TLS) <- colnames(TLS_023_014_data)
  aPD1_average_TLS <- colMeans(aPD1_average_TLS)
  
  colnames(aPD1_aIFNR_average_TLS) <- colnames(TLS_023_014_data)
  aPD1_aIFNR_average_TLS <- colMeans(aPD1_aIFNR_average_TLS)
  
  colnames(hrIFNy_average_TLS) <- colnames(TLS_023_014_data)
  hrIFNy_average_TLS <- colMeans(hrIFNy_average_TLS)
  
  #repeat for intermediate and tumour sections
  int_list <- list(int_023_013_data,int_023_014_data, int_035_013_data, int_035_015_data, int_052_010_data, int_052_014_data, int_066_015_data,
                   int_080_data, int_083_012_data, int_083_015_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data) #missing LU066 (TLS012)
  names(int_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014", "LU066_TLS015",
                       "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
  unstim_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Calculate condition averages per tumour
  for(i in seq_along(int_list)){
    data <- int_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
    unstim_average_int <- rbind(unstim_average_int, unstim)
    rownames(unstim_average_int)[[nrow(unstim_average_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_average_int <- rbind(aPD1_average_int, aPD1)
      rownames(aPD1_average_int)[[nrow(aPD1_average_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_average_int <- rbind(aPD1_aIFNR_average_int, aPD1_aIFNR)
      rownames(aPD1_aIFNR_average_int)[[nrow(aPD1_aIFNR_average_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_average_int <- rbind(hrIFNy_average_int, hrIFNy)
      rownames(hrIFNy_average_int)[[nrow(hrIFNy_average_int)]] <- names(int_list)[[i]]
    }
    
  }
  
  #Average across all tumours
  colnames(unstim_average_int) <- colnames(TLS_023_014_data)
  unstim_average_int <- colMeans(unstim_average_int)
  
  colnames(aPD1_average_int) <- colnames(TLS_023_014_data)
  aPD1_average_int <- colMeans(aPD1_average_int)
  
  colnames(aPD1_aIFNR_average_int) <- colnames(TLS_023_014_data)
  aPD1_aIFNR_average_int <- colMeans(aPD1_aIFNR_average_int)
  
  colnames(hrIFNy_average_int) <- colnames(TLS_023_014_data)
  hrIFNy_average_int <- colMeans(hrIFNy_average_int)
  
  #and tumour
  tumour_list <- list(tumour_023_013_data,tumour_023_014_data, tumour_035_013_data, tumour_035_015_data, tumour_052_010_data, tumour_052_014_data,tumour_066_012_data,
                      tumour_066_015_data, tumour_080_data, tumour_083_012_data, tumour_083_015_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data)
  names(tumour_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012", "LU066_TLS015",
                       "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
  unstim_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Calculate condition averages per tumour
  for(i in seq_along(tumour_list)){
    data <- tumour_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
    unstim_average_tumour <- rbind(unstim_average_tumour, unstim)
    rownames(unstim_average_tumour)[[nrow(unstim_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_average_tumour <- rbind(aPD1_average_tumour, aPD1)
      rownames(aPD1_average_tumour)[[nrow(aPD1_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_average_tumour <- rbind(aPD1_aIFNR_average_tumour, aPD1_aIFNR)
      rownames(aPD1_aIFNR_average_tumour)[[nrow(aPD1_aIFNR_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_average_tumour <- rbind(hrIFNy_average_tumour, hrIFNy)
      rownames(hrIFNy_average_tumour)[[nrow(hrIFNy_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
  }
  
  #Average across all tumours
  colnames(unstim_average_tumour) <- colnames(TLS_023_014_data)
  unstim_average_tumour <- colMeans(unstim_average_tumour)
  
  colnames(aPD1_average_tumour) <- colnames(TLS_023_014_data)
  aPD1_average_tumour <- colMeans(aPD1_average_tumour)
  
  colnames(aPD1_aIFNR_average_tumour) <- colnames(TLS_023_014_data)
  aPD1_aIFNR_average_tumour <- colMeans(aPD1_aIFNR_average_tumour)
  
  colnames(hrIFNy_average_tumour) <- colnames(TLS_023_014_data)
  hrIFNy_average_tumour <- colMeans(hrIFNy_average_tumour)
  
  
  
  
  #combine in a heatmap
  heatmap <- rbind(unstim_average_TLS,unstim_average_int,unstim_average_tumour,aPD1_average_TLS,aPD1_average_int,aPD1_average_tumour,
                   aPD1_aIFNR_average_TLS,aPD1_aIFNR_average_int,aPD1_aIFNR_average_tumour,hrIFNy_average_TLS,hrIFNy_average_int,hrIFNy_average_tumour)

  rownames(heatmap) <- c("unstim TLS" ,"unstim intermediate","unstim tumour bed","aPD1 TLS","aPD1 intermediate","aPD1 tumour bed",
                         "aPD1+aIFNyR1 TLS","aPD1+aIFNyR1 intermediate","aPD1+aIFNyR1 tumour bed","hrIFNy TLS","hrIFNy intermediate","hrIFNy tumour bed")  

  #z-scoring
  colSD <- apply(heatmap, 2, sd)
  colMean <- colMeans(heatmap)
  z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
  z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
  z_scores <- t(z_scores)

  rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y" #to make sure rownames can be seen
  rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"
  
  #annotation
  annotation <- data.frame(segment = rep(c("TLS", "intermediate", "tumor bed"),4), condition = rep(c("unstim","aPD1","aPD1+aIFNyR1","hrIFNy"), each = 3))
  rownames(annotation) <- colnames(z_scores)
  
  #heatmap
  heatmap <- pheatmap(z_scores, annotation_col = annotation,cluster_cols = F)
  
  save_pheatmap(heatmap, "~/Documents/Analyses/TLS/final_all_fragments/segment_condition_heatmap_normalised_windsor.pdf")
  
 #-------------------------------
  #Per tumour
  TLS_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% TLS_023_013,]
  int_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% int_023_013,]
  tumour_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% tumour_023_013,]
  
  TLS_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% TLS_023_014,]
  int_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% int_023_014,]
  tumour_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% tumour_023_014,]
  
  TLS_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% TLS_035_013,]
  int_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% int_035_013,]
  tumour_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% tumour_035_013,]
  
  TLS_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% TLS_035_015,]
  int_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% int_035_015,]
  tumour_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% tumour_035_015,]
  
  TLS_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% TLS_052_010,]
  int_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% int_052_010,]
  tumour_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% tumour_052_010,]
  #excluding IFNy_high fragments:
  TLS_052_010_data <- TLS_052_010_data[!grepl("hrIFNy_high", rownames(TLS_052_010_data)),]
  int_052_010_data <- int_052_010_data[!grepl("hrIFNy_high", rownames(int_052_010_data)),]
  tumour_052_010_data <- tumour_052_010_data[!grepl("hrIFNy_high", rownames(tumour_052_010_data)),]
  
  TLS_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% TLS_052_014,]
  int_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% int_052_014,]
  tumour_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% tumour_052_014,]
  
  #This tumour has no intermediate fragments. Therefore, this condition is skipped.
  TLS_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% TLS_066_012,]
  tumour_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% tumour_066_012,]
  
  TLS_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% TLS_066_015,]
  int_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% int_066_015,]
  tumour_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% tumour_066_015,]
  
  #LU080 has no TLS so we ignore this condition
  int_080_data <- lu080[rownames(lu080) %in% int_080,]
  tumour_080_data <- lu080[rownames(lu080) %in% tumour_080,]
  
  TLS_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% TLS_083_012,]
  int_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% int_083_012,]
  tumour_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% tumour_083_012,]
  
  TLS_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% TLS_083_015,]
  int_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% int_083_015,]
  tumour_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% tumour_083_015,]
  
  #the only TLS of LU086-1 were in RN440-15/RN440-16, which will not be considered. Therefore, we ignore this condition.
  int_086_1_data <- lu086_1[rownames(lu086_1) %in% int_086_1,]
  tumour_086_1_data <- lu086_1[rownames(lu086_1) %in% tumour_086_1,]
  
  #the only TLS of LU086-2 was in RN440-15, which was not taken along in the Legendplex. Therefore, we ignore this condition.
  int_086_2_data <- lu086_2[rownames(lu086_2) %in% int_086_2,]
  tumour_086_2_data <- lu086_2[rownames(lu086_2) %in% tumour_086_2,]
  
  TLS_088_data <- lu088[rownames(lu088) %in% TLS_088,]
  int_088_data <- lu088[rownames(lu088) %in% int_088,]
  tumour_088_data <- lu088[rownames(lu088) %in% tumour_088,]
  
  TLS_091_data <- lu091[rownames(lu091) %in% TLS_091,]
  int_091_data <- lu091[rownames(lu091) %in% int_091,]
  tumour_091_data <- lu091[rownames(lu091) %in% tumour_091,]
  
  
  TLS_per_tumour_data <- list(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                              TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data)
  names(TLS_per_tumour_data) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012","LU066_TLS015",
                                  "LU083_TLS012", "LU083_TLS015", "LU088", "LU091")  
  
  unstim_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  
  for(i in seq_along(TLS_per_tumour_data)){
    data <- TLS_per_tumour_data[[i]]
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
      unstim_per_tumour_TLS <- rbind(unstim_per_tumour_TLS, unstim)
      rownames(unstim_per_tumour_TLS)[[nrow(unstim_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_per_tumour_TLS <- rbind(aPD1_per_tumour_TLS, aPD1)
      rownames(aPD1_per_tumour_TLS)[[nrow(aPD1_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_per_tumour_TLS <- rbind(aPD1_aIFNR_per_tumour_TLS, aPD1_aIFNR)
      rownames(aPD1_aIFNR_per_tumour_TLS)[[nrow(aPD1_aIFNR_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_per_tumour_TLS <- rbind(hrIFNy_per_tumour_TLS, hrIFNy)
      rownames(hrIFNy_per_tumour_TLS)[[nrow(hrIFNy_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
    }
    
  }
  
  rownames(unstim_per_tumour_TLS) <- paste0("unstim TLS ", rownames(unstim_per_tumour_TLS))
  colnames(unstim_per_tumour_TLS) <- colnames(TLS_023_014_data)
  rownames(aPD1_per_tumour_TLS) <- paste0("aPD1 TLS ", rownames(aPD1_per_tumour_TLS))
  colnames(aPD1_per_tumour_TLS) <- colnames(TLS_023_014_data)
  rownames(aPD1_aIFNR_per_tumour_TLS) <- paste0("aPD1+aIFNyR1 TLS ", rownames(aPD1_aIFNR_per_tumour_TLS))
  colnames(aPD1_aIFNR_per_tumour_TLS) <- colnames(TLS_023_014_data)
  rownames(hrIFNy_per_tumour_TLS) <- paste0("hrIFNy TLS ", rownames(hrIFNy_per_tumour_TLS))
  colnames(hrIFNy_per_tumour_TLS) <- colnames(TLS_023_014_data)
  
  #for intermediate
  int_per_tumour_data <- list(int_023_013_data,int_023_014_data, int_035_013_data, int_035_015_data, int_052_010_data, int_052_014_data, int_066_015_data,
                              int_080_data, int_083_012_data, int_083_015_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data) #missing LU066 (TLS012)
  names(int_per_tumour_data) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014", "LU066_TLS015",
                       "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
  
  unstim_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))
  
  for(i in seq_along(int_per_tumour_data)){
    data <- int_per_tumour_data[[i]]
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
      unstim_per_tumour_int <- rbind(unstim_per_tumour_int, unstim)
      rownames(unstim_per_tumour_int)[[nrow(unstim_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_per_tumour_int <- rbind(aPD1_per_tumour_int, aPD1)
      rownames(aPD1_per_tumour_int)[[nrow(aPD1_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_per_tumour_int <- rbind(aPD1_aIFNR_per_tumour_int, aPD1_aIFNR)
      rownames(aPD1_aIFNR_per_tumour_int)[[nrow(aPD1_aIFNR_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_per_tumour_int <- rbind(hrIFNy_per_tumour_int, hrIFNy)
      rownames(hrIFNy_per_tumour_int)[[nrow(hrIFNy_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
    }
    
  }
  
  rownames(unstim_per_tumour_int) <- paste0("unstim intermediate ", rownames(unstim_per_tumour_int))
  colnames(unstim_per_tumour_int) <- colnames(TLS_023_014_data)
  rownames(aPD1_per_tumour_int) <- paste0("aPD1 intermediate ", rownames(aPD1_per_tumour_int))
  colnames(aPD1_per_tumour_int) <- colnames(TLS_023_014_data)
  rownames(aPD1_aIFNR_per_tumour_int) <- paste0("aPD1+aIFNyR1 intermediate ", rownames(aPD1_aIFNR_per_tumour_int))
  colnames(aPD1_aIFNR_per_tumour_int) <- colnames(TLS_023_014_data)
  rownames(hrIFNy_per_tumour_int) <- paste0("hrIFNy intermediate ", rownames(hrIFNy_per_tumour_int))
  colnames(hrIFNy_per_tumour_int) <- colnames(TLS_023_014_data)
  
  #for tumour
  tumour_per_tumour_data <- list(tumour_023_013_data,tumour_023_014_data, tumour_035_013_data, tumour_035_015_data, tumour_052_010_data, tumour_052_014_data,tumour_066_012_data,
                                 tumour_066_015_data, tumour_080_data, tumour_083_012_data, tumour_083_015_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data)
  names(tumour_per_tumour_data) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012", "LU066_TLS015",
                          "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
  
  unstim_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  
  for(i in seq_along(tumour_per_tumour_data)){
    data <- tumour_per_tumour_data[[i]]
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
      unstim_per_tumour_tumour <- rbind(unstim_per_tumour_tumour, unstim)
      rownames(unstim_per_tumour_tumour)[[nrow(unstim_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_per_tumour_tumour <- rbind(aPD1_per_tumour_tumour, aPD1)
      rownames(aPD1_per_tumour_tumour)[[nrow(aPD1_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_per_tumour_tumour <- rbind(aPD1_aIFNR_per_tumour_tumour, aPD1_aIFNR)
      rownames(aPD1_aIFNR_per_tumour_tumour)[[nrow(aPD1_aIFNR_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_per_tumour_tumour <- rbind(hrIFNy_per_tumour_tumour, hrIFNy)
      rownames(hrIFNy_per_tumour_tumour)[[nrow(hrIFNy_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
    }
    
  }
  
  rownames(unstim_per_tumour_tumour) <- paste0("unstim tumour bed ", rownames(unstim_per_tumour_tumour))
  colnames(unstim_per_tumour_tumour) <- colnames(TLS_023_014_data)
  rownames(aPD1_per_tumour_tumour) <- paste0("aPD1 tumour bed ", rownames(aPD1_per_tumour_tumour))
  colnames(aPD1_per_tumour_tumour) <- colnames(TLS_023_014_data)
  rownames(aPD1_aIFNR_per_tumour_tumour) <- paste0("aPD1+aIFNyR1 tumour bed ", rownames(aPD1_aIFNR_per_tumour_tumour))
  colnames(aPD1_aIFNR_per_tumour_tumour) <- colnames(TLS_023_014_data)
  rownames(hrIFNy_per_tumour_tumour) <- paste0("hrIFNy tumour bed ", rownames(hrIFNy_per_tumour_tumour))
  colnames(hrIFNy_per_tumour_tumour) <- colnames(TLS_023_014_data)
  
  #combine in a heatmap
  heatmap <- rbind(unstim_per_tumour_TLS,unstim_per_tumour_int,unstim_per_tumour_tumour,aPD1_per_tumour_TLS,aPD1_per_tumour_int,aPD1_per_tumour_tumour,
                   aPD1_aIFNR_per_tumour_TLS,aPD1_aIFNR_per_tumour_int,aPD1_aIFNR_per_tumour_tumour,hrIFNy_per_tumour_TLS,hrIFNy_per_tumour_int,hrIFNy_per_tumour_tumour)
  
  
  
  #z-scoring
  colSD <- apply(heatmap, 2, sd)
  colMean <- colMeans(heatmap)
  z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
  z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
  z_scores <- t(z_scores)
  
  #annotation (I gave up on this...too many samples..)
  #annotation <- data.frame(segment = c(rep("TLS",8), rep("intermediate",12), rep("tumour bed", 15)), condition = rep(c("unstim","aPD1","aPD1+aIFNyR1","hrIFNy"), each = 3))
  #rownames(annotation) <- colnames(z_scores)
  
  #heatmap
  
  rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y"
  rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"
  
  unclust <- pheatmap(z_scores,cluster_cols = F)
  save_pheatmap(unclust, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved_unclust.pdf")
  save_pheatmap(unclust, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved_unclust.png")
  
  clust <- pheatmap(z_scores,cluster_cols = T)
  save_pheatmap(unclust, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved_clust.pdf")
  save_pheatmap(unclust, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved_clust.png")
  
#---------------
#find fragments with high z-scores
#since this just concerns a dataframe of all fragments, paste all original data together:
lu023_TLS013_copy <- lu023_TLS013
rownames(lu023_TLS013_copy) <- gsub("LU023", "LU023_TLS013", rownames(lu023_TLS013_copy))

lu023_TLS014_copy <- lu023_TLS014
rownames(lu023_TLS014_copy) <- gsub("LU023", "LU023_TLS014", rownames(lu023_TLS014_copy))

lu035_TLS013_copy <- lu035_TLS013
rownames(lu035_TLS013_copy) <- gsub("TLS013 LU035", "LU035_TLS013", rownames(lu035_TLS013_copy))

lu035_TLS015_copy <- lu035_TLS015
rownames(lu035_TLS015_copy) <- gsub("TLS015 LU035", "LU035_TLS015", rownames(lu035_TLS015_copy))

lu052_TLS010_copy <- lu052_TLS010
rownames(lu052_TLS010_copy) <- gsub("TLS010 LU052", "LU052_TLS010", rownames(lu052_TLS010_copy))

lu052_TLS014_copy <- lu052_TLS014
rownames(lu052_TLS014_copy) <- gsub("TLS014 LU052", "LU052_TLS014", rownames(lu052_TLS014_copy))

lu066_TLS012_copy <- lu066_TLS012
rownames(lu066_TLS012_copy) <- gsub("LU066", "LU066_TLS012", rownames(lu066_TLS012_copy))

lu066_TLS015_copy <- lu066_TLS015
rownames(lu066_TLS015_copy) <- gsub("LU066", "LU066_TLS015", rownames(lu066_TLS015_copy))

lu066_TLS012_copy <- lu066_TLS012
rownames(lu066_TLS012_copy) <- gsub("LU066", "LU066_TLS012", rownames(lu066_TLS012_copy))

lu083_TLS012_copy <- lu083_TLS012
rownames(lu083_TLS012_copy) <- gsub("LU083", "LU083_TLS012", rownames(lu083_TLS012_copy))

lu083_TLS015_copy <- lu083_TLS015
rownames(lu083_TLS015_copy) <- gsub("LU083", "LU083_TLS015", rownames(lu083_TLS015_copy))


full_df <- rbind(lu023_TLS013_copy, lu023_TLS014_copy, lu035_TLS013_copy, lu035_TLS015_copy, lu052_TLS010_copy, lu052_TLS014_copy, lu066_TLS012_copy, lu066_TLS015_copy,
                 lu080, lu083_TLS012_copy, lu083_TLS015_copy, lu086_1, lu086_2, lu088, lu091)

full_df <- full_df[!grepl("average|RN440|hrIFNy_high", rownames(full_df)),]

#for IFNy-scoring:
full_df_IFNy <- full_df[!grepl("hrIFNy", rownames(full_df)),]


#z-scoring
colSD <- apply(full_df, 2, sd)
colMean <- colMeans(full_df)
outliers <- sweep(full_df, 2, colMean, "-") #subtract means per column
outliers <- sweep(outliers, 2, colSD, "/") #divide by standard deviation per column
outliers <- t(outliers)

outliers[outliers > 3 | outliers < -3] <- NA #filter out tumours with individual high z-scores

outlier_fragment <- data.frame(matrix(nrow = 25, ncol = 0))

for(i in 1:ncol(outliers)){
  if(any(is.na(outliers[,i]))){
    outlier_fragment <- cbind(outlier_fragment, outliers[,i])
    colnames(outlier_fragment)[[ncol(outlier_fragment)]] <- colnames(outliers)[i]
  } else {
    next
  }
}

#save result
write.csv(outlier_fragment, "~/Documents/Analyses/TLS/final_all_fragments/outliers.csv")


#calculate IFNy scores
colSD <- apply(full_df_IFNy, 2, sd)
colMean <- colMeans(full_df_IFNy)
ifny <- sweep(full_df_IFNy, 2, colMean, "-") #subtract means per column
ifny <- sweep(ifny, 2, colSD, "/") #divide by standard deviation per column
ifny <- t(ifny)
ifny <- ifny[c("CXCL9", "CXCL10", "CXCL11", "IFN-γ"),]
infy_score <- colSums(ifny, na.rm = T)

#redo heatmaps based on outlier information

#Get separate Legendplex scores per segment
TLS_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% TLS_023_013,]
TLS_023_013_data <- check_outliers_fragment(TLS_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")
int_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% int_023_013,]
int_023_013_data <- check_outliers_fragment(int_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")
tumour_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% tumour_023_013,]
tumour_023_013_data <- check_outliers_fragment(tumour_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")

TLS_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% TLS_023_014,]
TLS_023_014_data <- check_outliers_fragment(TLS_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")
int_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% int_023_014,]
int_023_014_data <- check_outliers_fragment(int_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")
tumour_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% tumour_023_014,]
tumour_023_014_data <- check_outliers_fragment(tumour_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")

TLS_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% TLS_035_013,]
TLS_035_013_data <- check_outliers_fragment(TLS_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")
int_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% int_035_013,]
int_035_013_data <- check_outliers_fragment(int_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")
tumour_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% tumour_035_013,]
tumour_035_013_data <- check_outliers_fragment(tumour_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")

TLS_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% TLS_035_015,]
TLS_035_015_data <- check_outliers_fragment(TLS_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")
int_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% int_035_015,]
int_035_015_data <- check_outliers_fragment(int_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")
tumour_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% tumour_035_015,]
tumour_035_015_data <- check_outliers_fragment(tumour_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")

#take out hrIFNy_high fragments for this sample
TLS_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% TLS_052_010,]
TLS_052_010_data <- check_outliers_fragment(TLS_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
TLS_052_010_data <- TLS_052_010_data[!grepl("hrIFNy_high", rownames(TLS_052_010_data)),]
int_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% int_052_010,]
int_052_010_data <- check_outliers_fragment(int_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
int_052_010_data <- int_052_010_data[!grepl("hrIFNy_high", rownames(int_052_010_data)),]
tumour_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% tumour_052_010,]
tumour_052_010_data <- check_outliers_fragment(tumour_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
tumour_052_010_data <- tumour_052_010_data[!grepl("hrIFNy_high", rownames(tumour_052_010_data)),]

TLS_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% TLS_052_014,]
TLS_052_014_data <- check_outliers_fragment(TLS_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")
int_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% int_052_014,]
int_052_014_data <- check_outliers_fragment(int_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")
tumour_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% tumour_052_014,]
tumour_052_014_data <- check_outliers_fragment(tumour_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")

#This tumour has no intermediate fragments. Therefore, this condition is skipped.
TLS_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% TLS_066_012,]
TLS_066_012_data <- check_outliers_fragment(TLS_066_012_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS012")
tumour_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% tumour_066_012,]
tumour_066_012_data <- check_outliers_fragment(tumour_066_012_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS012")

TLS_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% TLS_066_015,]
TLS_066_015_data <- check_outliers_fragment(TLS_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")
int_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% int_066_015,]
int_066_015_data <- check_outliers_fragment(int_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")
tumour_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% tumour_066_015,]
tumour_066_015_data <- check_outliers_fragment(tumour_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")

#LU080 has no TLS so we ignore this condition
int_080_data <- lu080[rownames(lu080) %in% int_080,]
int_080_data <- check_outliers_fragment(int_080_data, outliers = outlier_fragment, tumour = "LU080", experiment = NA)
#one cytokine (CXCL9) is NA/outlier for all fragments, set to NA:
int_080_data[int_080_data == "NaN"] <- NA

tumour_080_data <- lu080[rownames(lu080) %in% tumour_080,]
tumour_080_data <- check_outliers_fragment(tumour_080_data, outliers = outlier_fragment, tumour = "LU080", experiment = NA)

TLS_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% TLS_083_012,]
TLS_083_012_data <- check_outliers_fragment(TLS_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")
int_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% int_083_012,]
int_083_012_data <- check_outliers_fragment(int_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")
tumour_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% tumour_083_012,]
tumour_083_012_data <- check_outliers_fragment(tumour_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")

TLS_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% TLS_083_015,]
TLS_083_015_data <- check_outliers_fragment(TLS_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")
int_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% int_083_015,]
int_083_015_data <- check_outliers_fragment(int_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")
tumour_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% tumour_083_015,]
tumour_083_015_data <- check_outliers_fragment(tumour_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")

#the only TLS of LU086-1 were in RN440-15/RN440-16, which will not be considered. Therefore, we ignore this condition.
int_086_1_data <- lu086_1[rownames(lu086_1) %in% int_086_1,]
int_086_1_data <- check_outliers_fragment(int_086_1_data, outliers = outlier_fragment, tumour = "LU086-1", experiment = NA)
tumour_086_1_data <- lu086_1[rownames(lu086_1) %in% tumour_086_1,]
tumour_086_1_data <- check_outliers_fragment(tumour_086_1_data, outliers = outlier_fragment, tumour = "LU086-1", experiment = NA)

#the only TLS of LU086-2 was in RN440-15, which was not taken along in the Legendplex. Therefore, we ignore this condition.
int_086_2_data <- lu086_2[rownames(lu086_2) %in% int_086_2,]
int_086_2_data <- check_outliers_fragment(int_086_2_data, outliers = outlier_fragment, tumour = "LU086-2", experiment = NA)
tumour_086_2_data <- lu086_2[rownames(lu086_2) %in% tumour_086_2,]
tumour_086_2_data <- check_outliers_fragment(tumour_086_2_data, outliers = outlier_fragment, tumour = "LU086-2", experiment = NA)

TLS_088_data <- lu088[rownames(lu088) %in% TLS_088,]
TLS_088_data <- check_outliers_fragment(TLS_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)
int_088_data <- lu088[rownames(lu088) %in% int_088,]
int_088_data <- check_outliers_fragment(int_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)
tumour_088_data <- lu088[rownames(lu088) %in% tumour_088,]
tumour_088_data <- check_outliers_fragment(tumour_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)

TLS_091_data <- lu091[rownames(lu091) %in% TLS_091,]
TLS_091_data <- check_outliers_fragment(TLS_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)
int_091_data <- lu091[rownames(lu091) %in% int_091,]
int_091_data <- check_outliers_fragment(int_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)
tumour_091_data <- lu091[rownames(lu091) %in% tumour_091,]
tumour_091_data <- check_outliers_fragment(tumour_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)


#heatmap per tumour
generate_heatmap <- function(TLS_data, int_data, tumour_data, tumour_name, experiment){
  z_score <- rbind(TLS_data,int_data,tumour_data) #collect all outlier-removed data per tumour
  annotation <- data.frame(segment = c(rep("TLS", nrow(TLS_data)), rep("intermediate", nrow(int_data)), rep("tumour bed", nrow(tumour_data)))) #make the annotation based on the amount of fragments in each data frame
  rownames(annotation) <- rownames(z_score)
  
  #z-scoring
  colSD <- apply(z_score, 2, sd, na.rm = T)
  colSD[is.na(colSD)] <- 0 #happens when only one fragment has a value
  colMean <- colMeans(z_score, na.rm = T)
  z_score <- sweep(z_score, 2, colMean, "-") #subtract means per column
  z_score<- sweep(z_score, 2, colSD, "/") #divide by standard deviation per column
  z_score[z_score == "NaN"] <- 0 #happens when SD is zero
  z_score<- t(z_score)
  
  #one entry (LU080, CXCL9) has no values anywhere (except one fragment, but it makes no sense to calculate distances on one fragment): remove this mediator
  na_count <- rowSums(is.na(z_score)) #count the number of NA values in each row (mediator)
  z_score <- z_score[na_count < (ncol(z_score)-1),] #only keep the mediator if the amount of NAs is lower than the amount of samples -1 (2 or more non-NA values should exist)
  
  #make a heatmap
  heatmap<- pheatmap(z_score, annotation_col = annotation, cluster_cols = F)
  
  #store the result
  pdf(paste0("~/Documents/Analyses/TLS/final_all_fragments/heatmap_", tumour_name, "_", experiment, ".pdf"),width = 15)
  heatmap
  dev.off()
  
  save_pheatmap(heatmap, paste0("~/Documents/Analyses/TLS/final_all_fragments/heatmap_", tumour_name, "_", experiment, ".png"))
}

generate_heatmap(TLS_023_013_data,int_023_013_data,tumour_023_013_data, "LU023", "TLS013")
dev.off()
generate_heatmap(TLS_023_014_data,int_023_014_data,tumour_023_014_data, "LU023", "TLS014")
dev.off()
generate_heatmap(TLS_035_013_data,int_035_013_data,tumour_035_013_data, "LU035", "TLS013")
dev.off()
generate_heatmap(TLS_035_015_data,int_035_015_data,tumour_035_015_data, "LU035", "TLS015")
dev.off()
generate_heatmap(TLS_052_010_data,int_052_010_data,tumour_052_010_data, "LU052", "TLS010")
dev.off()
generate_heatmap(TLS_052_014_data,int_052_014_data,tumour_052_014_data, "LU052", "TLS014")
dev.off()
generate_heatmap(TLS_066_012_data,data.frame(),tumour_066_012_data, "LU066", "TLS012") #this entry has no intermediate fragments
dev.off()
generate_heatmap(TLS_066_015_data,int_066_015_data,tumour_066_015_data, "LU066", "TLS015")
dev.off()
generate_heatmap(TLS_data = data.frame(), int_data = int_080_data, tumour_data = tumour_080_data, tumour =  "LU080", experiment = "TLS011") #this entry has no TLS or CXCL9 (all outlier)
dev.off()
generate_heatmap(TLS_083_012_data,int_083_012_data,tumour_083_012_data, "LU083", "TLS012")
dev.off()
generate_heatmap(TLS_083_015_data,int_083_015_data,tumour_083_015_data, "LU083", "TLS015")
dev.off()
generate_heatmap(TLS_data = data.frame(),int_086_1_data,tumour_086_1_data, "LU086-1", "TLS011") #this entry has no TLS
dev.off()
generate_heatmap(TLS_data = data.frame(),int_086_2_data,tumour_086_2_data, "LU086-2", "TLS011") #this entry has no TLS
dev.off()
generate_heatmap(TLS_088_data,int_088_data,tumour_088_data, "LU088", "TLS016")
dev.off()
generate_heatmap(TLS_091_data,int_091_data,tumour_091_data, "LU091", "TLS016")
dev.off()


#calculate averages for all tumours
TLS_023_013_data <- colMeans(TLS_023_013_data, na.rm = T) #average score per cytokine across all conditions (ignore NA values)
int_023_013_data <- colMeans(int_023_013_data, na.rm = T)
tumour_023_013_data <- colMeans(tumour_023_013_data, na.rm = T)
TLS_023_014_data <- colMeans(TLS_023_014_data, na.rm = T) #average score per cytokine across all conditions
int_023_014_data <- colMeans(int_023_014_data, na.rm = T)
tumour_023_014_data <- colMeans(tumour_023_014_data, na.rm = T)
TLS_035_013_data <- colMeans(TLS_035_013_data, na.rm = T) #average score per cytokine across all conditions
int_035_013_data <- colMeans(int_035_013_data, na.rm = T)
tumour_035_013_data <- colMeans(tumour_035_013_data, na.rm = T)
TLS_035_015_data <- colMeans(TLS_035_015_data, na.rm = T) #average score per cytokine across all conditions
int_035_015_data <- colMeans(int_035_015_data, na.rm = T)
tumour_035_015_data <- colMeans(tumour_035_015_data, na.rm = T)
TLS_052_010_data <- colMeans(TLS_052_010_data, na.rm = T) #average score per cytokine across all conditions
int_052_010_data <- colMeans(int_052_010_data, na.rm = T)
tumour_052_010_data <- colMeans(tumour_052_010_data, na.rm = T)
TLS_052_014_data <- colMeans(TLS_052_014_data, na.rm = T) #average score per cytokine across all conditions
int_052_014_data <- colMeans(int_052_014_data, na.rm = T)
tumour_052_014_data <- colMeans(tumour_052_014_data, na.rm = T)
TLS_066_012_data <- colMeans(TLS_066_012_data, na.rm = T) #average score per cytokine across all conditions
tumour_066_012_data <- colMeans(tumour_066_012_data, na.rm = T)
TLS_066_015_data <- colMeans(TLS_066_015_data, na.rm = T) #average score per cytokine across all conditions
int_066_015_data <- colMeans(int_066_015_data, na.rm = T)
tumour_066_015_data <- colMeans(tumour_066_015_data, na.rm = T)
int_080_data <- colMeans(int_080_data, na.rm = T)
tumour_080_data <- colMeans(tumour_080_data, na.rm = T)
TLS_083_012_data <- colMeans(TLS_083_012_data, na.rm = T) #average score per cytokine across all conditions
int_083_012_data <- colMeans(int_083_012_data, na.rm = T)
tumour_083_012_data <- colMeans(tumour_083_012_data, na.rm = T)
TLS_083_015_data <- colMeans(TLS_083_015_data, na.rm = T) #average score per cytokine across all conditions
int_083_015_data <- colMeans(int_083_015_data, na.rm = T)
tumour_083_015_data <- colMeans(tumour_083_015_data, na.rm = T)
int_086_1_data <- colMeans(int_086_1_data, na.rm = T)
tumour_086_1_data <- colMeans(tumour_086_1_data, na.rm = T)
int_086_2_data <- colMeans(int_086_2_data, na.rm = T)
tumour_086_2_data <- colMeans(tumour_086_2_data, na.rm = T)
TLS_088_data <- colMeans(TLS_088_data, na.rm = T) #average score per cytokine across all conditions
int_088_data <- colMeans(int_088_data, na.rm = T)
tumour_088_data <- colMeans(tumour_088_data, na.rm = T)
TLS_091_data <- colMeans(TLS_091_data, na.rm = T) #average score per cytokine across all conditions
int_091_data <- colMeans(int_091_data, na.rm = T)
tumour_091_data <- colMeans(tumour_091_data, na.rm = T)

#generate heatmap per tumour and segment


#Average across all tumours
TLS_scores <- rbind(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                    TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data) #missing LU080, LU086-1, LU086-2
TLS_scores_mean <- colMeans(TLS_scores, na.rm = T)
rownames(TLS_scores) <- gsub("TLS_", "TLS LU", rownames(TLS_scores))
rownames(TLS_scores) <- gsub("_data", "", rownames(TLS_scores))
rownames(TLS_scores) <- gsub("_0", " TLS0", rownames(TLS_scores))


int_scores <- rbind(int_023_013_data,int_023_014_data, int_035_013_data, int_035_015_data, int_052_010_data, int_052_014_data, int_066_015_data,
                    int_080_data, int_083_012_data, int_083_015_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data) #missing LU066 (TLS012)
int_scores_mean <- colMeans(int_scores, na.rm = T)
rownames(int_scores) <- gsub("int_", "Intermediate LU", rownames(int_scores))
rownames(int_scores) <- gsub("_data", "", rownames(int_scores))
rownames(int_scores) <- gsub("_0", " TLS0", rownames(int_scores))

tumour_scores <- rbind(tumour_023_013_data,tumour_023_014_data, tumour_035_013_data, tumour_035_015_data, tumour_052_010_data, tumour_052_014_data, tumour_066_012_data,tumour_066_015_data,
                       tumour_080_data, tumour_083_012_data, tumour_083_015_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data)
tumour_scores_mean <- colMeans(tumour_scores, na.rm = T)
rownames(tumour_scores) <- gsub("tumour_", "Tumour_bed LU", rownames(tumour_scores))
rownames(tumour_scores) <- gsub("_data", "", rownames(tumour_scores))
rownames(tumour_scores) <- gsub("_0", " TLS0", rownames(tumour_scores))

heatmap_scores_tumour <- rbind(TLS_scores, int_scores, tumour_scores)
colSD <- apply(heatmap_scores_tumour, 2, sd, na.rm = T)
colMean <- colMeans(heatmap_scores_tumour, na.rm = T)
z_scores_tumour <- sweep(heatmap_scores_tumour, 2, colMean, "-") #subtract means per column
z_scores_tumour <- sweep(z_scores_tumour, 2, colSD, "/") #divide by standard deviation per column
z_scores_tumour <- t(z_scores_tumour)
rownames(z_scores_tumour)[rownames(z_scores_tumour) == "IFN-γ"] <- "IFN-y"
rownames(z_scores_tumour)[rownames(z_scores_tumour) == "TNF-α"] <- "TNF-a"

#add response score
annotation_tumour <- data.frame(PDTF_responder = rep(c("PDTF_NR","PDTF_R","PDTF_NR","None","PDTF_R","None","PDTF_NR","PDTF_NR","PDTF_NR","PDTF_R",
                                                        "PDTF_R","PDTF_NR","PDTF_NR","PDTF_NR","PDTF_NR"), times = 3))
#remove for segments that do not exist (LU066_012 int,LU080 TLS, etc.)
annotation_tumour <- as.data.frame(annotation_tumour[!(rownames(annotation_tumour) %in% c(22,12,13,9)),])
colnames(annotation_tumour) <- "PDTF_score"
rownames(annotation_tumour) <- colnames(z_scores_tumour)
annotation_tumour$PDTF_score <- factor(annotation_tumour$PDTF_score, levels = c("PDTF_NR", "PDTF_R", "None"))

ann_colors = list(
  PDTF_score = c(PDTF_NR = "pink",PDTF_R= "lightgreen",None= "grey"))

heatmap <- suppressWarnings(pheatmap(z_scores_tumour, annotation_col = annotation_tumour, annotation_colors = ann_colors))  
heatmap
save_pheatmap(heatmap, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved.pdf")
save_pheatmap(heatmap, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved.png")

heatmap_unclust <- suppressWarnings(pheatmap(z_scores_tumour, annotation_col = annotation_tumour, annotation_colors = ann_colors, cluster_cols = F))  
save_pheatmap(heatmap, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved_unclust.pdf")
save_pheatmap(heatmap, "~/Documents/Analyses/TLS/final_all_fragments/segment_tumours_normalised_outlierremoved_unclust.png")




#Generate a heatmap of z-scores per cytokine
#z-score is (value - average)/standard deviation. (is this the delta or is the average per condition?)
heatmap_scores <- rbind(TLS_scores_mean, int_scores_mean, tumour_scores_mean)
rownames(heatmap_scores) <- c("TLS", "intermediate", "tumour bed")
colSD <- apply(heatmap_scores, 2, sd)
colMean <- colMeans(heatmap_scores)
z_scores <- sweep(heatmap_scores, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)
rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y"
rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"

#heatmap
heatmap <- suppressWarnings(pheatmap(z_scores))  
heatmap

save_pheatmap(heatmap, "~/Documents/Analyses/TLS/final_all_fragments/segment_normalised_outlierremoved.pdf")




#split per condition
TLS_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% TLS_023_013,]
TLS_023_013_data <- check_outliers_fragment(TLS_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")
int_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% int_023_013,]
int_023_013_data <- check_outliers_fragment(int_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")
tumour_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% tumour_023_013,]
tumour_023_013_data <- check_outliers_fragment(tumour_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")

TLS_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% TLS_023_014,]
TLS_023_014_data <- check_outliers_fragment(TLS_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")
int_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% int_023_014,]
int_023_014_data <- check_outliers_fragment(int_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")
tumour_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% tumour_023_014,]
tumour_023_014_data <- check_outliers_fragment(tumour_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")

TLS_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% TLS_035_013,]
TLS_035_013_data <- check_outliers_fragment(TLS_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")
int_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% int_035_013,]
int_035_013_data <- check_outliers_fragment(int_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")
tumour_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% tumour_035_013,]
tumour_035_013_data <- check_outliers_fragment(tumour_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")

TLS_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% TLS_035_015,]
TLS_035_015_data <- check_outliers_fragment(TLS_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")
int_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% int_035_015,]
int_035_015_data <- check_outliers_fragment(int_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")
tumour_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% tumour_035_015,]
tumour_035_015_data <- check_outliers_fragment(tumour_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")

TLS_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% TLS_052_010,]
TLS_052_010_data <- check_outliers_fragment(TLS_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
int_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% int_052_010,]
int_052_010_data <- check_outliers_fragment(int_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
tumour_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% tumour_052_010,]
tumour_052_010_data <- check_outliers_fragment(tumour_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
#LU052_TLS010_tumour_bed has extremely high CXCL9/CXCL10, what does the heatmap look like if I take the IFN_high fragments out?
TLS_052_010_data <- TLS_052_010_data[!grepl("hrIFNy_high", rownames(TLS_052_010_data)),]
int_052_010_data <- int_052_010_data[!grepl("hrIFNy_high", rownames(int_052_010_data)),]
tumour_052_010_data <- tumour_052_010_data[!grepl("hrIFNy_high", rownames(tumour_052_010_data)),]

TLS_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% TLS_052_014,]
TLS_052_014_data <- check_outliers_fragment(TLS_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")
int_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% int_052_014,]
int_052_014_data <- check_outliers_fragment(int_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")
tumour_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% tumour_052_014,]
tumour_052_014_data <- check_outliers_fragment(tumour_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")

#This tumour has no intermediate fragments. Therefore, this condition is skipped.
TLS_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% TLS_066_012,]
TLS_066_012_data <- check_outliers_fragment(TLS_066_012_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS012")
tumour_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% tumour_066_012,]
tumour_066_012_data <- check_outliers_fragment(tumour_066_012_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS012")

TLS_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% TLS_066_015,]
TLS_066_015_data <- check_outliers_fragment(TLS_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")
int_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% int_066_015,]
int_066_015_data <- check_outliers_fragment(int_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")
tumour_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% tumour_066_015,]
tumour_066_015_data <- check_outliers_fragment(tumour_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")

#LU080 has no TLS so we ignore this condition
int_080_data <- lu080[rownames(lu080) %in% int_080,]
int_080_data <- check_outliers_fragment(int_080_data, outliers = outlier_fragment, tumour = "LU080", experiment = NA)
tumour_080_data <- lu080[rownames(lu080) %in% tumour_080,]
tumour_080_data <- check_outliers_fragment(tumour_080_data, outliers = outlier_fragment, tumour = "LU080", experiment = NA)

TLS_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% TLS_083_012,]
TLS_083_012_data <- check_outliers_fragment(TLS_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")
int_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% int_083_012,]
int_083_012_data <- check_outliers_fragment(int_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")
tumour_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% tumour_083_012,]
tumour_083_012_data <- check_outliers_fragment(tumour_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")

TLS_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% TLS_083_015,]
TLS_083_015_data <- check_outliers_fragment(TLS_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")
int_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% int_083_015,]
int_083_015_data <- check_outliers_fragment(int_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")
tumour_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% tumour_083_015,]
tumour_083_015_data <- check_outliers_fragment(tumour_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")

#the only TLS of LU086-1 were in RN440-15/RN440-16, which will not be considered. Therefore, we ignore this condition.
int_086_1_data <- lu086_1[rownames(lu086_1) %in% int_086_1,]
int_086_1_data <- check_outliers_fragment(int_086_1_data, outliers = outlier_fragment, tumour = "LU086-1", experiment = NA)
tumour_086_1_data <- lu086_1[rownames(lu086_1) %in% tumour_086_1,]
tumour_086_1_data <- check_outliers_fragment(tumour_086_1_data, outliers = outlier_fragment, tumour = "LU086-1", experiment = NA)

#the only TLS of LU086-2 was in RN440-15, which was not taken along in the Legendplex. Therefore, we ignore this condition.
int_086_2_data <- lu086_2[rownames(lu086_2) %in% int_086_2,]
int_086_2_data <- check_outliers_fragment(int_086_2_data, outliers = outlier_fragment, tumour = "LU086-2", experiment = NA)
tumour_086_2_data <- lu086_2[rownames(lu086_2) %in% tumour_086_2,]
tumour_086_2_data <- check_outliers_fragment(tumour_086_2_data, outliers = outlier_fragment, tumour = "LU086-2", experiment = NA)

TLS_088_data <- lu088[rownames(lu088) %in% TLS_088,]
TLS_088_data <- check_outliers_fragment(TLS_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)
int_088_data <- lu088[rownames(lu088) %in% int_088,]
int_088_data <- check_outliers_fragment(int_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)
tumour_088_data <- lu088[rownames(lu088) %in% tumour_088,]
tumour_088_data <- check_outliers_fragment(tumour_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)

TLS_091_data <- lu091[rownames(lu091) %in% TLS_091,]
TLS_091_data <- check_outliers_fragment(TLS_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)
int_091_data <- lu091[rownames(lu091) %in% int_091,]
int_091_data <- check_outliers_fragment(int_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)
tumour_091_data <- lu091[rownames(lu091) %in% tumour_091,]
tumour_091_data <- check_outliers_fragment(tumour_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)

#First for TLS
TLS_list <- list(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                 TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data) #missing LU080, LU086-1, LU086-2
names(TLS_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012","LU066_TLS015",
                     "LU083_TLS012", "LU083_TLS015", "LU088", "LU091")  

unstim_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(TLS_list)){
  data <- TLS_list[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_TLS <- rbind(unstim_average_TLS, unstim)
    rownames(unstim_average_TLS)[[nrow(unstim_average_TLS)]] <- names(TLS_list)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_TLS <- rbind(aPD1_average_TLS, aPD1)
    rownames(aPD1_average_TLS)[[nrow(aPD1_average_TLS)]] <- names(TLS_list)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_TLS <- rbind(aPD1_aIFNR_average_TLS, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_TLS)[[nrow(aPD1_aIFNR_average_TLS)]] <- names(TLS_list)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_TLS <- rbind(hrIFNy_average_TLS, hrIFNy)
    rownames(hrIFNy_average_TLS)[[nrow(hrIFNy_average_TLS)]] <- names(TLS_list)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_TLS) <- colnames(TLS_023_014_data)
unstim_average_TLS <- colMeans(unstim_average_TLS, na.rm = T)

colnames(aPD1_average_TLS) <- colnames(TLS_023_014_data)
aPD1_average_TLS <- colMeans(aPD1_average_TLS, na.rm = T)

colnames(aPD1_aIFNR_average_TLS) <- colnames(TLS_023_014_data)
aPD1_aIFNR_average_TLS <- colMeans(aPD1_aIFNR_average_TLS, na.rm = T)

colnames(hrIFNy_average_TLS) <- colnames(TLS_023_014_data)
hrIFNy_average_TLS <- colMeans(hrIFNy_average_TLS, na.rm = T)

#repeat for intermediate and tumour sections
int_list <- list(int_023_013_data,int_023_014_data, int_035_013_data, int_035_015_data, int_052_010_data, int_052_014_data, int_066_015_data,
                 int_080_data, int_083_012_data, int_083_015_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data) #missing LU066 (TLS012)
names(int_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014", "LU066_TLS015",
                     "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
unstim_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_int <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(int_list)){
  data <- int_list[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_int <- rbind(unstim_average_int, unstim)
    rownames(unstim_average_int)[[nrow(unstim_average_int)]] <- names(int_list)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_int <- rbind(aPD1_average_int, aPD1)
    rownames(aPD1_average_int)[[nrow(aPD1_average_int)]] <- names(int_list)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_int <- rbind(aPD1_aIFNR_average_int, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_int)[[nrow(aPD1_aIFNR_average_int)]] <- names(int_list)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_int <- rbind(hrIFNy_average_int, hrIFNy)
    rownames(hrIFNy_average_int)[[nrow(hrIFNy_average_int)]] <- names(int_list)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_int) <- colnames(TLS_023_014_data)
unstim_average_int <- colMeans(unstim_average_int, na.rm = T)

colnames(aPD1_average_int) <- colnames(TLS_023_014_data)
aPD1_average_int <- colMeans(aPD1_average_int, na.rm = T)

colnames(aPD1_aIFNR_average_int) <- colnames(TLS_023_014_data)
aPD1_aIFNR_average_int <- colMeans(aPD1_aIFNR_average_int, na.rm = T)

colnames(hrIFNy_average_int) <- colnames(TLS_023_014_data)
hrIFNy_average_int <- colMeans(hrIFNy_average_int, na.rm = T)

#and tumour
tumour_list <- list(tumour_023_013_data,tumour_023_014_data, tumour_035_013_data, tumour_035_015_data, tumour_052_010_data, tumour_052_014_data,tumour_066_012_data,
                    tumour_066_015_data, tumour_080_data, tumour_083_012_data, tumour_083_015_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data)
names(tumour_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012", "LU066_TLS015",
                        "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
unstim_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(tumour_list)){
  data <- tumour_list[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_tumour <- rbind(unstim_average_tumour, unstim)
    rownames(unstim_average_tumour)[[nrow(unstim_average_tumour)]] <- names(tumour_list)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_tumour <- rbind(aPD1_average_tumour, aPD1)
    rownames(aPD1_average_tumour)[[nrow(aPD1_average_tumour)]] <- names(tumour_list)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_tumour <- rbind(aPD1_aIFNR_average_tumour, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_tumour)[[nrow(aPD1_aIFNR_average_tumour)]] <- names(tumour_list)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_tumour <- rbind(hrIFNy_average_tumour, hrIFNy)
    rownames(hrIFNy_average_tumour)[[nrow(hrIFNy_average_tumour)]] <- names(tumour_list)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_tumour) <- colnames(TLS_023_014_data)
unstim_average_tumour <- colMeans(unstim_average_tumour, na.rm = T)

colnames(aPD1_average_tumour) <- colnames(TLS_023_014_data)
aPD1_average_tumour <- colMeans(aPD1_average_tumour, na.rm = T)

colnames(aPD1_aIFNR_average_tumour) <- colnames(TLS_023_014_data)
aPD1_aIFNR_average_tumour <- colMeans(aPD1_aIFNR_average_tumour, na.rm = T)

colnames(hrIFNy_average_tumour) <- colnames(TLS_023_014_data)
hrIFNy_average_tumour <- colMeans(hrIFNy_average_tumour, na.rm = T)




#combine in a heatmap
heatmap <- rbind(unstim_average_TLS,unstim_average_int,unstim_average_tumour,aPD1_average_TLS,aPD1_average_int,aPD1_average_tumour,
                 aPD1_aIFNR_average_TLS,aPD1_aIFNR_average_int,aPD1_aIFNR_average_tumour,hrIFNy_average_TLS,hrIFNy_average_int,hrIFNy_average_tumour)

rownames(heatmap) <- c("unstim TLS" ,"unstim intermediate","unstim tumour bed","aPD1 TLS","aPD1 intermediate","aPD1 tumour bed",
                       "aPD1+aIFNyR1 TLS","aPD1+aIFNyR1 intermediate","aPD1+aIFNyR1 tumour bed","hrIFNy TLS","hrIFNy intermediate","hrIFNy tumour bed")  

#z-scoring
colSD <- apply(heatmap, 2, sd)
colMean <- colMeans(heatmap)
z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y" #to make sure rownames can be seen
rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"

#annotation
annotation <- data.frame(segment = rep(c("TLS", "intermediate", "tumor bed"),4), condition = rep(c("unstim","aPD1","aPD1+aIFNyR1","hrIFNy"), each = 3))
rownames(annotation) <- colnames(z_scores)

#heatmap
correlations <- pheatmap(z_scores, annotation_col = annotation,cluster_cols = F)

correlations

#capture row order
row_order <- correlations$tree_row$order
values <- rownames(z_scores)[row_order]

#reorder z-scores based on row order
z_scores <- z_scores[order(match(rownames(z_scores), values)),]

#save z-scores
z_scores <- as.data.frame(z_scores)
write.xlsx(z_scores, "~/Documents/Analyses/TLS/final_all_fragments/z_scores_condition_segment_normalised_outlierremoved.xlsx", rowNames = T)

#-------------------------
#split up between R and NR (unstim-aPD1-combi)
TLS_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% TLS_023_013,]
TLS_023_013_data <- check_outliers_fragment(TLS_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")
int_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% int_023_013,]
int_023_013_data <- check_outliers_fragment(int_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")
tumour_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% tumour_023_013,]
tumour_023_013_data <- check_outliers_fragment(tumour_023_013_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS013")

TLS_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% TLS_023_014,]
TLS_023_014_data <- check_outliers_fragment(TLS_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")
int_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% int_023_014,]
int_023_014_data <- check_outliers_fragment(int_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")
tumour_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% tumour_023_014,]
tumour_023_014_data <- check_outliers_fragment(tumour_023_014_data, outliers = outlier_fragment, tumour = "LU023", experiment = "TLS014")

TLS_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% TLS_035_013,]
TLS_035_013_data <- check_outliers_fragment(TLS_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")
int_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% int_035_013,]
int_035_013_data <- check_outliers_fragment(int_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")
tumour_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% tumour_035_013,]
tumour_035_013_data <- check_outliers_fragment(tumour_035_013_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS013")

TLS_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% TLS_035_015,]
TLS_035_015_data <- check_outliers_fragment(TLS_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")
int_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% int_035_015,]
int_035_015_data <- check_outliers_fragment(int_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")
tumour_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% tumour_035_015,]
tumour_035_015_data <- check_outliers_fragment(tumour_035_015_data, outliers = outlier_fragment, tumour = "LU035", experiment = "TLS015")

TLS_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% TLS_052_010,]
TLS_052_010_data <- check_outliers_fragment(TLS_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
int_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% int_052_010,]
int_052_010_data <- check_outliers_fragment(int_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
tumour_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% tumour_052_010,]
tumour_052_010_data <- check_outliers_fragment(tumour_052_010_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS010")
#LU052_TLS010_tumour_bed has extremely high CXCL9/CXCL10, what does the heatmap look like if I take the IFN_high fragments out?
TLS_052_010_data <- TLS_052_010_data[!grepl("hrIFNy_high", rownames(TLS_052_010_data)),]
int_052_010_data <- int_052_010_data[!grepl("hrIFNy_high", rownames(int_052_010_data)),]
tumour_052_010_data <- tumour_052_010_data[!grepl("hrIFNy_high", rownames(tumour_052_010_data)),]

TLS_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% TLS_052_014,]
TLS_052_014_data <- check_outliers_fragment(TLS_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")
int_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% int_052_014,]
int_052_014_data <- check_outliers_fragment(int_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")
tumour_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% tumour_052_014,]
tumour_052_014_data <- check_outliers_fragment(tumour_052_014_data, outliers = outlier_fragment, tumour = "LU052", experiment = "TLS014")

#This tumour has no intermediate fragments. Therefore, this condition is skipped.
TLS_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% TLS_066_012,]
TLS_066_012_data <- check_outliers_fragment(TLS_066_012_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS012")
tumour_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% tumour_066_012,]
tumour_066_012_data <- check_outliers_fragment(tumour_066_012_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS012")

TLS_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% TLS_066_015,]
TLS_066_015_data <- check_outliers_fragment(TLS_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")
int_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% int_066_015,]
int_066_015_data <- check_outliers_fragment(int_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")
tumour_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% tumour_066_015,]
tumour_066_015_data <- check_outliers_fragment(tumour_066_015_data, outliers = outlier_fragment, tumour = "LU066", experiment = "TLS015")

#LU080 has no TLS so we ignore this condition
int_080_data <- lu080[rownames(lu080) %in% int_080,]
int_080_data <- check_outliers_fragment(int_080_data, outliers = outlier_fragment, tumour = "LU080", experiment = NA)
tumour_080_data <- lu080[rownames(lu080) %in% tumour_080,]
tumour_080_data <- check_outliers_fragment(tumour_080_data, outliers = outlier_fragment, tumour = "LU080", experiment = NA)

TLS_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% TLS_083_012,]
TLS_083_012_data <- check_outliers_fragment(TLS_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")
int_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% int_083_012,]
int_083_012_data <- check_outliers_fragment(int_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")
tumour_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% tumour_083_012,]
tumour_083_012_data <- check_outliers_fragment(tumour_083_012_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS012")

TLS_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% TLS_083_015,]
TLS_083_015_data <- check_outliers_fragment(TLS_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")
int_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% int_083_015,]
int_083_015_data <- check_outliers_fragment(int_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")
tumour_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% tumour_083_015,]
tumour_083_015_data <- check_outliers_fragment(tumour_083_015_data, outliers = outlier_fragment, tumour = "LU083", experiment = "TLS015")

#the only TLS of LU086-1 were in RN440-15/RN440-16, which will not be considered. Therefore, we ignore this condition.
int_086_1_data <- lu086_1[rownames(lu086_1) %in% int_086_1,]
int_086_1_data <- check_outliers_fragment(int_086_1_data, outliers = outlier_fragment, tumour = "LU086-1", experiment = NA)
tumour_086_1_data <- lu086_1[rownames(lu086_1) %in% tumour_086_1,]
tumour_086_1_data <- check_outliers_fragment(tumour_086_1_data, outliers = outlier_fragment, tumour = "LU086-1", experiment = NA)

#the only TLS of LU086-2 was in RN440-15, which was not taken along in the Legendplex. Therefore, we ignore this condition.
int_086_2_data <- lu086_2[rownames(lu086_2) %in% int_086_2,]
int_086_2_data <- check_outliers_fragment(int_086_2_data, outliers = outlier_fragment, tumour = "LU086-2", experiment = NA)
tumour_086_2_data <- lu086_2[rownames(lu086_2) %in% tumour_086_2,]
tumour_086_2_data <- check_outliers_fragment(tumour_086_2_data, outliers = outlier_fragment, tumour = "LU086-2", experiment = NA)

TLS_088_data <- lu088[rownames(lu088) %in% TLS_088,]
TLS_088_data <- check_outliers_fragment(TLS_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)
int_088_data <- lu088[rownames(lu088) %in% int_088,]
int_088_data <- check_outliers_fragment(int_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)
tumour_088_data <- lu088[rownames(lu088) %in% tumour_088,]
tumour_088_data <- check_outliers_fragment(tumour_088_data, outliers = outlier_fragment, tumour = "LU088", experiment = NA)

TLS_091_data <- lu091[rownames(lu091) %in% TLS_091,]
TLS_091_data <- check_outliers_fragment(TLS_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)
int_091_data <- lu091[rownames(lu091) %in% int_091,]
int_091_data <- check_outliers_fragment(int_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)
tumour_091_data <- lu091[rownames(lu091) %in% tumour_091,]
tumour_091_data <- check_outliers_fragment(tumour_091_data, outliers = outlier_fragment, tumour = "LU091", experiment = NA)

#First for TLS
TLS_list_R <- list(TLS_023_013_data,TLS_023_014_data, TLS_052_010_data, TLS_052_014_data,TLS_083_012_data, TLS_083_015_data)
names(TLS_list_R) <- c("LU023_TLS013","LU023_TLS014","LU052_TLS010","LU052_TLS014","LU083_TLS012", "LU083_TLS015")  

TLS_list_NR <- list(TLS_035_013_data, TLS_035_015_data, TLS_066_012_data,TLS_066_015_data, TLS_088_data, TLS_091_data)
names(TLS_list_NR) <- c("LU035_TLS013","LU035_TLS015","LU066_TLS012","LU066_TLS015","LU088", "LU091")  

unstim_average_TLS_R <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_TLS_R <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_TLS_R <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_TLS_R <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(TLS_list_R)){
  data <- TLS_list_R[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_TLS_R <- rbind(unstim_average_TLS_R, unstim)
    rownames(unstim_average_TLS_R)[[nrow(unstim_average_TLS_R)]] <- names(TLS_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_TLS_R <- rbind(aPD1_average_TLS_R, aPD1)
    rownames(aPD1_average_TLS_R)[[nrow(aPD1_average_TLS_R)]] <- names(TLS_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_TLS_R <- rbind(aPD1_aIFNR_average_TLS_R, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_TLS_R)[[nrow(aPD1_aIFNR_average_TLS_R)]] <- names(TLS_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_TLS_R <- rbind(hrIFNy_average_TLS_R, hrIFNy)
    rownames(hrIFNy_average_TLS_R)[[nrow(hrIFNy_average_TLS_R)]] <- names(TLS_list_R)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_TLS_R) <- colnames(TLS_023_014_data)
unstim_average_TLS_R <- colMeans(unstim_average_TLS_R, na.rm = T)

colnames(aPD1_average_TLS_R) <- colnames(TLS_023_014_data)
aPD1_average_TLS_R <- colMeans(aPD1_average_TLS_R, na.rm = T)

colnames(aPD1_aIFNR_average_TLS_R) <- colnames(TLS_023_014_data)
aPD1_aIFNR_average_TLS_R <- colMeans(aPD1_aIFNR_average_TLS_R, na.rm = T)

colnames(hrIFNy_average_TLS_R) <- colnames(TLS_023_014_data)
hrIFNy_average_TLS_R <- colMeans(hrIFNy_average_TLS_R, na.rm = T)

#repeat for non-responders
unstim_average_TLS_NR <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_TLS_NR <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_TLS_NR <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_TLS_NR <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(TLS_list_NR)){
  data <- TLS_list_NR[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_TLS_NR <- rbind(unstim_average_TLS_NR, unstim)
    rownames(unstim_average_TLS_NR)[[nrow(unstim_average_TLS_NR)]] <- names(TLS_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_TLS_NR <- rbind(aPD1_average_TLS_NR, aPD1)
    rownames(aPD1_average_TLS_NR)[[nrow(aPD1_average_TLS_NR)]] <- names(TLS_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_TLS_NR <- rbind(aPD1_aIFNR_average_TLS_NR, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_TLS_NR)[[nrow(aPD1_aIFNR_average_TLS_NR)]] <- names(TLS_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_TLS_NR <- rbind(hrIFNy_average_TLS_NR, hrIFNy)
    rownames(hrIFNy_average_TLS_NR)[[nrow(hrIFNy_average_TLS_NR)]] <- names(TLS_list_NR)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_TLS_NR) <- colnames(TLS_023_014_data)
unstim_average_TLS_NR <- colMeans(unstim_average_TLS_NR, na.rm = T)

colnames(aPD1_average_TLS_NR) <- colnames(TLS_023_014_data)
aPD1_average_TLS_NR <- colMeans(aPD1_average_TLS_NR, na.rm = T)

colnames(aPD1_aIFNR_average_TLS_NR) <- colnames(TLS_023_014_data)
aPD1_aIFNR_average_TLS_NR <- colMeans(aPD1_aIFNR_average_TLS_NR, na.rm = T)

colnames(hrIFNy_average_TLS_NR) <- colnames(TLS_023_014_data)
hrIFNy_average_TLS_NR <- colMeans(hrIFNy_average_TLS_NR, na.rm = T)



#repeat for intermediate and tumour sections
int_list_R <- list(int_023_013_data,int_023_014_data,int_052_010_data, int_052_014_data,int_083_012_data, int_083_015_data)
names(int_list_R) <- c("LU023_TLS013","LU023_TLS014","LU052_TLS010","LU052_TLS014", "LU083_TLS012", "LU083_TLS015")  

int_list_NR <- list(int_035_013_data, int_035_015_data, int_066_015_data,int_080_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data)
names(int_list_NR) <- c("LU035_TLS013","LU035_TLS015","LU066_TLS015","LU080", "LU086-1", "LU086-2", "LU088", "LU091")  

unstim_average_int_R <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_int_R <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_int_R <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_int_R <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(int_list_R)){
  data <- int_list_R[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_int_R <- rbind(unstim_average_int_R, unstim)
    rownames(unstim_average_int_R)[[nrow(unstim_average_int_R)]] <- names(int_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_int_R <- rbind(aPD1_average_int_R, aPD1)
    rownames(aPD1_average_int_R)[[nrow(aPD1_average_int_R)]] <- names(int_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_int_R <- rbind(aPD1_aIFNR_average_int_R, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_int_R)[[nrow(aPD1_aIFNR_average_int_R)]] <- names(int_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_int_R <- rbind(hrIFNy_average_int_R, hrIFNy)
    rownames(hrIFNy_average_int_R)[[nrow(hrIFNy_average_int_R)]] <- names(int_list_R)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_int_R) <- colnames(int_023_014_data)
unstim_average_int_R <- colMeans(unstim_average_int_R, na.rm = T)

colnames(aPD1_average_int_R) <- colnames(int_023_014_data)
aPD1_average_int_R <- colMeans(aPD1_average_int_R, na.rm = T)

colnames(aPD1_aIFNR_average_int_R) <- colnames(int_023_014_data)
aPD1_aIFNR_average_int_R <- colMeans(aPD1_aIFNR_average_int_R, na.rm = T)

colnames(hrIFNy_average_int_R) <- colnames(int_023_014_data)
hrIFNy_average_int_R <- colMeans(hrIFNy_average_int_R, na.rm = T)

#repeat for non-responders
unstim_average_int_NR <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_int_NR <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_int_NR <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_int_NR <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(int_list_NR)){
  data <- int_list_NR[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_int_NR <- rbind(unstim_average_int_NR, unstim)
    rownames(unstim_average_int_NR)[[nrow(unstim_average_int_NR)]] <- names(int_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_int_NR <- rbind(aPD1_average_int_NR, aPD1)
    rownames(aPD1_average_int_NR)[[nrow(aPD1_average_int_NR)]] <- names(int_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_int_NR <- rbind(aPD1_aIFNR_average_int_NR, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_int_NR)[[nrow(aPD1_aIFNR_average_int_NR)]] <- names(int_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_int_NR <- rbind(hrIFNy_average_int_NR, hrIFNy)
    rownames(hrIFNy_average_int_NR)[[nrow(hrIFNy_average_int_NR)]] <- names(int_list_NR)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_int_NR) <- colnames(int_023_014_data)
unstim_average_int_NR <- colMeans(unstim_average_int_NR, na.rm = T)

colnames(aPD1_average_int_NR) <- colnames(int_023_014_data)
aPD1_average_int_NR <- colMeans(aPD1_average_int_NR, na.rm = T)

colnames(aPD1_aIFNR_average_int_NR) <- colnames(int_023_014_data)
aPD1_aIFNR_average_int_NR <- colMeans(aPD1_aIFNR_average_int_NR, na.rm = T)

colnames(hrIFNy_average_int_NR) <- colnames(int_023_014_data)
hrIFNy_average_int_NR <- colMeans(hrIFNy_average_int_NR, na.rm = T)



#and tumour
tumour_list_R <- list(tumour_023_013_data,tumour_023_014_data,tumour_052_010_data, tumour_052_014_data,tumour_083_012_data, tumour_083_015_data)
names(tumour_list_R) <- c("LU023_TLS013","LU023_TLS014","LU052_TLS010","LU052_TLS014", "LU083_TLS012", "LU083_TLS015")  

tumour_list_NR <- list(tumour_035_013_data, tumour_035_015_data,tumour_066_012_data, tumour_066_015_data,tumour_080_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data)
names(tumour_list_NR) <- c("LU035_TLS013","LU035_TLS015","LU066_TLS012", "LU066_TLS015","LU080", "LU086-1", "LU086-2", "LU088", "LU091")  

unstim_average_tumour_R <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_tumour_R <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_tumour_R <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_tumour_R <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(tumour_list_R)){
  data <- tumour_list_R[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_tumour_R <- rbind(unstim_average_tumour_R, unstim)
    rownames(unstim_average_tumour_R)[[nrow(unstim_average_tumour_R)]] <- names(tumour_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_tumour_R <- rbind(aPD1_average_tumour_R, aPD1)
    rownames(aPD1_average_tumour_R)[[nrow(aPD1_average_tumour_R)]] <- names(tumour_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_tumour_R <- rbind(aPD1_aIFNR_average_tumour_R, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_tumour_R)[[nrow(aPD1_aIFNR_average_tumour_R)]] <- names(tumour_list_R)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_tumour_R <- rbind(hrIFNy_average_tumour_R, hrIFNy)
    rownames(hrIFNy_average_tumour_R)[[nrow(hrIFNy_average_tumour_R)]] <- names(tumour_list_R)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_tumour_R) <- colnames(tumour_023_014_data)
unstim_average_tumour_R <- colMeans(unstim_average_tumour_R, na.rm = T)

colnames(aPD1_average_tumour_R) <- colnames(tumour_023_014_data)
aPD1_average_tumour_R <- colMeans(aPD1_average_tumour_R, na.rm = T)

colnames(aPD1_aIFNR_average_tumour_R) <- colnames(tumour_023_014_data)
aPD1_aIFNR_average_tumour_R <- colMeans(aPD1_aIFNR_average_tumour_R, na.rm = T)

colnames(hrIFNy_average_tumour_R) <- colnames(tumour_023_014_data)
hrIFNy_average_tumour_R <- colMeans(hrIFNy_average_tumour_R, na.rm = T)

#repeat for non-responders
unstim_average_tumour_NR <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_average_tumour_NR <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_average_tumour_NR <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_average_tumour_NR <- data.frame(matrix(nrow = 0, ncol = 25))

#Calculate condition averages per tumour
for(i in seq_along(tumour_list_NR)){
  data <- tumour_list_NR[[i]]
  
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_average_tumour_NR <- rbind(unstim_average_tumour_NR, unstim)
    rownames(unstim_average_tumour_NR)[[nrow(unstim_average_tumour_NR)]] <- names(tumour_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_average_tumour_NR <- rbind(aPD1_average_tumour_NR, aPD1)
    rownames(aPD1_average_tumour_NR)[[nrow(aPD1_average_tumour_NR)]] <- names(tumour_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_average_tumour_NR <- rbind(aPD1_aIFNR_average_tumour_NR, aPD1_aIFNR)
    rownames(aPD1_aIFNR_average_tumour_NR)[[nrow(aPD1_aIFNR_average_tumour_NR)]] <- names(tumour_list_NR)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_average_tumour_NR <- rbind(hrIFNy_average_tumour_NR, hrIFNy)
    rownames(hrIFNy_average_tumour_NR)[[nrow(hrIFNy_average_tumour_NR)]] <- names(tumour_list_NR)[[i]]
  }
  
}

#Average across all tumours
colnames(unstim_average_tumour_NR) <- colnames(tumour_023_014_data)
unstim_average_tumour_NR <- colMeans(unstim_average_tumour_NR, na.rm = T)

colnames(aPD1_average_tumour_NR) <- colnames(tumour_023_014_data)
aPD1_average_tumour_NR <- colMeans(aPD1_average_tumour_NR, na.rm = T)

colnames(aPD1_aIFNR_average_tumour_NR) <- colnames(tumour_023_014_data)
aPD1_aIFNR_average_tumour_NR <- colMeans(aPD1_aIFNR_average_tumour_NR, na.rm = T)

colnames(hrIFNy_average_tumour_NR) <- colnames(tumour_023_014_data)
hrIFNy_average_tumour_NR <- colMeans(hrIFNy_average_tumour_NR, na.rm = T)




#combine in a heatmap
#first: responder unstim-aPD1-combi
heatmap <- rbind(unstim_average_TLS_R,unstim_average_int_R,unstim_average_tumour_R,aPD1_average_TLS_R,aPD1_average_int_R,aPD1_average_tumour_R,
                 aPD1_aIFNR_average_TLS_R,aPD1_aIFNR_average_int_R,aPD1_aIFNR_average_tumour_R)

rownames(heatmap) <- c("unstim TLS" ,"unstim intermediate","unstim tumour bed","aPD1 TLS","aPD1 intermediate","aPD1 tumour bed",
                       "aPD1+aIFNyR1 TLS","aPD1+aIFNyR1 intermediate","aPD1+aIFNyR1 tumour bed")  

#z-scoring
colSD <- apply(heatmap, 2, sd)
colMean <- colMeans(heatmap)
z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y" #to make sure rownames can be seen
rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"

#annotation
annotation <- data.frame(segment = rep(c("TLS", "intermediate", "tumor bed"),4), condition = rep(c("unstim","aPD1","aPD1+aIFNyR1","hrIFNy"), each = 3))
rownames(annotation) <- colnames(z_scores)

#heatmap
correlations <- pheatmap(z_scores, annotation_col = annotation,cluster_cols = F)

correlations

save_pheatmap(correlations, "~/Documents/Analyses/TLS/final_all_fragments/responder_triple.pdf")

#capture row order
row_order <- correlations$tree_row$order
values <- rownames(z_scores)[row_order]

#reorder z-scores based on row order
z_scores <- z_scores[order(match(rownames(z_scores), values)),]

#save z-scores
z_scores <- as.data.frame(z_scores)
write.xlsx(z_scores, "~/Documents/Analyses/TLS/final_all_fragments/z_scores_responder_triple.xlsx", rowNames = T)

#second: responder unstim vs hrIFNy
heatmap <- rbind(unstim_average_TLS_R,unstim_average_int_R,unstim_average_tumour_R,hrIFNy_average_TLS_R,hrIFNy_average_int_R,hrIFNy_average_tumour_R)

rownames(heatmap) <- c("unstim TLS" ,"unstim intermediate","unstim tumour bed","hrIFNy TLS","hrIFNy intermediate","hrIFNy tumour bed")  

#z-scoring
colSD <- apply(heatmap, 2, sd)
colMean <- colMeans(heatmap)
z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y" #to make sure rownames can be seen
rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"

#heatmap
correlations <- pheatmap(z_scores, cluster_cols = F)

correlations

#capture row order
row_order <- correlations$tree_row$order
values <- rownames(z_scores)[row_order]

#reorder z-scores based on row order
z_scores <- z_scores[order(match(rownames(z_scores), values)),]

#save z-scores
z_scores <- as.data.frame(z_scores)
write.xlsx(z_scores, "~/Documents/Analyses/TLS/final_all_fragments/z_scores_responder_double.xlsx", rowNames = T)

#third: non-responder unstim vs hrIFNy
heatmap <- rbind(unstim_average_TLS_NR,unstim_average_int_R,unstim_average_tumour_NR,hrIFNy_average_TLS_NR,hrIFNy_average_int_NR,hrIFNy_average_tumour_NR)

rownames(heatmap) <- c("unstim TLS" ,"unstim intermediate","unstim tumour bed","hrIFNy TLS","hrIFNy intermediate","hrIFNy tumour bed")  

#z-scoring
colSD <- apply(heatmap, 2, sd)
colMean <- colMeans(heatmap)
z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y" #to make sure rownames can be seen
rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"

#heatmap
correlations <- pheatmap(z_scores, cluster_cols = F)

correlations

#capture row order
row_order <- correlations$tree_row$order
values <- rownames(z_scores)[row_order]

#reorder z-scores based on row order
z_scores <- z_scores[order(match(rownames(z_scores), values)),]

#save z-scores
z_scores <- as.data.frame(z_scores)
write.xlsx(z_scores, "~/Documents/Analyses/TLS/final_all_fragments/z_scores_nonresponder_double.xlsx", rowNames = T)


#-----------
#per tumour
TLS_per_tumour_data <- list(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                            TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data)
names(TLS_per_tumour_data) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012","LU066_TLS015",
                                "LU083_TLS012", "LU083_TLS015", "LU088", "LU091")  

unstim_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_per_tumour_TLS <- data.frame(matrix(nrow = 0, ncol = 25))

for(i in seq_along(TLS_per_tumour_data)){
  data <- TLS_per_tumour_data[[i]]
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_per_tumour_TLS <- rbind(unstim_per_tumour_TLS, unstim)
    rownames(unstim_per_tumour_TLS)[[nrow(unstim_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_per_tumour_TLS <- rbind(aPD1_per_tumour_TLS, aPD1)
    rownames(aPD1_per_tumour_TLS)[[nrow(aPD1_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_per_tumour_TLS <- rbind(aPD1_aIFNR_per_tumour_TLS, aPD1_aIFNR)
    rownames(aPD1_aIFNR_per_tumour_TLS)[[nrow(aPD1_aIFNR_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_per_tumour_TLS <- rbind(hrIFNy_per_tumour_TLS, hrIFNy)
    rownames(hrIFNy_per_tumour_TLS)[[nrow(hrIFNy_per_tumour_TLS)]] <- names(TLS_per_tumour_data)[[i]]
  }
  
}

rownames(unstim_per_tumour_TLS) <- paste0("unstim TLS ", rownames(unstim_per_tumour_TLS))
colnames(unstim_per_tumour_TLS) <- colnames(TLS_023_014_data)
rownames(aPD1_per_tumour_TLS) <- paste0("aPD1 TLS ", rownames(aPD1_per_tumour_TLS))
colnames(aPD1_per_tumour_TLS) <- colnames(TLS_023_014_data)
rownames(aPD1_aIFNR_per_tumour_TLS) <- paste0("aPD1+aIFNyR1 TLS ", rownames(aPD1_aIFNR_per_tumour_TLS))
colnames(aPD1_aIFNR_per_tumour_TLS) <- colnames(TLS_023_014_data)
rownames(hrIFNy_per_tumour_TLS) <- paste0("hrIFNy TLS ", rownames(hrIFNy_per_tumour_TLS))
colnames(hrIFNy_per_tumour_TLS) <- colnames(TLS_023_014_data)

#for intermediate
int_per_tumour_data <- list(int_023_013_data,int_023_014_data, int_035_013_data, int_035_015_data, int_052_010_data, int_052_014_data, int_066_015_data,
                            int_080_data, int_083_012_data, int_083_015_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data) #missing LU066 (TLS012)
names(int_per_tumour_data) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014", "LU066_TLS015",
                                "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  

unstim_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_per_tumour_int <- data.frame(matrix(nrow = 0, ncol = 25))

for(i in seq_along(int_per_tumour_data)){
  data <- int_per_tumour_data[[i]]
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_per_tumour_int <- rbind(unstim_per_tumour_int, unstim)
    rownames(unstim_per_tumour_int)[[nrow(unstim_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_per_tumour_int <- rbind(aPD1_per_tumour_int, aPD1)
    rownames(aPD1_per_tumour_int)[[nrow(aPD1_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_per_tumour_int <- rbind(aPD1_aIFNR_per_tumour_int, aPD1_aIFNR)
    rownames(aPD1_aIFNR_per_tumour_int)[[nrow(aPD1_aIFNR_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_per_tumour_int <- rbind(hrIFNy_per_tumour_int, hrIFNy)
    rownames(hrIFNy_per_tumour_int)[[nrow(hrIFNy_per_tumour_int)]] <- names(int_per_tumour_data)[[i]]
  }
  
}

rownames(unstim_per_tumour_int) <- paste0("unstim intermediate ", rownames(unstim_per_tumour_int))
colnames(unstim_per_tumour_int) <- colnames(TLS_023_014_data)
rownames(aPD1_per_tumour_int) <- paste0("aPD1 intermediate ", rownames(aPD1_per_tumour_int))
colnames(aPD1_per_tumour_int) <- colnames(TLS_023_014_data)
rownames(aPD1_aIFNR_per_tumour_int) <- paste0("aPD1+aIFNyR1 intermediate ", rownames(aPD1_aIFNR_per_tumour_int))
colnames(aPD1_aIFNR_per_tumour_int) <- colnames(TLS_023_014_data)
rownames(hrIFNy_per_tumour_int) <- paste0("hrIFNy intermediate ", rownames(hrIFNy_per_tumour_int))
colnames(hrIFNy_per_tumour_int) <- colnames(TLS_023_014_data)

#for tumour
tumour_per_tumour_data <- list(tumour_023_013_data,tumour_023_014_data, tumour_035_013_data, tumour_035_015_data, tumour_052_010_data, tumour_052_014_data,tumour_066_012_data,
                               tumour_066_015_data, tumour_080_data, tumour_083_012_data, tumour_083_015_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data)
names(tumour_per_tumour_data) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012", "LU066_TLS015",
                                   "LU080","LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  

unstim_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
aPD1_aIFNR_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
hrIFNy_per_tumour_tumour <- data.frame(matrix(nrow = 0, ncol = 25))

for(i in seq_along(tumour_per_tumour_data)){
  data <- tumour_per_tumour_data[[i]]
  unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
  aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
  aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
  hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
  
  unstim <- colMeans(unstim, na.rm = T)
  aPD1 <- colMeans(aPD1, na.rm = T)
  aPD1_aIFNR <- colMeans(aPD1_aIFNR, na.rm = T)
  hrIFNy <- colMeans(hrIFNy, na.rm = T)
  
  if(!any(grepl("NaN", unstim))){
    unstim_per_tumour_tumour <- rbind(unstim_per_tumour_tumour, unstim)
    rownames(unstim_per_tumour_tumour)[[nrow(unstim_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1))){
    aPD1_per_tumour_tumour <- rbind(aPD1_per_tumour_tumour, aPD1)
    rownames(aPD1_per_tumour_tumour)[[nrow(aPD1_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", aPD1_aIFNR))){
    aPD1_aIFNR_per_tumour_tumour <- rbind(aPD1_aIFNR_per_tumour_tumour, aPD1_aIFNR)
    rownames(aPD1_aIFNR_per_tumour_tumour)[[nrow(aPD1_aIFNR_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
  }
  
  if(!any(grepl("NaN", hrIFNy))){
    hrIFNy_per_tumour_tumour <- rbind(hrIFNy_per_tumour_tumour, hrIFNy)
    rownames(hrIFNy_per_tumour_tumour)[[nrow(hrIFNy_per_tumour_tumour)]] <- names(tumour_per_tumour_data)[[i]]
  }
  
}

rownames(unstim_per_tumour_tumour) <- paste0("unstim tumour bed ", rownames(unstim_per_tumour_tumour))
colnames(unstim_per_tumour_tumour) <- colnames(TLS_023_014_data)
rownames(aPD1_per_tumour_tumour) <- paste0("aPD1 tumour bed ", rownames(aPD1_per_tumour_tumour))
colnames(aPD1_per_tumour_tumour) <- colnames(TLS_023_014_data)
rownames(aPD1_aIFNR_per_tumour_tumour) <- paste0("aPD1+aIFNyR1 tumour bed ", rownames(aPD1_aIFNR_per_tumour_tumour))
colnames(aPD1_aIFNR_per_tumour_tumour) <- colnames(TLS_023_014_data)
rownames(hrIFNy_per_tumour_tumour) <- paste0("hrIFNy tumour bed ", rownames(hrIFNy_per_tumour_tumour))
colnames(hrIFNy_per_tumour_tumour) <- colnames(TLS_023_014_data)

#combine in a heatmap
heatmap <- rbind(unstim_per_tumour_TLS,unstim_per_tumour_int,unstim_per_tumour_tumour,aPD1_per_tumour_TLS,aPD1_per_tumour_int,aPD1_per_tumour_tumour,
                 aPD1_aIFNR_per_tumour_TLS,aPD1_aIFNR_per_tumour_int,aPD1_aIFNR_per_tumour_tumour,hrIFNy_per_tumour_TLS,hrIFNy_per_tumour_int,hrIFNy_per_tumour_tumour)



#z-scoring
colSD <- apply(heatmap, 2, sd)
colMean <- colMeans(heatmap)
z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

#annotation (I gave up on this...too many samples..)
#annotation <- data.frame(segment = c(rep("TLS",8), rep("intermediate",12), rep("tumour bed", 15)), condition = rep(c("unstim","aPD1","aPD1+aIFNyR1","hrIFNy"), each = 3))
#rownames(annotation) <- colnames(z_scores)

rownames(z_scores)[rownames(z_scores) == "IFN-γ"] <- "IFN-y" #to make sure rownames can be seen
rownames(z_scores)[rownames(z_scores) == "TNF-α"] <- "TNF-a"

#heatmap
pheatmap(z_scores,cluster_cols = F)
pheatmap(z_scores,cluster_cols = T)


#--------------------------------
  ##TIME FOR FOLD CHANGES
  TLS_list <- list(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                   TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data) #missing LU080, LU086-1, LU086-2
  names(TLS_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012","LU066_TLS015",
                       "LU083_TLS012", "LU083_TLS015", "LU088", "LU091")  
  
  aPD1_FC_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_FC_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  colnames(aPD1_aIFNR_FC_TLS) <- colnames(TLS_023_014_data)
  hrIFNy_FC_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  colnames(hrIFNy_FC_TLS) <- colnames(TLS_023_014_data)
  
  #Calculate condition averages per tumour
  
  #Adding a pseudocount of 1 to avoid Inf (it messes up my averages)
  for(i in seq_along(TLS_list)){
    data <- TLS_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim_avg <- colMeans(unstim)
    aPD1_avg <- colMeans(aPD1)
    aPD1_aIFNR_avg <- colMeans(aPD1_aIFNR)
    hrIFNy_avg <- colMeans(hrIFNy)
    
    
    if(!any(grepl("NaN", aPD1_avg))){
      aPD1_FC <- (aPD1_avg+1)/(unstim_avg+1)
      aPD1_FC[aPD1_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      aPD1_FC_TLS <- rbind(aPD1_FC_TLS, aPD1_FC)
      rownames(aPD1_FC_TLS)[[nrow(aPD1_FC_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR_avg))){
      aPD1_aIFNR_FC <- (aPD1_aIFNR_avg+1)/(unstim_avg+1)
      aPD1_aIFNR_FC[aPD1_aIFNR_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      aPD1_aIFNR_FC_TLS <- rbind(aPD1_aIFNR_FC_TLS, aPD1_aIFNR_FC)
      rownames(aPD1_aIFNR_FC_TLS)[[nrow(aPD1_aIFNR_FC_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy_avg))){
      hrIFNy_FC <- (hrIFNy_avg+1)/(unstim_avg+1)
      hrIFNy_FC[hrIFNy_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      hrIFNy_FC_TLS <- rbind(hrIFNy_FC_TLS, hrIFNy_FC)
      rownames(hrIFNy_FC_TLS)[[nrow(hrIFNy_FC_TLS)]] <- names(TLS_list)[[i]]
    }
    
  }
  
  #Average over all tumours
  aPD1_FC_TLS <- colMeans(aPD1_FC_TLS)
  names(aPD1_FC_TLS) <- colnames(TLS_023_014_data)
  
  aPD1_aIFNR_FC_TLS <- colMeans(aPD1_aIFNR_FC_TLS)
  names(aPD1_aIFNR_FC_TLS) <- colnames(TLS_023_014_data)
  
  hrIFNy_FC_TLS <- colMeans(hrIFNy_FC_TLS)
  names(hrIFNy_FC_TLS) <- colnames(TLS_023_014_data)
  
  #repeat for intermediate
  aPD1_FC_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_FC_int <- data.frame(matrix(nrow = 0, ncol = 25))
  colnames(aPD1_aIFNR_FC_int) <- colnames(int_023_014_data)
  hrIFNy_FC_int <- data.frame(matrix(nrow = 0, ncol = 25))
  colnames(hrIFNy_FC_int) <- colnames(int_023_014_data)

  for(i in seq_along(int_list)){
    data <- int_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim_avg <- colMeans(unstim)
    aPD1_avg <- colMeans(aPD1)
    aPD1_aIFNR_avg <- colMeans(aPD1_aIFNR)
    hrIFNy_avg <- colMeans(hrIFNy)
    
    
    if(!any(grepl("NaN", aPD1_avg))){
      aPD1_FC <- (aPD1_avg+1)/(unstim_avg+1)
      aPD1_FC[aPD1_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      aPD1_FC_int <- rbind(aPD1_FC_int, aPD1_FC)
      rownames(aPD1_FC_int)[[nrow(aPD1_FC_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR_avg))){
      aPD1_aIFNR_FC <- (aPD1_aIFNR_avg+1)/(unstim_avg+1)
      aPD1_aIFNR_FC[aPD1_aIFNR_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      aPD1_aIFNR_FC_int <- rbind(aPD1_aIFNR_FC_int, aPD1_aIFNR_FC)
      rownames(aPD1_aIFNR_FC_int)[[nrow(aPD1_aIFNR_FC_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy_avg))){
      hrIFNy_FC <- (hrIFNy_avg+1)/(unstim_avg+1)
      hrIFNy_FC[hrIFNy_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      hrIFNy_FC_int <- rbind(hrIFNy_FC_int, hrIFNy_FC)
      rownames(hrIFNy_FC_int)[[nrow(hrIFNy_FC_int)]] <- names(int_list)[[i]]
    }
    
  }
  
  #Average over all tumours
  aPD1_FC_int <- colMeans(aPD1_FC_int)
  names(aPD1_FC_int) <- colnames(int_023_014_data)
  
  aPD1_aIFNR_FC_int <- colMeans(aPD1_aIFNR_FC_int)
  names(aPD1_aIFNR_FC_int) <- colnames(int_023_014_data)
  
  hrIFNy_FC_int <- colMeans(hrIFNy_FC_int)
  names(hrIFNy_FC_int) <- colnames(int_023_014_data)
  
  #repeat for tumour bed
  aPD1_FC_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_FC_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  colnames(aPD1_aIFNR_FC_tumour) <- colnames(tumour_023_014_data)
  hrIFNy_FC_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  colnames(hrIFNy_FC_tumour) <- colnames(tumour_023_014_data)

  for(i in seq_along(tumour_list)){
    data <- tumour_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim_avg <- colMeans(unstim)
    aPD1_avg <- colMeans(aPD1)
    aPD1_aIFNR_avg <- colMeans(aPD1_aIFNR)
    hrIFNy_avg <- colMeans(hrIFNy)
    
    
    if(!any(grepl("NaN", aPD1_avg))){
      aPD1_FC <- (aPD1_avg+1)/(unstim_avg+1)
      aPD1_FC[aPD1_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      aPD1_FC_tumour <- rbind(aPD1_FC_tumour, aPD1_FC)
      rownames(aPD1_FC_tumour)[[nrow(aPD1_FC_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR_avg))){
      aPD1_aIFNR_FC <- (aPD1_aIFNR_avg+1)/(unstim_avg+1)
      aPD1_aIFNR_FC[aPD1_aIFNR_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      aPD1_aIFNR_FC_tumour <- rbind(aPD1_aIFNR_FC_tumour, aPD1_aIFNR_FC)
      rownames(aPD1_aIFNR_FC_tumour)[[nrow(aPD1_aIFNR_FC_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy_avg))){
      hrIFNy_FC <- (hrIFNy_avg+1)/(unstim_avg+1)
      hrIFNy_FC[hrIFNy_FC == "NaN"] <- 0 #occurs when both aPD1 and unstim are zero
      hrIFNy_FC_tumour <- rbind(hrIFNy_FC_tumour, hrIFNy_FC)
      rownames(hrIFNy_FC_tumour)[[nrow(hrIFNy_FC_tumour)]] <- names(tumour_list)[[i]]
    }
    
  }
  
  #Average over all tumours
  aPD1_FC_tumour <- colMeans(aPD1_FC_tumour)
  names(aPD1_FC_tumour) <- colnames(tumour_023_014_data)
  
  aPD1_aIFNR_FC_tumour <- colMeans(aPD1_aIFNR_FC_tumour)
  names(aPD1_aIFNR_FC_tumour) <- colnames(tumour_023_014_data)
  
  hrIFNy_FC_tumour <- colMeans(hrIFNy_FC_tumour)
  names(hrIFNy_FC_tumour) <- colnames(tumour_023_014_data)
  
  #combine in a heatmap
  heatmap <- rbind(aPD1_FC_TLS,aPD1_FC_int,aPD1_FC_tumour,aPD1_aIFNR_FC_TLS,aPD1_aIFNR_FC_int,aPD1_aIFNR_FC_tumour,
                   hrIFNy_FC_TLS,hrIFNy_FC_int,hrIFNy_FC_tumour)
  
  library(gtools)
  for(i in 1:nrow(heatmap)){
    heatmap[i,] <- foldchange2logratio(heatmap[i,])
  }
  
  rownames(heatmap) <- c("aPD1 TLS","aPD1 intermediate","aPD1 tumour bed", "aPD1+aIFNyR1 TLS","aPD1+aIFNyR1 intermediate","aPD1+aIFNyR1 tumour bed",
                         "hrIFNy TLS","hrIFNy intermediate","hrIFNy tumour bed")  
  
  
  #z-scoring
  colSD <- apply(heatmap, 2, sd)
  colMean <- colMeans(heatmap)
  z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
  z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
  z_scores <- t(z_scores)
  pheatmap(z_scores, cluster_cols = F,legend_labels = c("-2","-1","0","1","2", "z_score\n") ,main = "")
  
  #without z-scoring:
  fold_change <- t(heatmap)
  annotation$sample <- NULL
  pheatmap(fold_change, annotation_col = annotation,cluster_cols = F)
  
  
  #--------------------------------
  #Heatmap based on deltas
  #read in all tumour values
  lu023_TLS014_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023_TLS014/LU023_TLS014.xlsx", sheet = 3)) #deltas
  rownames(lu023_TLS014_delta) <- lu023_TLS014_delta[,1]
  lu023_TLS014_delta[,1] <- NULL
  rownames(lu023_TLS014_delta) <- gsub("Unstim", "unstim", rownames(lu023_TLS014_delta))
  
  lu088_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023_TLS014/LU088.xlsx", sheet = 3)) #deltas
  rownames(lu088_delta) <- lu088_delta[,1]
  lu088_delta[,1] <- NULL
  rownames(lu088_delta) <- gsub("Unstim", "unstim", rownames(lu088_delta))
  
  lu091_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023_TLS014/LU091.xlsx", sheet = 3)) #deltas
  rownames(lu091_delta) <- lu091_delta[,1]
  lu091_delta[,1] <- NULL
  rownames(lu091_delta) <- gsub("Unstim", "unstim", rownames(lu091_delta))
  
  lu035_TLS013_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS013.xlsx", sheet = 3)) #deltas
  rownames(lu035_TLS013_delta) <- lu035_TLS013_delta[,1]
  lu035_TLS013_delta[,1] <- NULL
  
  lu035_TLS015_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS015.xlsx", sheet = 3)) #deltas
  rownames(lu035_TLS015_delta) <- lu035_TLS015_delta[,1]
  lu035_TLS015_delta[,1] <- NULL
  
  lu066_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU066.xlsx", sheet = 3)) #deltas
  rownames(lu066_delta) <- lu066_delta[,1]
  lu066_delta[,1] <- NULL
  
  lu052_TLS010_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS010.xlsx", sheet = 3)) #deltas
  rownames(lu052_TLS010_delta) <- lu052_TLS010_delta[,1]
  lu052_TLS010_delta[,1] <- NULL
  
  lu052_TLS014_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS014.xlsx", sheet = 3)) #deltas
  rownames(lu052_TLS014_delta) <- lu052_TLS014_delta[,1]
  lu052_TLS014_delta[,1] <- NULL
  
  #Get separate Legendplex scores per segment
  TLS_023_014_delta <- lu023_TLS014_delta[rownames(lu023_TLS014_delta) %in% TLS_023_014,]
  TLS_023_014_delta <- colMeans(TLS_023_014_delta) #average score per cytokine across all conditions
  int_023_014_delta <- lu023_TLS014_delta[rownames(lu023_TLS014_delta) %in% int_023_014,]
  int_023_014_delta <- colMeans(int_023_014_delta)
  tumour_023_014_delta <- lu023_TLS014_delta[rownames(lu023_TLS014_delta) %in% tumour_023_014,]
  tumour_023_014_delta <- colMeans(tumour_023_014_delta)
  
  TLS_088_delta <- lu088_delta[rownames(lu088_delta) %in% TLS_088,]
  TLS_088_delta <- colMeans(TLS_088_delta) #average score per cytokine across all conditions
  int_088_delta <- lu088_delta[rownames(lu088_delta) %in% int_088,]
  int_088_delta <- colMeans(int_088_delta)
  tumour_088_delta <- lu088_delta[rownames(lu088_delta) %in% tumour_088,]
  tumour_088_delta <- colMeans(tumour_088_delta)
  
  TLS_091_delta <- lu091_delta[rownames(lu091_delta) %in% TLS_091,]
  TLS_091_delta <- colMeans(TLS_091_delta) #average score per cytokine across all conditions
  int_091_delta <- lu091_delta[rownames(lu091_delta) %in% int_091,]
  int_091_delta <- colMeans(int_091_delta)
  tumour_091_delta <- lu091_delta[rownames(lu091_delta) %in% tumour_091,]
  tumour_091_delta <- colMeans(tumour_091_delta)
  
  TLS_035_013_delta <- lu035_TLS013_delta[rownames(lu035_TLS013_delta) %in% TLS_035_013,]
  TLS_035_013_delta <- colMeans(TLS_035_013_delta) #average score per cytokine across all conditions
  int_035_013_delta <- lu035_TLS013_delta[rownames(lu035_TLS013_delta) %in% int_035_013,]
  int_035_013_delta <- colMeans(int_035_013_delta)
  tumour_035_013_delta <- lu035_TLS013_delta[rownames(lu035_TLS013_delta) %in% tumour_035_013,]
  tumour_035_013_delta <- colMeans(tumour_035_013_delta)
  
  TLS_035_015_delta <- lu035_TLS015_delta[rownames(lu035_TLS015_delta) %in% TLS_035_015,]
  TLS_035_015_delta <- colMeans(TLS_035_015_delta) #average score per cytokine across all conditions
  int_035_015_delta <- lu035_TLS015_delta[rownames(lu035_TLS015_delta) %in% int_035_015,]
  int_035_015_delta <- colMeans(int_035_015_delta)
  tumour_035_015_delta <- lu035_TLS015_delta[rownames(lu035_TLS015_delta) %in% tumour_035_015,]
  tumour_035_015_delta <- colMeans(tumour_035_015_delta)
  
  TLS_066_015_delta <- lu066_delta[rownames(lu066_delta) %in% TLS_066_015,]
  TLS_066_015_delta <- colMeans(TLS_066_015_delta) #average score per cytokine across all conditions
  int_066_015_delta <- lu066_delta[rownames(lu066_delta) %in% int_066_015,]
  int_066_015_delta <- colMeans(int_066_015_delta)
  tumour_066_015_delta <- lu066_delta[rownames(lu066_delta) %in% tumour_066_015,]
  tumour_066_015_delta <- colMeans(tumour_066_015_delta)
  
  TLS_052_010_delta <- lu052_TLS010_delta[rownames(lu052_TLS010_delta) %in% TLS_052_010,]
  TLS_052_010_delta <- colMeans(TLS_052_010_delta) #average score per cytokine across all conditions
  int_052_010_delta <- lu052_TLS010_delta[rownames(lu052_TLS010_delta) %in% int_052_010,]
  int_052_010_delta <- colMeans(int_052_010_delta)
  tumour_052_010_delta <- lu052_TLS010_delta[rownames(lu052_TLS010_delta) %in% tumour_052_010,]
  tumour_052_010_delta <- colMeans(tumour_052_010_delta)
  
  TLS_052_014_delta <- lu052_TLS014_delta[rownames(lu052_TLS014_delta) %in% TLS_052_014,]
  TLS_052_014_delta <- colMeans(TLS_052_014_delta) #average score per cytokine across all conditions
  int_052_014_delta <- lu052_TLS014_delta[rownames(lu052_TLS014_delta) %in% int_052_014,]
  int_052_014_delta <- colMeans(int_052_014_delta)
  tumour_052_014_delta <- lu052_TLS014_delta[rownames(lu052_TLS014_delta) %in% tumour_052_014,]
  tumour_052_014_delta <- colMeans(tumour_052_014_delta)
  
  #Average across all tumours
  TLS_scores <- rbind(TLS_023_014_delta, TLS_035_013_delta, TLS_035_015_delta, TLS_088_delta, TLS_091_delta, TLS_066_015_delta, TLS_052_010_delta, TLS_052_014_delta)
  TLS_scores <- colMeans(TLS_scores)
  
  int_scores <- rbind(int_023_014_delta, int_035_013_delta, int_035_015_delta, int_088_delta, int_091_delta, int_066_015_delta, int_052_010_delta, int_052_014_delta)
  int_scores <- colMeans(int_scores)
  
  tumour_scores <- rbind(tumour_023_014_delta,tumour_035_013_delta,tumour_035_015_delta,tumour_088_delta,tumour_091_delta,tumour_066_015_delta, tumour_052_010_delta, tumour_052_014_delta)
  tumour_scores <- colMeans(tumour_scores)
  
  heatmap_scores <- rbind(TLS_scores, int_scores, tumour_scores)
  rownames(heatmap_scores) <- c("TLS", "intermediate", "tumour bed")
  
  #Generate a heatmap of z-scores per cytokine
  library(pheatmap)
  #z-score is (value - average)/standard deviation. (is this the delta or is the average per condition?)
  colSD <- apply(heatmap_scores, 2, sd)
  colMean <- colMeans(heatmap_scores)
  z_scores <- sweep(heatmap_scores, 2, colMean, "-") #subtract means per column
  z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
  z_scores <- t(z_scores)
  
  
  #heatmap
  heatmap <- suppressWarnings(pheatmap(z_scores))  
  
  
  #split per condition
  lu023_TLS014_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023_TLS014/LU023_TLS014.xlsx", sheet = 3)) #deltas
  rownames(lu023_TLS014_delta) <- lu023_TLS014_delta[,1]
  lu023_TLS014_delta[,1] <- NULL
  rownames(lu023_TLS014_delta) <- gsub("Unstim", "unstim", rownames(lu023_TLS014_delta))
  
  lu088_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023_TLS014/LU088.xlsx", sheet = 3)) #deltas
  rownames(lu088_delta) <- lu088_delta[,1]
  lu088_delta[,1] <- NULL
  rownames(lu088_delta) <- gsub("Unstim", "unstim", rownames(lu088_delta))
  
  lu091_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP063-LU088-LU091-LU023_TLS014/LU091.xlsx", sheet = 3)) #deltas
  rownames(lu091_delta) <- lu091_delta[,1]
  lu091_delta[,1] <- NULL
  rownames(lu091_delta) <- gsub("Unstim", "unstim", rownames(lu091_delta))
  
  lu035_TLS013_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS013.xlsx", sheet = 3)) #deltas
  rownames(lu035_TLS013_delta) <- lu035_TLS013_delta[,1]
  lu035_TLS013_delta[,1] <- NULL
  
  lu035_TLS015_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP064-Timm-LU035(2)/LU035_TLS015.xlsx", sheet = 3)) #deltas
  rownames(lu035_TLS015_delta) <- lu035_TLS015_delta[,1]
  lu035_TLS015_delta[,1] <- NULL
  
  lu066_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU066.xlsx", sheet = 3)) #deltas
  rownames(lu066_delta) <- lu066_delta[,1]
  lu066_delta[,1] <- NULL
  
  lu052_TLS010_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS010.xlsx", sheet = 3)) #deltas
  rownames(lu052_TLS010_delta) <- lu052_TLS010_delta[,1]
  lu052_TLS010_delta[,1] <- NULL
  
  lu052_TLS014_delta <- as.data.frame(read_excel("~/Documents/Experiments/TLS/LP065-LU066-LU052(2)/LU052_TLS014.xlsx", sheet = 3)) #deltas
  rownames(lu052_TLS014_delta) <- lu052_TLS014_delta[,1]
  lu052_TLS014_delta[,1] <- NULL
  
  #Get separate Legendplex scores per segment
  TLS_023_014_delta <- lu023_TLS014_delta[rownames(lu023_TLS014_delta) %in% TLS_023_014,]
  int_023_014_delta <- lu023_TLS014_delta[rownames(lu023_TLS014_delta) %in% int_023_014,]
  tumour_023_014_delta <- lu023_TLS014_delta[rownames(lu023_TLS014_delta) %in% tumour_023_014,]

  TLS_088_delta <- lu088_delta[rownames(lu088_delta) %in% TLS_088,]
  int_088_delta <- lu088_delta[rownames(lu088_delta) %in% int_088,]
  tumour_088_delta <- lu088_delta[rownames(lu088_delta) %in% tumour_088,]

  TLS_091_delta <- lu091_delta[rownames(lu091_delta) %in% TLS_091,]
  int_091_delta <- lu091_delta[rownames(lu091_delta) %in% int_091,]
  tumour_091_delta <- lu091_delta[rownames(lu091_delta) %in% tumour_091,]

  TLS_035_013_delta <- lu035_TLS013_delta[rownames(lu035_TLS013_delta) %in% TLS_035_013,]
  int_035_013_delta <- lu035_TLS013_delta[rownames(lu035_TLS013_delta) %in% int_035_013,]
  tumour_035_013_delta <- lu035_TLS013_delta[rownames(lu035_TLS013_delta) %in% tumour_035_013,]

  TLS_035_015_delta <- lu035_TLS015_delta[rownames(lu035_TLS015_delta) %in% TLS_035_015,]
  int_035_015_delta <- lu035_TLS015_delta[rownames(lu035_TLS015_delta) %in% int_035_015,]
  tumour_035_015_delta <- lu035_TLS015_delta[rownames(lu035_TLS015_delta) %in% tumour_035_015,]

  TLS_066_015_delta <- lu066_delta[rownames(lu066_delta) %in% TLS_066_015,]
  int_066_015_delta <- lu066_delta[rownames(lu066_delta) %in% int_066_015,]
  tumour_066_015_delta <- lu066_delta[rownames(lu066_delta) %in% tumour_066_015,]

  TLS_052_010_delta <- lu052_TLS010_delta[rownames(lu052_TLS010_delta) %in% TLS_052_010,]
  int_052_010_delta <- lu052_TLS010_delta[rownames(lu052_TLS010_delta) %in% int_052_010,]
  tumour_052_010_delta <- lu052_TLS010_delta[rownames(lu052_TLS010_delta) %in% tumour_052_010,]

  TLS_052_014_delta <- lu052_TLS014_delta[rownames(lu052_TLS014_delta) %in% TLS_052_014,]
  int_052_014_delta <- lu052_TLS014_delta[rownames(lu052_TLS014_delta) %in% int_052_014,]
  tumour_052_014_delta <- lu052_TLS014_delta[rownames(lu052_TLS014_delta) %in% tumour_052_014,]

  #First for TLS
  TLS_list <- list(TLS_023_014_data, TLS_035_013_data,TLS_035_015_data,TLS_088_data,TLS_091_data,TLS_066_015_data, TLS_052_010_data, TLS_052_014_data)
  names(TLS_list) <- c("LU023_TLS014","LU035_TLS013","LU035_TLS015","LU088", "LU091","LU066","LU052_TLS010","LU052_TLS014")  
  
  unstim_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Calculate condition averages per tumour
  for(i in seq_along(TLS_list)){
    data <- TLS_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
      unstim_average_TLS <- rbind(unstim_average_TLS, unstim)
      rownames(unstim_average_TLS)[[nrow(unstim_average_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_average_TLS <- rbind(aPD1_average_TLS, aPD1)
      rownames(aPD1_average_TLS)[[nrow(aPD1_average_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_average_TLS <- rbind(aPD1_aIFNR_average_TLS, aPD1_aIFNR)
      rownames(aPD1_aIFNR_average_TLS)[[nrow(aPD1_aIFNR_average_TLS)]] <- names(TLS_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_average_TLS <- rbind(hrIFNy_average_TLS, hrIFNy)
      rownames(hrIFNy_average_TLS)[[nrow(hrIFNy_average_TLS)]] <- names(TLS_list)[[i]]
    }
    
  }
  
  #Average across all tumours
  colnames(unstim_average_TLS) <- colnames(TLS_023_014_delta)
  unstim_average_TLS <- colMeans(unstim_average_TLS)
  
  colnames(aPD1_average_TLS) <- colnames(TLS_023_014_delta)
  aPD1_average_TLS <- colMeans(aPD1_average_TLS)
  
  colnames(aPD1_aIFNR_average_TLS) <- colnames(TLS_023_014_delta)
  aPD1_aIFNR_average_TLS <- colMeans(aPD1_aIFNR_average_TLS)
  
  colnames(hrIFNy_average_TLS) <- colnames(TLS_023_014_delta)
  hrIFNy_average_TLS <- colMeans(hrIFNy_average_TLS)
  
  #repeat for intermediate and tumour sections
  int_list <- list(int_023_014_delta, int_035_013_delta,int_035_015_delta,int_088_delta,int_091_delta, int_066_015_delta, int_052_010_delta, int_052_014_delta)
  names(int_list) <- c("LU023_TLS014","LU035_TLS013","LU035_TLS015","LU088", "LU091","LU066","LU052_TLS010","LU052_TLS014")  
  
  unstim_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Calculate condition averages per tumour
  for(i in seq_along(int_list)){
    data <- int_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
      unstim_average_int <- rbind(unstim_average_int, unstim)
      rownames(unstim_average_int)[[nrow(unstim_average_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_average_int <- rbind(aPD1_average_int, aPD1)
      rownames(aPD1_average_int)[[nrow(aPD1_average_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_average_int <- rbind(aPD1_aIFNR_average_int, aPD1_aIFNR)
      rownames(aPD1_aIFNR_average_int)[[nrow(aPD1_aIFNR_average_int)]] <- names(int_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_average_int <- rbind(hrIFNy_average_int, hrIFNy)
      rownames(hrIFNy_average_int)[[nrow(hrIFNy_average_int)]] <- names(int_list)[[i]]
    }
    
  }
  
  #Average across all tumours
  colnames(unstim_average_int) <- colnames(TLS_023_014_delta)
  unstim_average_int <- colMeans(unstim_average_int)
  
  colnames(aPD1_average_int) <- colnames(TLS_023_014_delta)
  aPD1_average_int <- colMeans(aPD1_average_int)
  
  colnames(aPD1_aIFNR_average_int) <- colnames(TLS_023_014_delta)
  aPD1_aIFNR_average_int <- colMeans(aPD1_aIFNR_average_int)
  
  colnames(hrIFNy_average_int) <- colnames(TLS_023_014_delta)
  hrIFNy_average_int <- colMeans(hrIFNy_average_int)
  
  #and tumour
  tumour_list <- list(tumour_023_014_delta, tumour_035_013_delta,tumour_035_015_delta,tumour_088_delta,tumour_091_delta, tumour_066_015_delta, tumour_052_010_delta, tumour_052_014_delta)
  names(tumour_list) <- c("LU023_TLS014","LU035_TLS013","LU035_TLS015","LU088", "LU091","LU066","LU052_TLS010","LU052_TLS014")  
  
  unstim_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Calculate condition averages per tumour
  for(i in seq_along(tumour_list)){
    data <- tumour_list[[i]]
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    unstim <- colMeans(unstim)
    aPD1 <- colMeans(aPD1)
    aPD1_aIFNR <- colMeans(aPD1_aIFNR)
    hrIFNy <- colMeans(hrIFNy)
    
    if(!any(grepl("NaN", unstim))){
      unstim_average_tumour <- rbind(unstim_average_tumour, unstim)
      rownames(unstim_average_tumour)[[nrow(unstim_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1))){
      aPD1_average_tumour <- rbind(aPD1_average_tumour, aPD1)
      rownames(aPD1_average_tumour)[[nrow(aPD1_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", aPD1_aIFNR))){
      aPD1_aIFNR_average_tumour <- rbind(aPD1_aIFNR_average_tumour, aPD1_aIFNR)
      rownames(aPD1_aIFNR_average_tumour)[[nrow(aPD1_aIFNR_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
    if(!any(grepl("NaN", hrIFNy))){
      hrIFNy_average_tumour <- rbind(hrIFNy_average_tumour, hrIFNy)
      rownames(hrIFNy_average_tumour)[[nrow(hrIFNy_average_tumour)]] <- names(tumour_list)[[i]]
    }
    
  }
  
  #Average across all tumours
  colnames(unstim_average_tumour) <- colnames(TLS_023_014_delta)
  unstim_average_tumour <- colMeans(unstim_average_tumour)
  
  colnames(aPD1_average_tumour) <- colnames(TLS_023_014_delta)
  aPD1_average_tumour <- colMeans(aPD1_average_tumour)
  
  colnames(aPD1_aIFNR_average_tumour) <- colnames(TLS_023_014_delta)
  aPD1_aIFNR_average_tumour <- colMeans(aPD1_aIFNR_average_tumour)
  
  colnames(hrIFNy_average_tumour) <- colnames(TLS_023_014_delta)
  hrIFNy_average_tumour <- colMeans(hrIFNy_average_tumour)
  
  
  
  
  #combine in a heatmap
  heatmap <- rbind(unstim_average_TLS,unstim_average_int,unstim_average_tumour,aPD1_average_TLS,aPD1_average_int,aPD1_average_tumour,
                   aPD1_aIFNR_average_TLS,aPD1_aIFNR_average_int,aPD1_aIFNR_average_tumour,hrIFNy_average_TLS,hrIFNy_average_int,hrIFNy_average_tumour)
  
  rownames(heatmap) <- c("unstim TLS" ,"unstim intermediate","unstim tumour bed","aPD1 TLS","aPD1 intermediate","aPD1 tumour bed",
                         "aPD1+aIFNyR1 TLS","aPD1+aIFNyR1 intermediate","aPD1+aIFNyR1 tumour bed","hrIFNy TLS","hrIFNy intermediate","hrIFNy tumour bed")  
  
  #z-scoring
  colSD <- apply(heatmap, 2, sd)
  colMean <- colMeans(heatmap)
  z_scores <- sweep(heatmap, 2, colMean, "-") #subtract means per column
  z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
  z_scores <- t(z_scores)
  
  #annotation
  annotation <- data.frame(segment = rep(c("TLS", "intermediate", "tumor bed"),4), condition = rep(c("unstim","aPD1","aPD1+aIFNyR1","hrIFNy"), each = 3))
  rownames(annotation) <- colnames(z_scores)
  
  #heatmap
  pheatmap(z_scores, annotation_col = annotation,cluster_cols = F)
  
  #------------
  #Generate cumulative scores per fragment (unstim, compare segments)
  TLS_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% TLS_023_013,]
  int_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% int_023_013,]
  tumour_023_013_data <- lu023_TLS013[rownames(lu023_TLS013) %in% tumour_023_013,]
  
  TLS_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% TLS_023_014,]
  int_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% int_023_014,]
  tumour_023_014_data <- lu023_TLS014[rownames(lu023_TLS014) %in% tumour_023_014,]
  
  TLS_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% TLS_035_013,]
  int_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% int_035_013,]
  tumour_035_013_data <- lu035_TLS013[rownames(lu035_TLS013) %in% tumour_035_013,]
  
  TLS_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% TLS_035_015,]
  int_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% int_035_015,]
  tumour_035_015_data <- lu035_TLS015[rownames(lu035_TLS015) %in% tumour_035_015,]
  
  TLS_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% TLS_052_010,]
  int_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% int_052_010,]
  tumour_052_010_data <- lu052_TLS010[rownames(lu052_TLS010) %in% tumour_052_010,]
  #excluding IFNy_high fragments:
  TLS_052_010_data <- TLS_052_010_data[!grepl("hrIFNy_high", rownames(TLS_052_010_data)),]
  int_052_010_data <- int_052_010_data[!grepl("hrIFNy_high", rownames(int_052_010_data)),]
  tumour_052_010_data <- tumour_052_010_data[!grepl("hrIFNy_high", rownames(tumour_052_010_data)),]
  
  TLS_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% TLS_052_014,]
  int_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% int_052_014,]
  tumour_052_014_data <- lu052_TLS014[rownames(lu052_TLS014) %in% tumour_052_014,]
  
  #This tumour has no intermediate fragments. Therefore, this condition is skipped.
  TLS_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% TLS_066_012,]
  tumour_066_012_data <- lu066_TLS012[rownames(lu066_TLS012) %in% tumour_066_012,]
  
  TLS_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% TLS_066_015,]
  int_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% int_066_015,]
  tumour_066_015_data <- lu066_TLS015[rownames(lu066_TLS015) %in% tumour_066_015,]
  
  #LU080 has no TLS so we ignore this condition
  int_080_data <- lu080[rownames(lu080) %in% int_080,]
  tumour_080_data <- lu080[rownames(lu080) %in% tumour_080,]
  
  TLS_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% TLS_083_012,]
  int_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% int_083_012,]
  tumour_083_012_data <- lu083_TLS012[rownames(lu083_TLS012) %in% tumour_083_012,]
  
  TLS_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% TLS_083_015,]
  int_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% int_083_015,]
  tumour_083_015_data <- lu083_TLS015[rownames(lu083_TLS015) %in% tumour_083_015,]
  
  #the only TLS of LU086-1 were in RN440-15/RN440-16, which will not be considered. Therefore, we ignore this condition.
  int_086_1_data <- lu086_1[rownames(lu086_1) %in% int_086_1,]
  tumour_086_1_data <- lu086_1[rownames(lu086_1) %in% tumour_086_1,]
  
  #the only TLS of LU086-2 was in RN440-15, which was not taken along in the Legendplex. Therefore, we ignore this condition.
  int_086_2_data <- lu086_2[rownames(lu086_2) %in% int_086_2,]
  tumour_086_2_data <- lu086_2[rownames(lu086_2) %in% tumour_086_2,]
  
  TLS_088_data <- lu088[rownames(lu088) %in% TLS_088,]
  int_088_data <- lu088[rownames(lu088) %in% int_088,]
  tumour_088_data <- lu088[rownames(lu088) %in% tumour_088,]
  
  TLS_091_data <- lu091[rownames(lu091) %in% TLS_091,]
  int_091_data <- lu091[rownames(lu091) %in% int_091,]
  tumour_091_data <- lu091[rownames(lu091) %in% tumour_091,]
  
  
  #First for TLS
  TLS_list <- list(TLS_023_013_data,TLS_023_014_data, TLS_035_013_data, TLS_035_015_data, TLS_052_010_data, TLS_052_014_data, TLS_066_012_data,TLS_066_015_data,
                   TLS_083_012_data, TLS_083_015_data, TLS_088_data, TLS_091_data) #missing LU080, LU086-1, LU086-2
  names(TLS_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012","LU066_TLS015",
                       "LU083_TLS012", "LU083_TLS015", "LU088", "LU091")  
  
  unstim_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_TLS <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Separate per condition
  for(i in seq_along(TLS_list)){
    data <- TLS_list[[i]]
    
    
    if(grepl("TLS", names(TLS_list)[[i]])){ #if a tumour was spread over multiple experiments
      if(!grepl("LU035|LU052", names(TLS_list)[[i]])){ #and if that tumour was not present in the same legendplex (bc those already have the TLS number in front)
        number <- strsplit(names(TLS_list)[[i]], "_")
        number <- number[[1]][[2]]
        rownames(data) <- paste0(number, " ", rownames(data)) #add the experiment number in front of the rownames
      }
    }
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    
    if(nrow(unstim)!=0){
      unstim_average_TLS <- rbind(unstim_average_TLS, unstim)
      colnames(unstim_average_TLS) <- colnames(data)
    }
    
    if(nrow(aPD1)!=0){
      aPD1_average_TLS <- rbind(aPD1_average_TLS, aPD1)
      colnames(aPD1_average_TLS) <- colnames(data)
    }
    
    if(nrow(aPD1_aIFNR)!=0){
      aPD1_aIFNR_average_TLS <- rbind(aPD1_aIFNR_average_TLS, aPD1_aIFNR)
      colnames(aPD1_aIFNR_average_TLS) <- colnames(data)
    }
    
    if(nrow(hrIFNy)!=0){
      hrIFNy_average_TLS <- rbind(hrIFNy_average_TLS, hrIFNy)
      colnames(hrIFNy_average_TLS) <- colnames(data)
    }
    
  }
  

  #for intermediate
  int_list <- list(int_023_013_data,int_023_014_data, int_035_013_data, int_035_015_data, int_052_010_data, int_052_014_data,int_066_015_data, int_080_data,
                   int_083_012_data, int_083_015_data, int_086_1_data, int_086_2_data, int_088_data, int_091_data) #missing LU066
  names(int_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS015","LU080",
                       "LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
  
  unstim_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_int <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Separate per condition
  for(i in seq_along(int_list)){
    data <- int_list[[i]]
    
    
    if(grepl("TLS", names(int_list)[[i]])){ #if a tumour was spread over multiple experiments
      if(!grepl("LU035|LU052", names(int_list)[[i]])){ #and if that tumour was not present in the same legendplex (bc those already have the TLS number in front)
        number <- strsplit(names(int_list)[[i]], "_")
        number <- number[[1]][[2]]
        rownames(data) <- paste0(number, " ", rownames(data)) #add the experiment number in front of the rownames
      }
    }
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    
    if(nrow(unstim)!=0){
      unstim_average_int <- rbind(unstim_average_int, unstim)
      colnames(unstim_average_int) <- colnames(data)
    }
    
    if(nrow(aPD1)!=0){
      aPD1_average_int <- rbind(aPD1_average_int, aPD1)
      colnames(aPD1_average_int) <- colnames(data)
    }
    
    if(nrow(aPD1_aIFNR)!=0){
      aPD1_aIFNR_average_int <- rbind(aPD1_aIFNR_average_int, aPD1_aIFNR)
      colnames(aPD1_aIFNR_average_int) <- colnames(data)
    }
    
    if(nrow(hrIFNy)!=0){
      hrIFNy_average_int <- rbind(hrIFNy_average_int, hrIFNy)
      colnames(hrIFNy_average_int) <- colnames(data)
    }
    
  }
  
  #for tumour bed
  tumour_list <- list(tumour_023_013_data,tumour_023_014_data, tumour_035_013_data, tumour_035_015_data, tumour_052_010_data, tumour_052_014_data,tumour_066_012_data, tumour_066_015_data,
                      tumour_080_data,tumour_083_012_data, tumour_083_015_data, tumour_086_1_data, tumour_086_2_data, tumour_088_data, tumour_091_data) #missing LU066
  names(tumour_list) <- c("LU023_TLS013","LU023_TLS014","LU035_TLS013","LU035_TLS015","LU052_TLS010","LU052_TLS014","LU066_TLS012", "LU066_TLS015","LU080",
                       "LU083_TLS012", "LU083_TLS015","LU086-1", "LU086-2", "LU088", "LU091")  
  
  unstim_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  aPD1_aIFNR_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  hrIFNy_average_tumour <- data.frame(matrix(nrow = 0, ncol = 25))
  
  #Separate per condition. #NEED TO ADD EXPERIMENT NUMBER IN THERE SOMEHOW BC IT IS ADJUSTING FRAGMENT NAMING FOR ME WITH DUPLICATE ROWNAMES...
  for(i in seq_along(tumour_list)){
    data <- tumour_list[[i]]
    
    if(grepl("TLS", names(tumour_list)[[i]])){ #if a tumour was spread over multiple experiments
      if(!grepl("LU035|LU052", names(tumour_list)[[i]])){ #and if that tumour was not present in the same legendplex (bc those already have the TLS number in front)
        number <- strsplit(names(tumour_list)[[i]], "_")
        number <- number[[1]][[2]]
        rownames(data) <- paste0(number, " ", rownames(data)) #add the experiment number in front of the rownames
      }
    }
    
    unstim <- data[grepl("unstim", rownames(data), ignore.case = T),]
    aPD1 <- data[grepl("aPD1 fragment", rownames(data)),] #the "fragment" is there to prevent grabbing the combi condition
    aPD1_aIFNR <- data[grepl("aPD1\\+aIFNyR1", rownames(data)),]
    hrIFNy <- data[grepl("hrIFNy", rownames(data)),]
    
    
    if(nrow(unstim)!=0){
      unstim_average_tumour <- rbind(unstim_average_tumour, unstim)
      colnames(unstim_average_tumour) <- colnames(data)
    }
    
    if(nrow(aPD1)!=0){
      aPD1_average_tumour <- rbind(aPD1_average_tumour, aPD1)
      colnames(aPD1_average_tumour) <- colnames(data)
    }
    
    if(nrow(aPD1_aIFNR)!=0){
      aPD1_aIFNR_average_tumour <- rbind(aPD1_aIFNR_average_tumour, aPD1_aIFNR)
      colnames(aPD1_aIFNR_average_tumour) <- colnames(data)
    }
    
    if(nrow(hrIFNy)!=0){
      hrIFNy_average_tumour <- rbind(hrIFNy_average_tumour, hrIFNy)
      colnames(hrIFNy_average_tumour) <- colnames(data)
    }
    
  }
  
  
#generate cumulative scores (unstim, all 3 segments)
  
#to separate between TLS from one tumour across two experiments, add the experiment number to the names
bind <- list(TLS_023_013,TLS_023_014, TLS_035_013, TLS_035_015, TLS_052_010, TLS_052_014, TLS_066_012, TLS_066_015, TLS_083_012, TLS_083_015, TLS_088, TLS_091,
             int_023_013,int_023_014, int_035_013, int_035_015, int_052_010, int_052_014, int_066_015, int_080, int_083_012, int_083_015, int_086_1, int_086_2, int_088, int_091,
             tumour_023_013,tumour_023_014, tumour_035_013, tumour_035_015, tumour_052_010, tumour_052_014, tumour_066_012, tumour_066_015, tumour_080, tumour_083_012, tumour_083_015, tumour_086_1, tumour_086_2, tumour_088, tumour_091)

names(bind) <- c("TLS_023_013","TLS_023_014", "TLS_035_013", "TLS_035_015", "TLS_052_010", "TLS_052_014", "TLS_066_012", "TLS_066_015", "TLS_083_012", "TLS_083_015", "TLS_088", "TLS_091",
                 "int_023_013","int_023_014", "int_035_013", "int_035_015", "int_052_010", "int_052_014", "int_066_015", "int_080", "int_083_012", "int_083_015", "int_086_1", "int_086_2", "int_088", "int_091",
                 "tumour_023_013","tumour_023_014", "tumour_035_013", "tumour_035_015", "tumour_052_010", "tumour_052_014", "tumour_066_012", "tumour_066_015", "tumour_080", "tumour_083_012", "tumour_083_015", "tumour_086_1", "tumour_086_2", "tumour_088", "tumour_091")

for(i in seq_along(bind)){
    if(!grepl("035|052|080|086_1|086_2|088|091", names(bind)[[i]])){ #if a tumour was not present in the same legendplex (bc those already have the TLS number in front) or was not repeated, skip
      number <- strsplit(names(bind)[[i]], "_") #the TLS number is the second number in the name
      number <- number[[1]][[3]]
      bind[[i]] <- paste0("TLS", number, " ", bind[[i]]) #add the experiment number in front of the rownames
      
  }
}

TLS_new <- unlist(bind[c(1:12)])
int_new <- unlist(bind[c(13:26)])
tumour_new <- unlist(bind[c(27:41)])

unstim_segment <- rbind(unstim_average_TLS,unstim_average_int,unstim_average_tumour)
colSD <- apply(unstim_segment, 2, sd, na.rm = T)
colMean <- colMeans(unstim_segment, na.rm = T)
z_scores <- sweep(unstim_segment, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

cumulative_score <- as.data.frame(colSums(z_scores)) #calculate cumulative scores by adding up z-scores per fragment
rownames(cumulative_score) <- colnames(z_scores)
colnames(cumulative_score) <- "cumulative_score_fragment"
cumulative_score$segment <- NA
cumulative_score$segment[rownames(cumulative_score) %in% TLS_new] <- "TLS" #assign segments back to fragments based on whether they were in the segment-separated data frames
cumulative_score$segment[rownames(cumulative_score) %in% int_new] <- "Intermediate"
cumulative_score$segment[rownames(cumulative_score) %in% tumour_new] <- "Tumour bed"
tumours <- str_extract_part(rownames(cumulative_score), before = FALSE, pattern = "LU") #grab the characters that come after "LU"
tumours <- strsplit(tumours, " ") #grab only the first characters that come after "LU" (which is the tumour number)
tumours <- sapply(tumours, "[[", 1) 
tumours <- paste0("LU", tumours) #put the "LU" back in front of the tumour number
cumulative_score$tumour <- tumours
cumulative_score$segment <- factor(cumulative_score$segment, levels = c("TLS", "Intermediate", "Tumour bed"))

ggplot(cumulative_score, aes(x=segment, y= cumulative_score_fragment, col = tumour)) + geom_point() + ggtitle("unstim") + xlab("segment") + ylab("cumulative score") + theme_bw()

#export to prism
colnames(cumulative_score) <- c("value", "segment", "tumour")
cumulative_score$name <- rownames(cumulative_score)
test <- cumulative_score #since reshape does not like to have an additional column in there and cannot work with duplicated tumour names: remove it
test$tumour <- NULL
wide <- reshape(test, idvar = "name", timevar = "segment", direction = "wide") #convert to wide
wide$name <- NULL
colnames(wide) <- gsub("value\\.", "", colnames(wide))
#since the row order is the same as the row order in the original data frame, we can just add the tumour name back in:
wide$tumours <- cumulative_score$tumour
#re-order
wide_unstim <- wide[,c(4,1,2,3)]

#export to excel
write.xlsx(wide_unstim, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_unstim_segments.xlsx",rowNames=T, colNames = T)


#repeat for aPD1
aPD1_segment <- rbind(aPD1_average_TLS,aPD1_average_int,aPD1_average_tumour)
colSD <- apply(aPD1_segment, 2, sd, na.rm = T)
colMean <- colMeans(aPD1_segment, na.rm = T)
z_scores <- sweep(aPD1_segment, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

cumulative_score <- as.data.frame(colSums(z_scores)) #calculate cumulative scores by adding up z-scores per fragment
rownames(cumulative_score) <- colnames(z_scores)
colnames(cumulative_score) <- "cumulative_score_fragment"
cumulative_score$segment <- NA
cumulative_score$segment[rownames(cumulative_score) %in% TLS_new] <- "TLS" #assign segments back to fragments based on whether they were in the segment-separated data frames
cumulative_score$segment[rownames(cumulative_score) %in% int_new] <- "Intermediate"
cumulative_score$segment[rownames(cumulative_score) %in% tumour_new] <- "Tumour bed"
tumours <- str_extract_part(rownames(cumulative_score), before = FALSE, pattern = "LU") #grab the characters that come after "LU"
tumours <- strsplit(tumours, " ") #grab only the first characters that come after "LU" (which is the tumour number)
tumours <- sapply(tumours, "[[", 1) 
tumours <- paste0("LU", tumours) #put the "LU" back in front of the tumour number
cumulative_score$tumour <- tumours
cumulative_score$segment <- factor(cumulative_score$segment, levels = c("TLS", "Intermediate", "Tumour bed"))
colnames(cumulative_score) <- c("value", "segment", "tumour")
cumulative_score$name <- rownames(cumulative_score)
test <- cumulative_score #since reshape does not like to have an additional column in there and cannot work with duplicated tumour names: remove it
test$tumour <- NULL
wide <- reshape(test, idvar = "name", timevar = "segment", direction = "wide") #convert to wide
wide$name <- NULL
colnames(wide) <- gsub("value\\.", "", colnames(wide))
#since the row order is the same as the row order in the original data frame, we can just add the tumour name back in:
wide$tumours <- cumulative_score$tumour
#re-order
wide_PD1 <- wide[,c(4,1,2,3)]

#export to excel
write.xlsx(wide_PD1, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_aPD1_segments.xlsx",rowNames=T, colNames = T)

#per tumour
full <- rbind(unstim_segment, aPD1_segment)
colSD <- apply(full, 2, sd, na.rm = T)
colMean <- colMeans(full, na.rm = T)
z_scores <- sweep(full, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)
cumulative_score <- as.data.frame(colSums(z_scores)) #calculate cumulative scores by adding up z-scores per fragment
rownames(cumulative_score) <- colnames(z_scores)
colnames(cumulative_score) <- "cumulative_score_fragment"
tumours <- str_extract_part(rownames(cumulative_score), before = FALSE, pattern = "LU") #grab the characters that come after "LU"
tumours <- strsplit(tumours, " ") #grab only the first characters that come after "LU" (which is the tumour number)
tumours <- sapply(tumours, "[[", 1) 
tumours <- paste0("LU", tumours) #put the "LU" back in front of the tumour number
cumulative_score$tumour <- tumours
cumulative_score$condition <- NA
cumulative_score$condition[grepl("unstim", rownames(cumulative_score))] <- "unstim"
cumulative_score$condition[grepl("PD1", rownames(cumulative_score))] <- "aPD1"
cumulative_score$tumour[grepl("LU086", rownames(cumulative_score))] <- "LU086"

cumulative_score_tumour <- split(cumulative_score, cumulative_score$tumour)

for(i in seq_along(cumulative_score_tumour)){
  tumour <- names(cumulative_score_tumour)[[i]]
  cumulative_score_tumour[[i]]$tumour <- NULL
  cumulative_score_tumour[[i]]$name <- rownames(cumulative_score_tumour[[i]])
  
  wide <- reshape(cumulative_score_tumour[[i]], idvar = "name", timevar = "condition", direction = "wide") #convert to wide
  colnames(wide) <- c("fragment","unstim", "aPD1")
  
  write.csv(wide, paste0("~/Documents/Analyses/TLS/mass_data_export/", tumour, "_cumulative_score_unstim_aPD1.csv"))
  
}

#now compare unstim to aPD1 in the three different segments
#TLS
TLS_condition <- rbind(unstim_average_TLS,aPD1_average_TLS)
colSD <- apply(TLS_condition, 2, sd, na.rm = T)
colMean <- colMeans(TLS_condition, na.rm = T)
z_scores <- sweep(TLS_condition, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

cumulative_score <- as.data.frame(colSums(z_scores)) #calculate cumulative scores by adding up z-scores per fragment
rownames(cumulative_score) <- colnames(z_scores)
colnames(cumulative_score) <- "cumulative_score_fragment"
cumulative_score$condition <- NULL
cumulative_score$condition[rownames(cumulative_score) %in% rownames(unstim_average_TLS)] <- "unstim" #assign segments back to fragments based on whether they were in the segment-separated data frames
cumulative_score$condition[rownames(cumulative_score) %in% rownames(aPD1_average_TLS)] <- "aPD1"
tumours <- str_extract_part(rownames(cumulative_score), before = FALSE, pattern = "LU") #grab the characters that come after "LU"
tumours <- strsplit(tumours, " ") #grab only the first characters that come after "LU" (which is the tumour number)
tumours <- sapply(tumours, "[[", 1) 
tumours <- paste0("LU", tumours) #put the "LU" back in front of the tumour number
cumulative_score$tumour <- tumours
cumulative_score$condition <- factor(cumulative_score$condition, levels = c("unstim", "aPD1"))
colnames(cumulative_score) <- c("value", "condition", "tumour")
cumulative_score$name <- rownames(cumulative_score)
test <- cumulative_score #since reshape does not like to have an additional column in there and cannot work with duplicated tumour names: remove it
test$tumour <- NULL
wide <- reshape(test, idvar = "name", timevar = "condition", direction = "wide") #convert to wide
wide$name <- NULL
colnames(wide) <- gsub("value\\.", "", colnames(wide))
#since the row order is the same as the row order in the original data frame, we can just add the tumour name back in:
wide$tumours <- cumulative_score$tumour
#re-order
wide_TLS <- wide[,c(3,1,2)]

#export to excel
write.xlsx(wide_TLS, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_TLS_conditions.xlsx",rowNames=T, colNames = T)

#Intermediate
int_condition <- rbind(unstim_average_int,aPD1_average_int)
colSD <- apply(int_condition, 2, sd, na.rm = T)
colMean <- colMeans(int_condition, na.rm = T)
z_scores <- sweep(int_condition, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

cumulative_score <- as.data.frame(colSums(z_scores)) #calculate cumulative scores by adding up z-scores per fragment
rownames(cumulative_score) <- colnames(z_scores)
colnames(cumulative_score) <- "cumulative_score_fragment"
cumulative_score$condition <- NULL
cumulative_score$condition[rownames(cumulative_score) %in% rownames(unstim_average_int)] <- "unstim" #assign segments back to fragments based on whether they were in the segment-separated data frames
cumulative_score$condition[rownames(cumulative_score) %in% rownames(aPD1_average_int)] <- "aPD1"
tumours <- str_extract_part(rownames(cumulative_score), before = FALSE, pattern = "LU") #grab the characters that come after "LU"
tumours <- strsplit(tumours, " ") #grab only the first characters that come after "LU" (which is the tumour number)
tumours <- sapply(tumours, "[[", 1) 
tumours <- paste0("LU", tumours) #put the "LU" back in front of the tumour number
cumulative_score$tumour <- tumours
cumulative_score$condition <- factor(cumulative_score$condition, levels = c("unstim", "aPD1"))
colnames(cumulative_score) <- c("value", "condition", "tumour")
cumulative_score$name <- rownames(cumulative_score)
test <- cumulative_score #since reshape does not like to have an additional column in there and cannot work with duplicated tumour names: remove it
test$tumour <- NULL
wide <- reshape(test, idvar = "name", timevar = "condition", direction = "wide") #convert to wide
wide$name <- NULL
colnames(wide) <- gsub("value\\.", "", colnames(wide))
#since the row order is the same as the row order in the original data frame, we can just add the tumour name back in:
wide$tumours <- cumulative_score$tumour
#re-order
wide_int <- wide[,c(3,1,2)]

#export to excel
write.xlsx(wide_int, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_int_conditions.xlsx",rowNames=T, colNames = T)


#tumour bed
tumour_condition <- rbind(unstim_average_tumour,aPD1_average_tumour)
colSD <- apply(tumour_condition, 2, sd, na.rm = T)
colMean <- colMeans(tumour_condition, na.rm = T)
z_scores <- sweep(tumour_condition, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)

cumulative_score <- as.data.frame(colSums(z_scores)) #calculate cumulative scores by adding up z-scores per fragment
rownames(cumulative_score) <- colnames(z_scores)
colnames(cumulative_score) <- "cumulative_score_fragment"
cumulative_score$condition <- NULL
cumulative_score$condition[rownames(cumulative_score) %in% rownames(unstim_average_tumour)] <- "unstim" #assign segments back to fragments based on whether they were in the segment-separated data frames
cumulative_score$condition[rownames(cumulative_score) %in% rownames(aPD1_average_tumour)] <- "aPD1"
tumours <- str_extract_part(rownames(cumulative_score), before = FALSE, pattern = "LU") #grab the characters that come after "LU"
tumours <- strsplit(tumours, " ") #grab only the first characters that come after "LU" (which is the tumour number)
tumours <- sapply(tumours, "[[", 1) 
tumours <- paste0("LU", tumours) #put the "LU" back in front of the tumour number
cumulative_score$tumour <- tumours
cumulative_score$condition <- factor(cumulative_score$condition, levels = c("unstim", "aPD1"))
colnames(cumulative_score) <- c("value", "condition", "tumour")
cumulative_score$name <- rownames(cumulative_score)
test <- cumulative_score #since reshape does not like to have an additional column in there and cannot work with duplicated tumour names: remove it
test$tumour <- NULL
wide <- reshape(test, idvar = "name", timevar = "condition", direction = "wide") #convert to wide
wide$name <- NULL
colnames(wide) <- gsub("value\\.", "", colnames(wide))
#since the row order is the same as the row order in the original data frame, we can just add the tumour name back in:
wide$tumours <- cumulative_score$tumour
#re-order
wide_tumour <- wide[,c(3,1,2)]

#export to excel
write.xlsx(wide_tumour, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_tumour_conditions.xlsx",rowNames=T, colNames = T)


#repeat, but only take proven responders
wide_unstim_R <- wide_unstim[wide_unstim$tumours %in% c("LU023","LU052","LU083"),]
wide_PD1_R <- wide_PD1[wide_PD1$tumours %in% c("LU023","LU052","LU083"),]
wide_TLS_R <- wide_TLS[wide_TLS$tumours %in% c("LU023","LU052","LU083"),]
wide_int_R <- wide_int[wide_int$tumours %in% c("LU023","LU052","LU083"),]
wide_tumour_R <- wide_tumour[wide_tumour$tumours %in% c("LU023","LU052","LU083"),]

write.xlsx(wide_unstim_R, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_unstim_responder_segments.xlsx",rowNames=T, colNames = T)
write.xlsx(wide_PD1_R, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_aPD1_responder_segments.xlsx",rowNames=T, colNames = T)
write.xlsx(wide_TLS_R, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_TLS_responder_conditions.xlsx",rowNames=T, colNames = T)
write.xlsx(wide_int_R, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_int_responder_conditions.xlsx",rowNames=T, colNames = T)
write.xlsx(wide_tumour_R, "~/Documents/Analyses/TLS/final_all_fragments/cumulative_score_tumour_responder_conditions.xlsx",rowNames=T, colNames = T)

#plot cytokine changes for responders
cytokines_unstim <- rbind(unstim_average_TLS, unstim_average_int, unstim_average_tumour) #grab unstim cytokine values
cytokines_aPD1 <- rbind(aPD1_average_TLS, aPD1_average_int, aPD1_average_tumour) #grab aPD1 cytokine values

cytokines_unstim <- cytokines_unstim[grepl("LU083|LU052|LU023", rownames(cytokines_unstim)),] #only take responders
cytokines_aPD1 <- cytokines_aPD1[grepl("LU083|LU052|LU023", rownames(cytokines_aPD1)),]

cytokines_unstim <- as.data.frame(colMeans(cytokines_unstim, na.rm = T)) #calculate averages per cytokine
cytokines_aPD1 <- as.data.frame(colMeans(cytokines_aPD1, na.rm = T)) #calculate averages per cytokine

colnames(cytokines_unstim) <- "unstim"
colnames(cytokines_aPD1) <- "aPD1"

cytokines <- cbind(cytokines_unstim, cytokines_aPD1)
cytokines$log2fc <- log2(cytokines_aPD1/cytokines_unstim)
colnames(cytokines) <- c("unstim","aPD1","log2fc")
cytokines[] <- lapply(cytokines, function(x) if(class(x) == 'data.frame') unlist(x) else x) #log2fc returns a column that is a data.frame somehow (class(cytokines[,3])) so change this
write.xlsx(cytokines, "~/Documents/Analyses/TLS/final_all_fragments/responder_average_cytokines_unstim_aPD1.xlsx",rowNames = T)

#bubble plot attempt
cytokines$mediator <- rownames(cytokines)
cytokines$position <- "x"
library(viridis)
ggplot(cytokines, aes(x=position, y = mediator, col = log2fc, size = unstim)) + geom_point() + scale_size_continuous(range=c(5,12)) +scale_color_viridis() + theme_bw() + xlab("")

#separate by segment (they all have at least 2 fragments)
unstim_average_TLS_small <- unstim_average_TLS[grepl("LU083|LU052|LU023", rownames(unstim_average_TLS)),] #only take responders
aPD1_average_TLS_small <- aPD1_average_TLS[grepl("LU083|LU052|LU023", rownames(aPD1_average_TLS)),] #only take responders

unstim_average_int_small <- unstim_average_int[grepl("LU083|LU052|LU023", rownames(unstim_average_int)),] #only take responders
aPD1_average_int_small <- aPD1_average_int[grepl("LU083|LU052|LU023", rownames(aPD1_average_int)),] #only take responders

unstim_average_tumour_small <- unstim_average_tumour[grepl("LU083|LU052|LU023", rownames(unstim_average_tumour)),] #only take responders
aPD1_average_tumour_small <- aPD1_average_tumour[grepl("LU083|LU052|LU023", rownames(aPD1_average_tumour)),] #only take responders


unstim_average_TLS_small <- colMeans(unstim_average_TLS_small, na.rm = T) #calculate averages per cytokine
aPD1_average_TLS_small <- colMeans(aPD1_average_TLS_small, na.rm = T) #calculate averages per cytokine
unstim_average_int_small <- colMeans(unstim_average_int_small, na.rm = T) #calculate averages per cytokine
aPD1_average_int_small <- colMeans(aPD1_average_int_small, na.rm = T) #calculate averages per cytokine
unstim_average_tumour_small <- colMeans(unstim_average_tumour_small, na.rm = T) #calculate averages per cytokine
aPD1_average_tumour_small <- colMeans(aPD1_average_tumour_small, na.rm = T) #calculate averages per cytokine

TLS_log2FC <- log2(aPD1_average_TLS_small/unstim_average_TLS_small)
int_log2FC <- log2(aPD1_average_int_small/unstim_average_int_small)
tumour_log2FC <- log2(aPD1_average_tumour_small/unstim_average_tumour_small)

cytokines <- data.frame(cbind(TLS_log2FC, int_log2FC, tumour_log2FC))
colnames(cytokines) <- c("TLS","Intermediate","Tumour bed")
write.xlsx(cytokines, "~/Documents/Analyses/TLS/final_all_fragments/responder_log2FC_segment.xlsx",rowNames = T)


#show heterogeneity in LU083 and LU052 individual heatmaps
lu083 <- rbind(lu083_TLS012_copy, lu083_TLS015_copy) #use the copy versions as they also have experiment name in their fragment names: prevents duplicate fragment names and renaming
lu083 <- lu083[!grepl("average|hrIFNy|RN440-15|RN440-16|aPD1\\+aIFNyR1", rownames(lu083)),]
colSD <- apply(lu083, 2, sd, na.rm = T)
colMean <- colMeans(lu083, na.rm = T)
z_scores <- sweep(lu083, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)
z_scores[z_scores == "NaN"] <- 0 #column SD is zero in some cases, likely due to normalisation

#rename vectors with fragment names for TLS etc. to match experiment tags
TLS_083_012_copy <- gsub("LU083", "LU083_TLS012", TLS_083_012)
TLS_083_015_copy <- gsub("LU083", "LU083_TLS015", TLS_083_015)

int_083_012_copy <- gsub("LU083", "LU083_TLS012", int_083_012)
int_083_015_copy <- gsub("LU083", "LU083_TLS015", int_083_015)

tumour_083_012_copy <- gsub("LU083", "LU083_TLS012", tumour_083_012)
tumour_083_015_copy <- gsub("LU083", "LU083_TLS015", tumour_083_015)

#reorder columns in this order: TLS-Intermediate-Tumour bed
TLS_083 <- z_scores[,colnames(z_scores) %in% c(TLS_083_012_copy, TLS_083_015_copy)]
int_083 <- z_scores[,colnames(z_scores) %in% c(int_083_012_copy, int_083_015_copy)]
tumour_083 <- z_scores[,colnames(z_scores) %in% c(tumour_083_012_copy, tumour_083_015_copy)]

z_scores <- as.data.frame(cbind(TLS_083, int_083, tumour_083))

heatmap <- pheatmap(z_scores)

row_order <- heatmap$tree_row$order

z_scores <- z_scores[order(row_order),]

write.xlsx(z_scores, "~/Documents/Analyses/TLS/final_all_fragments/lu083_z_scores.xlsx", rowNames= T)

#for LU052
lu052 <- rbind(lu052_TLS010, lu052_TLS014) #no need to add experiment names; they are already there
lu052 <- lu052[!grepl("average|hrIFNy|aPD1\\+aIFNyR1", rownames(lu052)),]
colSD <- apply(lu052, 2, sd, na.rm = T)
colMean <- colMeans(lu052, na.rm = T)
z_scores <- sweep(lu052, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)
z_scores[z_scores == "NaN"] <- 0 #column SD is zero in some cases, likely due to normalisation

#reorder columns in this order: TLS-Intermediate-Tumour bed
TLS_052 <- z_scores[,colnames(z_scores) %in% c(TLS_052_010, TLS_052_014)]
int_052 <- z_scores[,colnames(z_scores) %in% c(int_052_010, int_052_014)]
tumour_052 <- z_scores[,colnames(z_scores) %in% c(tumour_052_010, tumour_052_014)]

z_scores <- as.data.frame(cbind(TLS_052, int_052, tumour_052))
heatmap <- pheatmap(z_scores)

row_order <- heatmap$tree_row$order
z_scores <- z_scores[order(row_order),]

write.xlsx(z_scores, "~/Documents/Analyses/TLS/final_all_fragments/lu052_z_scores.xlsx", rowNames= T)


#for LU023
#show heterogeneity in LU083 and LU052 individual heatmaps
lu023 <- rbind(lu023_TLS013_copy, lu023_TLS014_copy) #use the copy versions as they also have experiment name in their fragment names: prevents duplicate fragment names and renaming
lu023 <- lu023[!grepl("average|hrIFNy|RN440-15|RN440-16|aPD1\\+aIFNyR1", rownames(lu023)),]
colSD <- apply(lu023, 2, sd, na.rm = T)
colMean <- colMeans(lu023, na.rm = T)
z_scores <- sweep(lu023, 2, colMean, "-") #subtract means per column
z_scores <- sweep(z_scores, 2, colSD, "/") #divide by standard deviation per column
z_scores <- t(z_scores)
z_scores[z_scores == "NaN"] <- 0 #column SD is zero in some cases, likely due to normalisation

#rename vectors with fragment names for TLS etc. to match experiment tags
TLS_023_013_copy <- gsub("LU023", "LU023_TLS013", TLS_023_013)
TLS_023_014_copy <- gsub("LU023", "LU023_TLS014", TLS_023_014)

int_023_013_copy <- gsub("LU023", "LU023_TLS013", int_023_013)
int_023_014_copy <- gsub("LU023", "LU023_TLS014", int_023_014)

tumour_023_013_copy <- gsub("LU023", "LU023_TLS013", tumour_023_013)
tumour_023_014_copy <- gsub("LU023", "LU023_TLS014", tumour_023_014)

#reorder columns in this order: TLS-Intermediate-Tumour bed
TLS_023 <- z_scores[,colnames(z_scores) %in% c(TLS_023_013_copy, TLS_023_014_copy)]
int_023 <- z_scores[,colnames(z_scores) %in% c(int_023_013_copy, int_023_014_copy)]
tumour_023 <- z_scores[,colnames(z_scores) %in% c(tumour_023_013_copy, tumour_023_014_copy)]

z_scores <- as.data.frame(cbind(TLS_023, int_023, tumour_023))

heatmap <- pheatmap(z_scores)

row_order <- heatmap$tree_row$order

z_scores <- z_scores[order(row_order),]

write.xlsx(z_scores, "~/Documents/Analyses/TLS/final_all_fragments/lu023_z_scores.xlsx", rowNames= T)
