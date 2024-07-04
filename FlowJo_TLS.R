#Load libraries
library(readxl)
library(openxlsx)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(ggsignif)

#Only run the below section if you want to regenerate the file
#---------------------
#read in data and add experiment number
flow_file_TLS016 <- read_excel("~/Documents/Experiments/TLS/TLS016_LU088_LU091_22042024_SZ_NS/TLS016.xls")
flow_file_TLS016$experiment <- "TLS016"
flow_file_TLS015 <- read_excel("~/Documents/Experiments/TLS/TLS015_LU066_LU083_15042024_SZ_NS/TLS015.xls")
flow_file_TLS015$experiment <- "TLS015"
flow_file_TLS014 <- read_excel("~/Documents/Experiments/TLS/TLS014_LU023_LU052_10042024_SZ_NS/TLS014.xls")
flow_file_TLS014$experiment <- "TLS014"
flow_file_TLS013 <- read_excel("~/Documents/Experiments/TLS/TLS013_LU023_LU035_28032024_SZ_NS/TLS013.xls")
flow_file_TLS013$experiment <- "TLS013"
flow_file_TLS012 <- read_excel("~/Documents/Experiments/TLS/TLS012_LU066_LU083_22032024_SZ_NS_PK_GC_RW/TLS012.xls")
flow_file_TLS012$experiment <- "TLS012"
flow_file_TLS011 <- read_excel("~/Documents/Experiments/TLS/TLS011_21032024_SZ_NS_PK_GC_RW/TLS011.xls")
flow_file_TLS011$experiment <- "TLS011"
flow_file_TLS010 <- read_excel("~/Documents/Experiments/TLS/TLS010 SZ NS 14032024/TLS010_AF.xls")
flow_file_TLS010$experiment <- "TLS010"
flow_file_TLS009 <- read_excel("~/Documents/Experiments/TLS/TLS009_LU086_SZ_NS_29022024/TLS009.xls")
flow_file_TLS009$experiment <- "TLS009"


flow_file <- rbind(flow_file_TLS009, flow_file_TLS010, flow_file_TLS011, flow_file_TLS012,
                   flow_file_TLS013, flow_file_TLS014, flow_file_TLS015, flow_file_TLS016)
colnames(flow_file) <- c("Name", "Statistic", "Cells","Experiment")

#grab individual fragments, distinguished as #1, #2 etc.
fragments <- flow_file[grepl("#|fragment",flow_file$Name),]

#remove entries that do not contain full fragment or B cell data
fragments <- fragments[!grepl("CD3|T cell|Myeloid|CD20+|Cell Cycle|unknown",fragments$Name,ignore.case = T),]

#remove PBMC data
fragments <- fragments[!grepl("PBMC",fragments$Name,ignore.case = T),]

#remove single stains
fragments <- fragments[grepl("LU",fragments$Name,ignore.case = T),]

#The median of a marker is always below the population it refers to. Use this:
MFI_expression_file <- data.frame()

for(i in 1:nrow(flow_file)){
  if(grepl("Median", flow_file$Name[[i]])){
    
    #grab median expression value
    MFI_row_stat <- strsplit(flow_file$Name[[i]],"=")
    MFI_row_marker <- MFI_row_stat[[1]][1]
    MFI_row_marker <- gsub("Median : Comp-|-A| " , "", MFI_row_marker)
    
    MFI_row_stat <- MFI_row_stat[[1]][2]
    MFI_row_stat <- gsub(" ", "", MFI_row_stat) #remove spaces
    
    #find the row it corresponds to
    for(j in 1:i){
      #look in all rows above the entry
      if(grepl("LU", flow_file$Name[[i-j]])){
        row <- flow_file$Name[[i-j]] #if there is a tumour sample found, couple the value to this sample and exit the loop
        break
      } else if (grepl("stim|unstim", flow_file$Name[[i-j]])) {
        row <- NA #if there is a PBMC sample found, do not enter a sample and exit the loop
        break
      } else {
        next #if no sample is found on this row, go to the row above it
      }
    } 
    
    #quick check: if no B cells are there, the MFI will be n/a, so leave this entry out
    if(is.na(MFI_row_stat)){
      next
    }
    
    if(is.na(row)){
      next #if the row is NA, this was a PBMC sample, so do not couple it
    }
    
    #make the data entry
    row <- gsub("/B cells", "", row)
    row_full <- data.frame(Name = row, marker = MFI_row_marker, MFI = MFI_row_stat)
    
    #append to data.frame
    MFI_expression_file <- rbind(MFI_expression_file, row_full)
  }
}


#remove NA values
MFI_expression_file <- MFI_expression_file[!is.na(MFI_expression_file$Name),]

#Couple the fluorochromes to an actual marker
MFI_expression_file$marker[MFI_expression_file$marker == "BB515"] <- "HLA-DR"
MFI_expression_file$marker[MFI_expression_file$marker == "BUV563"] <- "CD86"
MFI_expression_file$marker[MFI_expression_file$marker == "PE"] <- "CD80"

#reshaping
MFI_expression_file_wide <- spread(MFI_expression_file, key = marker, value = MFI)
MFI_expression_file_wide[,-1] <- sapply(MFI_expression_file_wide[,-1],as.numeric)

#add to table
fragments <- left_join(fragments, MFI_expression_file_wide, by = "Name")

#add B cell information present in second row
fragments_noB <- fragments[!grepl("B cell", fragments$Name),]
B_cell_info <- data.frame()

for (i in 1:nrow(fragments_noB)){
  name <- fragments_noB$Name[[i]]
  name_B_cell <- paste0(name, "/B cells") #grab flowjo path of B cell gate
  B_cell <- fragments[fragments$Name %in% name_B_cell,]
  
  colnames(B_cell)[colnames((B_cell)) == "Cells"] <- "B_cells"
  colnames(B_cell)[colnames((B_cell)) == "Statistic"] <- "B_cell_percentage"
  
  B_cell$Name <- gsub("/B cells", "", B_cell$Name)
  
  B_cell$CD80 <- NULL
  B_cell$CD86 <- NULL
  B_cell$`HLA-DR` <- NULL
  
  B_cell_info <- rbind(B_cell_info,B_cell)
  
}

#add B cell information to data frame
B_cell_info$Experiment <- NULL
B_cell_info[,-1] <- sapply(B_cell_info[,-1],as.numeric)
fragments_noB <- left_join(fragments_noB, B_cell_info, by = "Name")

#note: NA values for B cell counts/percentages is because I did not gate fragments split up over multiple clusters (like 1 and 1_2)

#rename columns for clarity
colnames(fragments_noB) <- c("Name","percent_of_CD45","cells", "experiment", "CD80_B_cell_MFI","CD86_B_cell_MFI","HLA-DR_B_cell_MFI", "B_cell_percentage",
                             "B_cells")
fragments_noB$percent_of_CD45 <- as.numeric(fragments_noB$percent_of_CD45)


#add T cell information
T_cell_data <- flow_file[grepl("CD3\\+ T", flow_file$Name),]
T_cell_data <- T_cell_data[!grepl("PBMC",T_cell_data$Name),]
test <- strsplit(T_cell_data$Name, "/")

T_cell_info <- data.frame()

for (i in seq_along(test)){
  result <- test[[i]][grepl("\\+",test[[i]])] #specific population on this line of the file
  
  sample <- test[[i]][grepl("\\.fcs",test[[i]])] #sample in which the population was found
  
  if(grepl("unstim|stim",sample)){
    if(!grepl("LU",sample)){ #if it is a PBMC sample, skip
      next
    }
  }
  
  #give the T cell population a name:
  if(any(grepl("CD4",result))){
    T_cell_name <- paste0("CD4+ ", tail(result,1))
  } else if(any(grepl("CD8",result))){
    T_cell_name <- paste0("CD8+ ", tail(result,1))
  } else if(any(grepl("DN", test[[i]]))){ #check DN and DP populations in original list as they do not contain a + sign
    T_cell_name <- "CD3+ DN"
  } else if(any(grepl("DP", test[[i]]))){
    T_cell_name <- "CD3+ DP"
  } else {
    T_cell_name <- "CD3+"
  }
  
  #grab the number of cells and percentages for that population
  T_cell_count <- T_cell_data[i,"Cells"] #since the list and the dataframe are in the same order, we can just work with the list index for the data frame
  T_cell_count <- as.numeric(T_cell_count)
  T_cell_percentage <- T_cell_data[i,"Statistic"]
  T_cell_percentage <- as.numeric(T_cell_percentage)
  
  #add the real sample name (without all the T cell gating)
  sample_name <- T_cell_data[i,"Name"]
  sample_name <- gsub("/CD3\\+ T cells", ":",sample_name)
  sample_name <- strsplit(sample_name, ":")
  sample_name <- sample_name[[1]][1]
  
  #make and export the data
  T_cell_row <- data.frame(Name = sample_name, population = T_cell_name, percentage = T_cell_percentage, count = T_cell_count)
  
  T_cell_info <- rbind(T_cell_info, T_cell_row)
}

#filter out some artefacts
T_cell_info$population[T_cell_info$population == "CD4+ CD4+ T"] <- "CD4+"
T_cell_info$population[T_cell_info$population == "CD4+ CD4+"] <- "CD4+"
T_cell_info$population[T_cell_info$population == "CD8+ CD8+ T"] <- "CD8+"
T_cell_info$population[T_cell_info$population == "CD8+ CD8+"] <- "CD8+"

#filter out any percentages that exceed 100 (this might be a FlowJo artefact?)
T_cell_info <- T_cell_info[T_cell_info$percentage <= 100,] #does not filter anything out anymore so artefacts are gone

#convert to wide, separately for percentage and count
T_cell_percentages <- T_cell_info[,c("Name","population","percentage")]
T_cell_counts <- T_cell_info[,c("Name","population","count")]

T_cell_percentages <- spread(T_cell_percentages, key = population, value = percentage)
colnames(T_cell_percentages)[-1] <- paste0(colnames(T_cell_percentages)[-1], "_percentage")
T_cell_counts <- spread(T_cell_counts, key = population, value = count)
colnames(T_cell_counts)[-1] <- paste0(colnames(T_cell_counts)[-1], "_counts")

#get the condition and tumour out by strsplit
fragments_noB_condition <- strsplit(fragments_noB$Name, "/")
fragments_noB_condition <- sapply(fragments_noB_condition,"[[",1)

#few gsub commands to get all condition naming the same
fragments_noB_condition <- gsub("_"," ",fragments_noB_condition,ignore.case = T)
fragments_noB_condition <- gsub("fragment","",fragments_noB_condition,ignore.case = T)
fragments_noB_condition <- gsub("pool","",fragments_noB_condition,ignore.case = T)
fragments_noB_condition <- gsub("\\.fcs", "", fragments_noB_condition)
fragments_noB_condition <- gsub("Tumor|Tumor ", "", fragments_noB_condition)
fragments_noB_condition <- gsub("aPD1\\+ahrIFNyR1|aIFNgR1\\+aPD1|aIFNgR \\+ aPD1|aIFNgR1 \\+ aPD1|aPD1\\+aIFNR1|aIFNgR1\\+aPD1|aIFNgR \\+ aPD1|ahrIFNyR1\\+aPD1|aPD1\\+aIFNR|aPD1\\+aIFNR|hrIFNyR1 \\+ aPD1|IFNgR1 \\+ aPD1", "aIFNyR1+aPD1", fragments_noB_condition)
fragments_noB_condition <- gsub("rhIFNg|IFNg", "hrIFNy", fragments_noB_condition)
fragments_noB_condition <- gsub("rhIFNg|IFNg", "hrIFNy", fragments_noB_condition)
fragments_noB_condition <- gsub("hrhrIFNy", "hrIFNy", fragments_noB_condition)
fragments_noB_condition <- gsub("hrhrIFNy", "hrIFNy", fragments_noB_condition) #artefact caused by the line above
fragments_noB_condition <- gsub("Unstim", "unstim", fragments_noB_condition)
fragments_noB_condition <- gsub(" 1| 2| 3| 4| 5| 6| 7| 8| 9| 10| 11", "", fragments_noB_condition)
fragments_noB_condition <- gsub("LU86", "LU086", fragments_noB_condition)
fragments_noB_condition <- gsub("unstim ", "unstim", fragments_noB_condition)

#Split the tumour, number and condition
test <- strsplit(fragments_noB_condition, " ")

#grab tumour numbers out
tumours <- unlist(lapply(test, grep, pattern="LU", value=TRUE)) #works!

#grab the condition (which is always the last entry in the list)
conditions <- sapply(test,tail,1)

conditions[conditions == "aPD1+ahrIFNyR1"] <- "aIFNyR1+aPD1"


#add the information back to the sample list
fragments_noB$tumour <- tumours
fragments_noB$condition <- conditions

#add fragment information
names <- strsplit(fragments_noB$Name, "#|fragment_")
names <- sapply(names,"[[",2)
names <- gsub("\\?","",names)
names <- gsub("_1.fcs|.fcs","",names)

fragments_noB$fragment <- names

#remove fragments that are split up into two (marked by a _2); their concatenated sample should be present already in the file
duplicate <- fragments_noB[grepl("_2", fragments_noB$fragment),]

concat_final <- data.frame()
delete <- data.frame()

for(i in 1:nrow(duplicate)){
  #grab all entries: the concatenated entry but also the individual entries
  row <- duplicate[i,]
  fragment <- gsub("_2|_3", "", row$fragment)
  
  #for entries concerning conditional replicates (aPD1 vs aPD1 2), the other fragment will also be identified, so fix this:
  if(grepl("2.fcs", row$Name)){
    row$condition <- paste0(row$condition, " 2")
  }
  
  corresponding <- fragments_noB[fragments_noB$experiment %in% row$experiment & fragments_noB$tumour %in% row$tumour & fragments_noB$condition %in% row$condition & fragments_noB$fragment %in% fragment,]
  
  if(any(grepl("2.fcs", corresponding$Name)) & !grepl("2.fcs", row$Name)){
    corresponding <- corresponding[!grepl("2.fcs", corresponding$Name),] #if a fragment was grabbed from a conditional replicate where the original duplicate fragment was
    #not from a replicate, discard this entry
  }
  
  #this is to identify a fragment ending on _3:
  corresponding3 <- fragments_noB[fragments_noB$experiment %in% row$experiment & fragments_noB$tumour %in% row$tumour & fragments_noB$condition %in% row$condition & fragments_noB$fragment %in% paste0(fragment,"_3"),]
  
  full_dup <- rbind(row,corresponding)
  
  if(nrow(corresponding3) != 0){
    full_dup <- rbind(full_dup,corresponding3)
  }
  
  #for entries concerning LU086, this will also identify fragments coming from the other replicate (LU086_1 and LU086_2).
  if(grepl("LU086", row$Name) & row$experiment == "TLS011"){
    if(grepl("LU086_1", row$Name)){
      full_dup <- full_dup[grepl("LU086_1", full_dup$Name),] #exclude LU086_2 entries
    } else if(grepl("LU086_2", row$Name)){
      full_dup <- full_dup[grepl("LU086_2", full_dup$Name),] #exclude LU086_2 entries
    }
  }
  
 
  #the trick is this: the concatenated entry containing both clusters will not have a "percent of CD45" statistic
  #as they are not gated from anything
  #Select the concatenated entry using this, calculate the percentage of CD45 stat by adding up the individual percentages
  concat <- full_dup[is.na(full_dup$percent_of_CD45),]
  individual <- full_dup[!is.na(full_dup$percent_of_CD45),]
  
  concat$percent_of_CD45 <- sum(individual$percent_of_CD45)
  
  #save the concatenated entry to later add this back into the full data frame
  concat_final <- rbind(concat_final, concat)
  
  #the individual entries and the unedited concatenated entry can now be removed from the full data frame:
  delete <- rbind(delete, full_dup)
  }

#now delete the individual and unfinished concat fragments
fragments_noB <- fragments_noB[!(fragments_noB$Name %in% delete$Name),]

#and add back the finished concat fragments
fragments_noB <- rbind(fragments_noB, concat_final)



#add information on whether the fragment has passed the 100 cell threshold
fragments_noB$cell_threshold_100 <- ifelse(fragments_noB$cells >= 100, "YES", "NO")

#add information on whether the fragment is a TLS, normal or tumour bed
fragments_noB$segment <- NA

for (i in 1:nrow(fragments_noB)){
  if (fragments_noB$B_cell_percentage[[i]] < 10 | is.na(fragments_noB$B_cell_percentage[[i]])){
    fragments_noB$segment[[i]] <- "tumour_bed"
  } else if (fragments_noB$B_cell_percentage[[i]] >= 10 & fragments_noB$B_cell_percentage[[i]] < 25){
    fragments_noB$segment[[i]] <- "intermediate"
  } else if (fragments_noB$B_cell_percentage[[i]] >= 25){
    fragments_noB$segment[[i]] <- "TLS"
  }
}


#separate samples with an extremely high IFN concentration
fragments_noB$condition[fragments_noB$condition == "hrIFNy" & fragments_noB$experiment %in% c("TLS009", "TLS010")] <- "hrIFNy_high"

#correction: for LU091, the last two columns (fragments #10 and #11) are actually LU088
fragments_noB$fragment[fragments_noB$tumour == "LU091" & fragments_noB$fragment == "10"] <- "12"
fragments_noB$fragment[fragments_noB$tumour == "LU091" & fragments_noB$fragment == "11"] <- "13"

fragments_noB$tumour[fragments_noB$tumour == "LU091" & fragments_noB$fragment %in% c("12","13")] <- "LU088"

#another correction: in TLS013, no hrIFNy was added for LU035, making this an unstim
fragments_noB$fragment[fragments_noB$tumour == "LU035" & fragments_noB$condition == "hrIFNy" & fragments_noB$experiment == "TLS013"] <- as.numeric(fragments_noB$fragment[fragments_noB$tumour == "LU035" & fragments_noB$condition == "hrIFNy" & fragments_noB$experiment == "TLS013"]) + 11
fragments_noB$condition[fragments_noB$tumour == "LU035" & fragments_noB$condition == "hrIFNy" & fragments_noB$experiment == "TLS013"] <- "unstim"

#last correction: sometimes, we fill multiple rows with the same condition, 
#making it such that there will be multiple fragments named 1 (and 2, and so on) in the same tumour in the same condition.
fragments_backup <- fragments_noB
fragments_noB$fragment[grepl("hrIFNy 2|aPD1 2|aPD1\\+aIFNgR1 2|aPD1\\+aIFNR 2",fragments_noB$Name)] <- as.numeric(fragments_noB$fragment[grepl("hrIFNy 2|aPD1 2|aPD1\\+aIFNgR1 2|aPD1\\+aIFNR 2",fragments_noB$Name)]) + 11

#add T cell information
T_cell_data_full <- full_join(T_cell_percentages,T_cell_counts, by = "Name")
fragments_noB <- left_join(fragments_noB, T_cell_data_full, by = "Name") #fixed the issue with T cell entries going missing in fragment duplication removal :)

#correction for an incorrectly named condition
fragments_noB$condition[fragments_noB$condition == "aPD1+ahrIFNyR1"] <- "aIFNyR1+aPD1"


#save the table
write.csv(fragments_noB, "~/Documents/Analyses/TLS/fragment_table.csv")


#--------------------
#read in table
fragments_noB <- read.csv("~/Documents/Analyses/TLS/fragment_table.csv", check.names = F)
fragments_noB[,1] <- NULL

#continue without RN samples
fragments_small <- fragments_noB[!(fragments_noB$condition %in% c("RN440-15","RN440-16")),]
#also remove fragments that did not pass the threshold
fragments_small <- fragments_small[fragments_small$cell_threshold_100 == "YES",]

#plot TLS number
counts <- fragments_noB %>%
  group_by(tumour, cell_threshold_100, condition, segment) %>%
  filter(cell_threshold_100 == "YES") %>%
  dplyr::count()

ggplot(counts,aes(x=factor(condition,levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high","RN440-15","RN440-16")),y=n,fill=segment)) + geom_col(position = "stack") + facet_wrap(~tumour, scales='free_x') +
  ylab("count") + xlab("condition") + scale_fill_manual(values = c("#CC9FFF","#99FFCC","#99CCFF")) +theme_minimal() + theme(text=element_text(size=20),axis.text.x = element_text(angle = 45, hjust = 1))

#export to prism
counts <- fragments_small %>%
  group_by(tumour, cell_threshold_100, condition, segment) %>%
  filter(cell_threshold_100 == "YES") %>%
  dplyr::count()

counts$cell_threshold_100 <- NULL

LU023 <- counts[counts$tumour == "LU023",]
LU035 <- counts[counts$tumour == "LU035",]
LU052 <- counts[counts$tumour == "LU052",]
LU066 <- counts[counts$tumour == "LU066",]
LU080 <- counts[counts$tumour == "LU080",]
LU083 <- counts[counts$tumour == "LU083",]
LU086 <- counts[counts$tumour == "LU086",]
LU088 <- counts[counts$tumour == "LU088",]
LU091 <- counts[counts$tumour == "LU091",]

count_list <- list(LU023, LU035, LU052, LU066, LU080, LU083, LU086, LU088, LU091)
names(count_list) <- c("LU023", "LU035", "LU052", "LU066", "LU080", "LU083", "LU086", "LU088", "LU091")

for (i in seq_along(count_list)){
  count_list[[i]]$tumour <- NULL
  count_list[[i]] <- as.data.frame(count_list[[i]])
  colnames(count_list[[i]]) <- c("condition", "segment", "value")
  wide <- reshape(count_list[[i]], idvar = "condition", timevar = "segment", direction = "wide")
  wide[is.na(wide)] <- 0
  
  
  conditions <- c("unstim", "aPD1", "aIFNyR1+aPD1", "hrIFNy")
  missing <- conditions[!(conditions %in% wide$condition)]
  if(length(missing) > 0){
    missing_add <- c(missing,0,0,0) #assuming no more than 1 condition missing
    suppressWarnings(wide <- rbind(wide, missing_add))
  }

  
  if(any(grepl("TLS",colnames(wide)))){
    colnames(wide)[[which(grepl("TLS",colnames(wide)))]] <- "TLS"
  } else{
    TLS <- c(0,0,0,0)
    wide$TLS <- TLS
  }
  
  if(any(grepl("intermediate",colnames(wide)))){
    colnames(wide)[[which(grepl("intermediate",colnames(wide)))]] <- "intermediate"
  } else{
    intermediate <- c(0,0,0,0)
    wide$intermediate <- intermediate
  }
  
  if(any(grepl("tumour_bed",colnames(wide)))){
    colnames(wide)[[which(grepl("tumour_bed",colnames(wide)))]] <- "tumour_bed"
  } else{
    tumour_bed <- c(0,0,0,0)
    wide$tumour_bed <- tumour_bed
  }

  wide <- wide[order(match(wide$condition, conditions)),]
  segments <- c("TLS","intermediate","tumour_bed")
  wide <- wide[,order(match(colnames(wide), segments))]
  
  count_list[[i]] <- wide
}

write.csv(count_list[["LU023"]], "~/Documents/Analyses/TLS/LU023_TLS_counts.csv")
write.csv(count_list[["LU035"]], "~/Documents/Analyses/TLS/LU035_TLS_counts.csv")
write.csv(count_list[["LU052"]], "~/Documents/Analyses/TLS/LU052_TLS_counts.csv")
write.csv(count_list[["LU066"]], "~/Documents/Analyses/TLS/LU066_TLS_counts.csv")
write.csv(count_list[["LU080"]], "~/Documents/Analyses/TLS/LU080_TLS_counts.csv")
write.csv(count_list[["LU083"]], "~/Documents/Analyses/TLS/LU083_TLS_counts.csv")
write.csv(count_list[["LU086"]], "~/Documents/Analyses/TLS/LU086_TLS_counts.csv")
write.csv(count_list[["LU088"]], "~/Documents/Analyses/TLS/LU088_TLS_counts.csv")
write.csv(count_list[["LU091"]], "~/Documents/Analyses/TLS/LU091_TLS_counts.csv")



#also export B cell percentages to Prism
b_cell <- write.csv(fragments_small$B_cell_percentage, "~/Documents/Analyses/TLS/b_cell_percentages.csv")


#count the amount of fragments, threshold-passed and TLS in each tumour
tumour_fragment_counts <- fragments_noB %>%
  group_by(tumour) %>%
  dplyr::count()

cell_threshold_counts <- fragments_noB %>%
  group_by(tumour,cell_threshold_100) %>%
  filter(cell_threshold_100 == "YES") %>%
  dplyr::count()

TLS_counts <- fragments_noB %>%
  group_by(tumour,cell_threshold_100,segment) %>%
  filter(cell_threshold_100 == "YES") %>%
  filter(segment == "TLS") %>%
  dplyr::count()

LU080 <- data.frame(tumour = "LU080", cell_threshold_100 = "YES", segment = "TLS", TLS = 0) #LU080 has no TLS and needs to be added manually
TLS_counts <- rbind(TLS_counts,LU080)
TLS_counts <- TLS_counts[order(TLS_counts$tumour),]

colnames(tumour_fragment_counts) <- c("tumour", "fragments_found")
colnames(cell_threshold_counts) <- c("tumour", "cell_threshold_100", "fragments_passed_threshold")
colnames(TLS_counts) <- c("tumour", "cell_threshold_100", "segment", "TLS")

tumour_fragment_counts <- cbind(tumour_fragment_counts, cell_threshold_counts$fragments_passed_threshold)
tumour_fragment_counts <- cbind(tumour_fragment_counts, TLS_counts$TLS)

colnames(tumour_fragment_counts) <- c("tumour", "fragments_found","fragments_passed_threshold","TLS")

#reshaping for proper plotting
long <- melt(setDT(tumour_fragment_counts), id.vars = c("tumour"), variable.name = "count")

ggplot(long,aes(x=tumour,y=value,fill=factor(count,levels = c("fragments_found","fragments_passed_threshold","TLS")))) + geom_col(position = "dodge") +
  ylab("count") + xlab("tumour") + scale_fill_manual(values = c("grey","#FF69B4", "#99FFCC")) +theme_minimal() +guides(fill=guide_legend(title="fragment")) + theme(text=element_text(size=20))




#is there a baseline difference between HLA-DR expression in tumour and TLS?
TLS_MFI <- fragments_small[fragments_small$segment == "TLS", "HLA-DR_B_cell_MFI"]
tumour_MFI <- fragments_small[fragments_small$segment == "tumour_bed", "HLA-DR_B_cell_MFI"]
tumour_MFI <- tumour_MFI[!is.na(tumour_MFI)]

t.test(TLS_MFI, tumour_MFI, alternative = "two.sided") #no difference at all

mean(TLS_MFI)
mean(tumour_MFI)



#remove intermediate segments
fragments_smaller <- fragments_small[fragments_small$segment != "intermediate",]
write.csv(fragments_smaller, "~/Documents/Analyses/TLS/fragment_smaller.csv")


#plot HLA-DR MFI
ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `HLA-DR_B_cell_MFI`, fill = segment)) +
  geom_boxplot() + geom_jitter(alpha = 0.2)  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("HLA-DR MFI (B cells)") #+ geom_signif(comparisons = list(c("unstim", "aPD1"),c("aPD1", "aIFNyR1+aPD1"),c("aIFNyR1+aPD1", "hrIFNy")), map_signif_level=TRUE,tip_length = 0)

library(ggprism)
ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `HLA-DR_B_cell_MFI`, fill = segment)) +
  geom_boxplot() + geom_jitter(alpha = 0.2)  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_prism() + theme(text=element_text(size=20)) + xlab("condition") +ylab("HLA-DR MFI (B cells)") #+ geom_signif(comparisons = list(c("unstim", "aPD1"),c("aPD1", "aIFNyR1+aPD1"),c("aIFNyR1+aPD1", "hrIFNy")), map_signif_level=TRUE,tip_length = 0)

#export to prism
test <- fragments_smaller[,c("condition","segment","HLA-DR_B_cell_MFI")]
test <- test[order(test$condition, test$segment),]
test$indx <-with(test, ave(seq_along(condition), segment, condition, FUN=seq_along))
test2 <- dcast(test, condition~segment+indx, value.var="HLA-DR_B_cell_MFI")
test2<- rbind(test2, colnames(test2))
test2[6,] <- ifelse(grepl("TLS",test2[6,]), "TLS","Tumour bed")
test2[6,1] <- "condition"
export <- rbind(test2[6,], test2[1:5,])
write.xlsx(export, "~/Documents/Analyses/TLS/HLA_DR_MFI.xlsx")

test <- fragments_smaller[,c("condition","segment","HLA-DR_B_cell_MFI")]
test <- test[order(test$condition, test$segment),]
write.csv(test, "~/Documents/Analyses/TLS/HLA_DR_MFI.csv")






#compare all conditions
one.way <- aov(`HLA-DR_B_cell_MFI`~ condition, data = fragments_smaller)
summary(one.way)

#include the segment status (TLS/non-TLS)
two.way <- aov(`HLA-DR_B_cell_MFI` ~ condition+segment, data = fragments_smaller)
summary(two.way)

#include segment status as an interaction with the condition
interaction <- aov(`HLA-DR_B_cell_MFI` ~ condition*segment, data = fragments_smaller)
summary(interaction)

#use tumour as a blocking variable
block <- aov(`HLA-DR_B_cell_MFI` ~ condition*segment + tumour, data = fragments_smaller)
summary(block)


#individual t-tests
TLS_unstim_MFI <- fragments_small[fragments_small$segment == "TLS" & fragments_small$condition == "unstim", "HLA-DR_B_cell_MFI"]
TLS_aPD1_MFI <- fragments_small[fragments_small$segment == "TLS" & fragments_small$condition == "aPD1", "HLA-DR_B_cell_MFI"]

t.test(TLS_unstim_MFI, TLS_aPD1_MFI, alternative = "two.sided")

#only including tumours with both unstim and aPD1 observations
fragments_unstim_aPD1 <- fragments_small[fragments_small$tumour %in% c("LU035","LU052"),]
TLS_unstim_MFI <- fragments_unstim_aPD1[fragments_unstim_aPD1$segment == "TLS" & fragments_unstim_aPD1$condition == "unstim", "HLA-DR_B_cell_MFI"]
TLS_aPD1_MFI <- fragments_unstim_aPD1[fragments_unstim_aPD1$segment == "TLS" & fragments_unstim_aPD1$condition == "aPD1", "HLA-DR_B_cell_MFI"]

t.test(TLS_unstim_MFI, TLS_aPD1_MFI, alternative = "two.sided")

#comparing aPD1 to aPD1+aIFNR
TLS_combi_MFI <- fragments_small[fragments_small$segment == "TLS" & fragments_small$condition == "aIFNyR1+aPD1", "HLA-DR_B_cell_MFI"]
TLS_aPD1_MFI <- fragments_small[fragments_small$segment == "TLS" & fragments_small$condition == "aPD1", "HLA-DR_B_cell_MFI"]

t.test(TLS_combi_MFI, TLS_aPD1_MFI, alternative = "two.sided")


#calculate difference in CD25+ CD4+ T cells
TLS_reg <- fragments_small[fragments_small$segment == "TLS", "CD4+ CD25+_percentage"]
norm_reg <- fragments_small[fragments_small$segment == "tumour_bed", "CD4+ CD25+_percentage"]

t.test(TLS_reg, norm_reg, alternative = "two.sided")

ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD4+ CD25+_percentage`, fill = segment)) +
  geom_boxplot()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD4+ CD25+ T cells (% of CD4)") #add scatter


ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD4+ CD25+_percentage`, fill = segment)) +
  geom_violin()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD4+ CD25+ T cells (% of CD4)") #add scatter


ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD8+ CD39+ CD137+_percentage`, fill = segment)) +
  geom_boxplot()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition")# +ylab("CD4+ CD25+ T cells (% of CD4)") #add scatter


#COMPARISONS
#CD80 MFI
ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD80_B_cell_MFI`, fill = segment)) +
  geom_boxplot()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD80 B cell MFI") #add scatter

#CD86 MFI
ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD86_B_cell_MFI`, fill = segment)) +
  geom_boxplot()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD86 B cell MFI") #add scatter

#CD4+ PD-1+
ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD4+ PD-1+_percentage`, fill = segment)) +
  geom_boxplot()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD4+ PD-1+ (of CD4+)") #add scatter

#CD8+ PD-1+
ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD8+ PD-1+_percentage`, fill = segment)) +
  geom_boxplot()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD8+ PD-1+ (of CD8+)") #add scatter

#CD8+ CD25+
ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `CD8+ CD25+_percentage`, fill = segment)) +
  geom_boxplot()  + scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD8+ CD25+ (of CD8+)") #add scatter

#Collect data for prism
fragments_small$segment <- factor(fragments_small$segment, levels = c("TLS", "intermediate","tumour_bed"))
fragments_summarised <- split(fragments_small, fragments_small$segment)

fragments_summarised_final <- list()
condition <- c("unstim", "aPD1", "aIFNyR1+aPD1", "hrIFNy", "hrIFNy_high")

for(i in seq_along(fragments_summarised)){
  fragments_summarized_data <- fragments_summarised[[i]] %>%
    group_by(tumour, condition) %>%
    summarise(count = n(), across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) #take colMeans in all numeric columns in each group
  
  fragments_summarized_data[fragments_summarized_data == "NaN"] <- NA
  
  fragments_summarized_data <- fragments_summarized_data[order(fragments_summarized_data$tumour, match(fragments_summarized_data$condition, condition)),]
  fragments_summarized_data$condition <- factor(fragments_summarized_data$condition, levels = condition)
  
  fragments_summarized_data <- fragments_summarized_data[fragments_summarized_data$condition != "hrIFNy_high",]
  
  fragments_summarised_final[[length(fragments_summarised_final) + 1]] <- fragments_summarized_data
  names(fragments_summarised_final)[[length(fragments_summarised_final)]] <- names(fragments_summarised)[[i]]
}


#individual plot test
#broken on column_y, but only for lists
plot_condition <- function(fragments_summarised_final, variable){
  plot_list <- list()
  for (i in seq_along(fragments_summarised_final)){
    column_y <- fragments_summarised_final[[i]][,colnames(fragments_summarised_final[[i]]) == variable]
   column_y <- c(column_y)[[1]]
    plot <- ggplot(fragments_summarised_final[[i]], aes(x = condition, y = column_y, fill = tumour, group=factor(tumour))) +
      geom_point(position = position_dodge(width=0.4), aes(fill = tumour)) + geom_line(position = position_dodge(width=0.4), aes(col=tumour)) + theme_bw() +
      ggtitle(label = names(fragments_summarised_final)[[i]])
    
    plot_list[[length(plot_list) + 1]] <- plot
  }
  cowplot::plot_grid(plotlist = plot_list, nrow = 1)
  return(plot_list)
}

hla_dr <- plot_condition(fragments_summarised_final, "HLA-DR_B_cell_MFI")


#mass export to prism
for(i in seq_along(fragments_summarised_final)){
  data <- as.data.frame(fragments_summarised_final[[i]])
  for (j in 1:ncol(fragments_summarised_final[[i]])){
    if(is.numeric(data[,j])){
      colname <- colnames(data)[[j]]
      long <- data[,c("tumour","condition", colname)]
      colnames(long) <- c("tumour", "condition", "value")
      wide <- reshape(long, idvar = "tumour", timevar = "condition", direction = "wide")
      colnames(wide)[colnames(wide) == "value.unstim"] <- "unstim"
      colnames(wide)[colnames(wide) == "value.aPD1"] <- "aPD1"
      colnames(wide)[colnames(wide) == "value.aIFNyR1+aPD1"] <- "aIFNyR1+aPD1"
      colnames(wide)[colnames(wide) == "value.hrIFNy"] <- "hrIFNy"
      
      file_location <- paste0("~/Documents/Analyses/TLS/mass_data_export/",names(fragments_summarised_final)[[i]],"_",colname,".csv")
      
      write.csv(wide,file_location)
    }
  }
}

#calculate cellular compositions
composition <- fragments_small[,c("Name","cells","experiment","tumour","condition","fragment", "segment","B_cell_percentage","CD3+_percentage","CD3+ DN_percentage",
                                  "CD3+ DP_percentage","CD4+_percentage","CD8+_percentage")]
composition$CD3_CD19_DN <- 100-(composition$B_cell_percentage + composition$`CD3+_percentage`)

library(data.table)
composition_long <- melt(as.data.table(composition),id.vars = c("Name","cells","experiment","tumour","condition","fragment", "segment"), value.name = "percentage", variable.name = "population")
composition_long$population <- gsub("_percentage", "", composition_long$population)

ggplot(composition_long, aes(x=factor(fragment), y = percentage, fill = population)) + geom_col(position = "stack") +
  scale_fill_manual(values = c("pink","orange","red","lightblue","lightgreen","darkgreen","purple")) +theme_bw() + facet_wrap(~tumour)

##----------------
#correct the MFI: calculate the MFI as a ratio of the MFI in PBMCs with the same fragment number
#prepare PBMC data
PBMC_TLS009 <- flow_file_TLS009[grepl("3. pool unstim|Median",flow_file_TLS009$Name),]
PBMC_TLS009 <- PBMC_TLS009[1:403,] #remove median values for non-PBMC samples, change number if amount of gates is different

PBMC_TLS010 <- flow_file_TLS010[grepl("3. pool unstim|Median",flow_file_TLS010$Name),]
PBMC_TLS010 <- PBMC_TLS010[1:396,] #remove median values for non-PBMC samples, change number if amount of gates is different

PBMC_TLS011 <- flow_file_TLS011[grepl("PBMC unstim|Median",flow_file_TLS011$Name),]
PBMC_TLS011 <- PBMC_TLS011[1:357,] #remove median values for non-PBMC samples, change number if amount of gates is different

PBMC_TLS012 <- flow_file_TLS012[grepl("PBMC unstim|Median",flow_file_TLS012$Name),]
PBMC_TLS012 <- PBMC_TLS012[170:492,] #remove median values for non-PBMC samples, change number if amount of gates is different

PBMC_TLS013 <- flow_file_TLS013[grepl("PBMC pool unstim|Median",flow_file_TLS013$Name),]
PBMC_TLS013 <- PBMC_TLS013[283:745,] #remove median values for non-PBMC samples, change number if amount of gates is different

#make a list for easier processing
PBMC_list <- list(PBMC_TLS009,PBMC_TLS010,PBMC_TLS011,PBMC_TLS012,PBMC_TLS013)

names(PBMC_list) <- c("TLS009","TLS010","TLS011","TLS012","TLS013")


#couple the B cell data as done for the tumour fragments
for(j in seq_along(PBMC_list)){
  #grab the entries that are needed
  colnames(PBMC_list[[j]]) <- c("Name","percent_of_live","cells")
  PBMC_list[[j]] <- PBMC_list[[j]][!grepl("CD3|T cell|Myeloid|CD20",PBMC_list[[j]]$Name,ignore.case = T),]
  PBMC_list[[j]] <- PBMC_list[[j]][grepl("#|Median",PBMC_list[[j]]$Name,ignore.case = T),]
  
  #setup for loop
  MFI_expression_file <- data.frame()
  flow_file <- PBMC_list[[j]]
  
  for(o in 1:nrow(flow_file)){
    if(grepl("Median", flow_file$Name[[o]])){
      
      #grab median expression value
      MFI_row_stat <- strsplit(flow_file$Name[[o]],"=")
      MFI_row_marker <- MFI_row_stat[[1]][1]
      MFI_row_marker <- gsub("Median : Comp-|-A| " , "", MFI_row_marker)
      
      MFI_row_stat <- MFI_row_stat[[1]][2]
      MFI_row_stat <- gsub(" ", "", MFI_row_stat) #remove spaces
      
      #find the row it corresponds to
      for(x in (1:o)-1){
        #look in all rows above the entry
        if(grepl("PBMC|stim|unstim", flow_file$Name[[o-x]])){
          row <- flow_file$Name[[o-x]] #if there is a PBMC sample found, couple the value to this sample and exit the loop
          break
        } else {
          next #if no sample is found on this row, go to the row above it. No more failsaves are needed as it should always find the correct sample above it
        }
      }  #end of x loop
      
      #quick check: if no B cells are there, the MFI will be n/a, so leave this entry out
      if(is.na(MFI_row_stat)){
        next
      }
      
      
      #make the data entry
      row_full <- data.frame(Name = row, marker = MFI_row_marker, MFI = MFI_row_stat)
      
      #append to data.frame
      MFI_expression_file <- rbind(MFI_expression_file, row_full)
    } 
  } #end of o loop
  
  #append to data.frame
  MFI_expression_file <- MFI_expression_file[!is.na(MFI_expression_file$Name),]
  MFI_expression_file$MFI <- as.numeric(MFI_expression_file$MFI)
  MFI_expression_file$marker[MFI_expression_file$marker == "BB515"] <- "HLA-DR"
  MFI_expression_file$marker[MFI_expression_file$marker == "BUV563"] <- "CD86"
  MFI_expression_file$marker[MFI_expression_file$marker == "PE"] <- "CD80"
  
  #reshaping
  MFI_expression_file_wide <- spread(MFI_expression_file, key = marker, value = MFI)
  MFI_expression_file_wide$Name <- gsub("/B cells", "", MFI_expression_file_wide$Name)
  
  PBMC_list[[j]] <- left_join(PBMC_list[[j]],MFI_expression_file_wide,by="Name")
  
  #add B cell information
  PBMC <- PBMC_list[[j]][!grepl("B cell|Median", PBMC_list[[j]]$Name),]
  PBMC$B_cell_percentage <- NA
  PBMC$B_cells <- NA
  
  for (s in 1:nrow(PBMC)){
    name <- PBMC$Name[[s]]
    name_B_cell <- paste0(name, "/B cells") #grab flowjo path of B cell gate
    B_cell <- PBMC_list[[j]][PBMC_list[[j]]$Name %in% name_B_cell,]
    B_cell <- B_cell[1,]
    
    #add information to sample
    PBMC$B_cell_percentage[[s]] <- B_cell$percent_of_live
    PBMC$B_cells[[s]] <- B_cell$cells
    
  }
  
  PBMC_list[[j]] <- PBMC
  
  #get the fragment number out
  names <- strsplit(PBMC_list[[j]]$Name,"#")
  names <- sapply(names,"[[",2)
  PBMC_list[[j]]$fragment <- names
  
  #column renaming
  colnames(PBMC_list[[j]])[colnames(PBMC_list[[j]]) == "HLA-DR"] <- "HLA.DR_B_cell_MFI"
  colnames(PBMC_list[[j]])[colnames(PBMC_list[[j]]) == "CD80"] <- "CD80_B_cell_MFI"
  colnames(PBMC_list[[j]])[colnames(PBMC_list[[j]]) == "CD86"] <- "CD86_B_cell_MFI"
  
}


#assign experiment numbers to each tumour
fragments_small$experiment <- NA
for (i in 1:nrow(fragments_small)){
  if(fragments_small$tumour[[i]] %in% c("LU086")){
    fragments_small$experiment[[i]] <- "TLS009"
  }
  if(fragments_small$tumour[[i]] %in% c("LU052")){
    fragments_small$experiment[[i]] <- "TLS010"
  }
  if(fragments_small$tumour[[i]] %in% c("LU086_1","LU086_2","LU080")){
    fragments_small$experiment[[i]] <- "TLS011"
  }
  if(fragments_small$tumour[[i]] %in% c("LU066","LU083")){
    fragments_small$experiment[[i]] <- "TLS012"
  }
  if(fragments_small$tumour[[i]] %in% c("LU023","LU035")){
    fragments_small$experiment[[i]] <- "TLS013"
  }
}


#compare the MFI to PBMCs of the same experiment
fragments_small$CD80_B_cell_MFI_ratio <- NA
fragments_small$CD86_B_cell_MFI_ratio <- NA
fragments_small$HLA_DR_B_cell_MFI_ratio <- NA


for (i in 1:nrow(fragments_small)){
    PBMC_ref <- PBMC_list[[fragments_small$experiment[[i]]]] #select the correct PBMC reference

    fragment <- strsplit(fragments_small$fragment[[i]],"_") #to ensure fragments like 1_2 are handled correctly
    fragment <- fragment[[1]][1]
    
    if(fragments_small$experiment[[i]] == "TLS011" & fragment == "9"){ #this fragment was not found in PBMCs
      next
    }
    
    #grab the correct PBMC fragment and 
    PBMC_values <- PBMC_ref[PBMC_ref$fragment == fragment, colnames(PBMC_ref) %in% c("CD80_B_cell_MFI","CD86_B_cell_MFI","HLA.DR_B_cell_MFI")]
    fragment_values <- fragments_small[i,colnames(fragments_small[i,]) %in% c("CD80_B_cell_MFI","CD86_B_cell_MFI","HLA.DR_B_cell_MFI")]
    
    if("HLA.DR_B_cell_MFI" %in% colnames(PBMC_values)){
      HLA_DR_ratio <- fragment_values$HLA.DR_B_cell_MFI/PBMC_values$HLA.DR_B_cell_MFI
      fragments_small$HLA_DR_B_cell_MFI_ratio[[i]] <- HLA_DR_ratio
    }
    if("CD80_B_cell_MFI" %in% colnames(PBMC_values)){
      CD80_ratio <- fragment_values$CD80_B_cell_MFI/PBMC_values$CD80_B_cell_MFI
      fragments_small$CD80_B_cell_MFI_ratio[[i]] <- CD80_ratio
    }
    if("CD86_B_cell_MFI" %in% colnames(PBMC_values)){
      CD86_ratio <- fragment_values$CD86_B_cell_MFI/PBMC_values$CD86_B_cell_MFI
      fragments_small$CD86_B_cell_MFI_ratio[[i]] <- CD86_ratio
    }
}

#missing ratios are either due to a lack of B cells or no PBMC reference existing

#plot ratios
fragments_smaller <- fragments_small[fragments_small$segment != "intermediate",]
write.csv(fragments_smaller, "~/Documents/Analyses/TLS/fragment_smaller.csv")

ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= `HLA_DR_B_cell_MFI_ratio`, fill = segment)) +
  geom_boxplot()+ scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("HLA-DR MFI on B cells (ratio fragment/PBMC)") 

ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= CD80_B_cell_MFI_ratio, fill = segment)) +
  geom_boxplot()+ scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD80 MFI on B cells (ratio fragment/PBMC)") 

ggplot(fragments_smaller, aes(x=factor(condition, levels = c("unstim","aPD1","aIFNyR1+aPD1","hrIFNy","hrIFNy_high")), y= CD86_B_cell_MFI_ratio, fill = segment)) +
  geom_boxplot()+ scale_fill_manual(values = c("#99FFCC","#99CCFF")) +
  theme_minimal() + theme(text=element_text(size=20)) + xlab("condition") +ylab("CD86 MFI on B cells (ratio fragment/PBMC)") 


#Averaging across TLS/non-TLS
TLS <- fragments_small[fragments_small$segment == "TLS",]
tumour <- fragments_small[fragments_small$segment == "tumour_bed",]

TLS_unstim <- TLS[TLS$condition == "unstim",]
TLS_PD1 <- TLS[TLS$condition == "aPD1",]
TLS_PD1_IFN <- TLS[TLS$condition == "aIFNyR1+aPD1",]
TLS_IFN <- TLS[TLS$condition == "hrIFNy",]

#take mean of every statistic per tumour
TLS_list <- list(TLS_unstim, TLS_PD1, TLS_PD1_IFN, TLS_IFN)
names(TLS_list) <- c("TLS unstim", "TLS aPD1", "TLS aPD1_aIFNR", "TLS hrIFNy")
TLS_concat <- list()

for (i in seq_along(TLS_list)){
  data <- aggregate(TLS_list[[i]][, c(2:3,5:9,15:90)], list(TLS_list[[i]]$tumour), mean)
  TLS_concat[[i]] <- data
}
names(TLS_concat) <- names(TLS_list)

#for tumour bed fragments:
tumour_unstim <- tumour[tumour$condition == "unstim",]
tumour_PD1 <- tumour[tumour$condition == "aPD1",]
tumour_PD1_IFN <- tumour[tumour$condition == "aIFNyR1+aPD1",]
tumour_IFN <- tumour[tumour$condition == "hrIFNy",]

#take mean of every statistic per tumour
tumour_list <- list(tumour_unstim, tumour_PD1, tumour_PD1_IFN, tumour_IFN)
names(tumour_list) <- c("tumour unstim", "tumour aPD1", "tumour aPD1_aIFNR", "tumour hrIFNy")
tumour_concat <- list()

for (i in seq_along(tumour_list)){
  data <- aggregate(tumour_list[[i]][, c(2:3,5:9,15:90)], list(tumour_list[[i]]$tumour), mean)
  tumour_concat[[i]] <- data
}
names(tumour_concat) <- names(tumour_list)

#get the data out
TLS_data <- do.call(rbind.data.frame, TLS_concat)
TLS_data$condition <- NA
TLS_data$condition[(grepl("unstim", rownames(TLS_data)))] <- "unstim"
TLS_data$condition[(grepl("aPD1", rownames(TLS_data))) &!grepl("IFN", rownames(TLS_data))] <- "aPD1"
TLS_data$condition[(grepl("aPD1_aIFNR", rownames(TLS_data)))] <- "aIFNyR1+aPD1"
TLS_data$condition[(grepl("hrIFNy", rownames(TLS_data)))] <- "hrIFNy"
colnames(TLS_data)[[1]] <- "tumour"
TLS_data$segment <- "TLS"

tumour_data <- do.call(rbind.data.frame, tumour_concat)
tumour_data$condition <- NA
tumour_data$condition[(grepl("unstim", rownames(tumour_data)))] <- "unstim"
tumour_data$condition[(grepl("aPD1", rownames(tumour_data))) &!grepl("IFN", rownames(tumour_data))] <- "aPD1"
tumour_data$condition[(grepl("aPD1_aIFNR", rownames(tumour_data)))] <- "aIFNyR1+aPD1"
tumour_data$condition[(grepl("hrIFNy", rownames(tumour_data)))] <- "hrIFNy"
colnames(tumour_data)[[1]] <- "tumour"
tumour_data$segment <- "tumour_bed"

data_final <- rbind(TLS_data, tumour_data)
data_final$category <- paste0(data_final$segment, "_", data_final$condition, "_", data_final$tumour)

CD8_test <- data_final[,c("tumour","condition", "segment", "category","CD8+ CD39+ CD137+_percentage")]
CD8_test_wide <- reshape(CD8_test, idvar = "tumour", timevar = "condition", direction = "wide")
write.csv(CD8_test, "~/Documents/Analyses/TLS/CD8_CD39_CD137_test.csv",row.names = F)

#------------------
#cell proportions test: can we plot different CD4+ and CD8+ populations?
cd8 <- fragments_small[,c("Name","experiment","tumour", "condition", "fragment","segment","CD8+_percentage","CD8+ CD25+_percentage","CD8+ CD39+_percentage","CD8+ CD103+_percentage","CD8+ CD137+_percentage","CD8+ OX40+_percentage","CD8+ PD-1+_percentage","CD8+ IL7R+_percentage")]
summary_cd8 <- cd8 %>%
  group_by(condition, segment) %>%
  summarize(across(6:12, \(x) mean(x, na.rm = TRUE)))

colnames(summary_cd8) <- gsub("_percentage", "", colnames(summary_cd8))

summary_cd8_long <- melt(summary_cd8)

summary_cd8_long$condition <- factor(summary_cd8_long$condition, levels = c("unstim", "aPD1", "aIFNyR1+aPD1", "hrIFNy", "hrIFNy_high"))
summary_cd8_long <- summary_cd8_long[summary_cd8_long$condition != "hrIFNy_high",]

summary_cd8_long$segment <- factor(summary_cd8_long$segment, levels = c("TLS", "intermediate", "tumour_bed"))

ggplot(summary_cd8_long, aes(x=segment, y = value, fill = variable)) + geom_col() + facet_grid(~condition) + theme_minimal() +
  scale_fill_manual(values = c("darkgoldenrod1", "darkred", "pink","darkcyan","chartreuse4","cornflowerblue","purple"))



cd4 <- fragments_small[,c("Name","experiment","tumour", "condition", "fragment","segment","CD4+_percentage","CD4+ CD25+_percentage","CD4+ CD39+_percentage","CD4+ CD103+_percentage","CD4+ CD137+_percentage","CD4+ OX40+_percentage","CD4+ PD-1+_percentage","CD4+ IL7R+_percentage")]
summary_cd4 <- cd4 %>%
  group_by(condition, segment) %>%
  summarize(across(6:12, \(x) mean(x, na.rm = TRUE)))

colnames(summary_cd4) <- gsub("_percentage", "", colnames(summary_cd4))

summary_cd4_long <- melt(summary_cd4)

summary_cd4_long$condition <- factor(summary_cd4_long$condition, levels = c("unstim", "aPD1", "aIFNyR1+aPD1", "hrIFNy", "hrIFNy_high"))
summary_cd4_long <- summary_cd4_long[summary_cd4_long$condition != "hrIFNy_high",]

summary_cd4_long$segment <- factor(summary_cd4_long$segment, levels = c("TLS", "intermediate", "tumour_bed"))

ggplot(summary_cd4_long, aes(x=segment, y = value, fill = variable)) + geom_col() + facet_grid(~condition) + theme_minimal() +
  scale_fill_manual(values = c("darkgoldenrod1", "darkred", "pink","darkcyan","chartreuse4","cornflowerblue","purple"))


#correlations
fragment_noPD1 <- fragments_small[fragments_small$condition %in% c("unstim", "hrIFNy"),]

#first test aPD1, which can only be done in conditions where it is not blocked
cor_coefs <- cor.test(fragment_noPD1$B_cell_percentage, fragment_noPD1$`CD8+ PD-1+_percentage`)
ggplot(fragment_noPD1, aes(x= B_cell_percentage, y = `CD8+ PD-1+_percentage`)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  annotate("text", x = 0, y = 112, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 15, y = 112, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))

cor_coefs <- cor.test(fragment_noPD1$B_cell_percentage, fragment_noPD1$`CD4+ PD-1+_percentage`)
ggplot(fragment_noPD1, aes(x= B_cell_percentage, y = `CD4+ PD-1+_percentage`)) + geom_point() + geom_smooth(method='lm', formula= y~x) +
  annotate("text", x = 0, y = 112, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 15, y = 112, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))


plot_correlation <- function(x, y){
  column_x <- fragments_small[,colnames(fragments_small) == x]
  column_y <- fragments_small[,colnames(fragments_small) == y]
  
  cor_coefs <- cor.test(x = column_x, y = column_y)
  ggplot(fragments_small, aes(x=column_x, y = column_y)) + geom_point() + geom_smooth(method='lm', formula= y~x) + xlab(x) + ylab(y) +
    annotate("text", x = 0, y = 112, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
    annotate("text", x = 15, y = 112, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))
}


plot_correlation("B_cell_percentage", "CD8+ CD39+ CD137+_percentage")
plot_correlation("B_cell_percentage", "CD4+ IL7R+_percentage")
plot_correlation("B_cell_percentage", "CD8+ IL7R+_percentage")
plot_correlation("B_cell_percentage", "CD8+ CD25+_percentage")


plot_correlation("B_cell_percentage", "CD4+_percentage")
plot_correlation("B_cell_percentage", "CD8+_percentage")

#plot the correlation between CD4 and CD8 percentage (of CD45) and the B cell percentage
cd4 <- fragments_small[,c("Name", "B_cell_percentage", "CD3+_percentage", "CD4+_percentage", "CD8+_percentage")]
cd4$`CD4+_of_CD45` <- (cd4$`CD3+_percentage`/100) * (cd4$`CD4+_percentage`/100) *100
cd4$`CD8+_of_CD45` <- (cd4$`CD3+_percentage`/100) * (cd4$`CD8+_percentage`/100) *100

cd4$cd4_cd8_ratio <- cd4$`CD4+_percentage`/cd4$`CD8+_percentage`
cd4$cd4_cd8_ratio_cd45 <- cd4$`CD4+_of_CD45`/cd4$`CD8+_of_CD45`


write.csv(cd4, "~/Documents/Analyses/TLS/mass_data_export/CD4_of_CD45.csv")

cor_coefs <- cor.test(x = cd4$B_cell_percentage, y = cd4$`CD4+_of_CD45`)
ggplot(cd4, aes(x=B_cell_percentage, y = `CD4+_of_CD45`)) + geom_point() + geom_smooth(method='lm', formula= y~x) + xlab("B cell percentage") + ylab("CD4+ (%of CD45+)") +
  annotate("text", x = 0, y = 112, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 15, y = 112, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))

cor_coefs <- cor.test(x = cd4$B_cell_percentage, y = cd4$`CD8+_of_CD45`)
ggplot(cd4, aes(x=B_cell_percentage, y = `CD8+_of_CD45`)) + geom_point() + geom_smooth(method='lm', formula= y~x) + xlab("B cell percentage") + ylab("CD8+ (%of CD45+)") +
  annotate("text", x = 0, y = 112, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 15, y = 112, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))

cor_coefs <- cor.test(x = cd4$B_cell_percentage, y = cd4$cd4_cd8_ratio)
ggplot(cd4, aes(x=B_cell_percentage, y = cd4_cd8_ratio)) + geom_point() + geom_smooth(method='lm', formula= y~x) + xlab("B cell percentage") + ylab("Rato of CD4+/CD8+") +
  annotate("text", x = 0, y = 112, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 15, y = 112, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))

cor_coefs <- cor.test(x = cd4$B_cell_percentage, y = cd4$cd4_cd8_ratio_cd45)
ggplot(cd4, aes(x=B_cell_percentage, y = cd4_cd8_ratio_cd45)) + geom_point() + geom_smooth(method='lm', formula= y~x) + xlab("B cell percentage") + ylab("Rato of CD4+/CD8+") +
  annotate("text", x = 0, y = 112, label = paste0("R: ", round(cor_coefs$estimate, 2))) +
  annotate("text", x = 15, y = 112, label = paste0("p-value: ", round(cor_coefs$p.value, 10)))

ggplot(fragments_small, aes(x=segment, y = percent_of_CD45)) + geom_boxplot() +geom_jitter()

#calculate whether TLS segments have more cells than other segments
aov_cd45 <- aov(percent_of_CD45 ~ segment, data = fragments_small)
summary(aov_cd45)
TukeyHSD(aov_cd45)

cd45_small <- fragments_small[,c("Name", "percent_of_CD45", "segment")]
colnames(cd45_small) <- c("Name", "value", "segment")
cd45_small$segment <- factor(cd45_small$segment, levels = c("TLS", "intermediate", "tumour_bed"))
cd45_wide <- reshape(cd45_small, idvar = "Name", timevar = "segment", direction = "wide")
colnames(cd45_wide) <- c("name", "tumour bed", "intermediate", "TLS")

write.csv(cd45_wide, "~/Documents/Analyses/TLS/mass_data_export/CD45_percentage_segment.csv")

ggplot(cd45_small, aes(x=segment, y = value)) + geom_boxplot() +geom_jitter()
