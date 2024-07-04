#Load libraries
library(readxl)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(ggsignif)
library(cowplot)

#----------------
#read in data and add experiment number
flow_file_SLICE004 <- read_excel("~/Documents/Experiments/SLICE/SLICE004_RE181_slice_culture_24h_48h_22042024/SLICE004.xls")
flow_file_SLICE004$experiment <- "SLICE004"
flow_file_SLICE004$Name <- gsub("RE181 unstim 0h", "Fragment unstim 0h", flow_file_SLICE004$Name)

flow_file_SLICE005 <- read_excel("~/Documents/Experiments/SLICE/SLICE005_RE186_72h_03052024/SLICE005.xls")
flow_file_SLICE005$experiment <- "SLICE005"

flow_file_SLICE006 <- read_excel("~/Documents/Experiments/SLICE/SLICE006_RE187_24h_48h_culture_06052024/SLICE006.xls")
flow_file_SLICE006$experiment <- "SLICE006"
flow_file_SLICE006$Depth <- NULL

flow_file_SLICE007 <- read_excel("~/Documents/Experiments/SLICE/SLICE007_cryo_RE187_48h_SZ_NS_22052024/SLICE007.xls")
flow_file_SLICE007$experiment <- "SLICE007"
flow_file_SLICE007$Name <- gsub("cryo ", "", flow_file_SLICE007$Name)

flow_file <- rbind(flow_file_SLICE004, flow_file_SLICE005, flow_file_SLICE006,flow_file_SLICE007)
colnames(flow_file) <- c("Name", "Statistic", "Cells","Experiment")


#grab entries for slices and fragments
fragments <- flow_file[grepl("Fragment|Slice",flow_file$Name, ignore.case = T),]
fragments$Name <- paste0(fragments$Experiment, "_", fragments$Name)

#grab the lines referencing the original file without any gating
fragments_default <- fragments[is.na(fragments$Statistic),]


#remove the LD entry
fragments_default <- fragments_default[!grepl("LD",fragments_default$Name),]


#add gate information
gates <- fragments[!(fragments$Name %in% fragments_default$Name),]
gates <- gates[grepl("CD45|live", gates$Name, ignore.case = T),]
test <- strsplit(gates$Name, "/")

gate_info <- data.frame()

for (i in seq_along(test)){
  gate_name <- NA
  sample <- test[[i]][grepl("\\.fcs",test[[i]])] #sample in which the population was found
  
  if(tail(test[[i]],1) == "CD3\\+ T cell"){ #add numbers of CD3+ T cells
    gate_name <- " CD3+ T cell"
  }
  if(tail(test[[i]],1) == "T cells_Myeloids"){ #add numbers of CD3+ T cells
    gate_name <- "T cell+Myeloid"
  }
  
  if(any(grepl("CD4\\+|CD8\\+", test[[i]]))){ #if this is a T cell entry, feed into the T cell pipeline
    result <- test[[i]][grepl("\\+",test[[i]])] #specific population on this line of the file
    
    #give the T cell population a name:
    if(any(grepl("CD4\\+",result))){
      gate_name <- paste0("CD4+ ", tail(result,1))
    } else if(any(grepl("CD8+",result))){
      gate_name <- paste0("CD8+ ", tail(result,1))
    } else if(any(grepl("DN", test[[i]], ignore.case = F))){ #check DN and DP populations in original list as they do not contain a + sign
      gate_name <- "CD3+ DN"
    } else if(any(grepl("DP", test[[i]], ignore.case = F))){
      gate_name <- "CD3+ DP"
    } else {
      gate_name <- "CD3+"
    }
    
  }
  
  if(any(grepl("B cell", test[[i]]))){ #if it is a B cell gate, name the gate as such
    gate_name <- "B cell"
  }
  
  if(any(grepl("CD45", test[[i]])) & !any(grepl("CD3\\+|B cells|Myeloid", test[[i]]))){ #CD45 gates (live/dead)
    if(tail(test[[i]],1) == "CD45- total"){
      gate_name <- "CD45_negative_total"
    } else if(tail(test[[i]],1) == "CD45- live"){
      gate_name <- "CD45_negative_live"
    } else if(tail(test[[i]],1) == "CD45+ total"){
      gate_name <- "CD45_positive_total"
    } else if(tail(test[[i]],1) == "CD45+ live"){
      gate_name <- "CD45_positive_live"
    } else if(tail(test[[i]],1) == "CD45+"){
      gate_name <- "CD45+_of_live"
    }
  }
  
  if(tail(test[[i]], 1) == "Live"){ #Gate for Live cells from Single Cell
    gate_name <- "live_of_ss"
  }
  
  if(any(grepl("Myeloid", test[[i]]))){ #myeloid gating check
    if(any(grepl("NK", test[[i]]))){
      if(tail(test[[i]],1) == "CD16+ NK"){
        gate_name <- "CD16_positive_NK"
      } else if(tail(test[[i]],1) == "CD16- NK"){
        gate_name <- "CD16_negative_NK"
      }
    } else if (tail(test[[i]],1) == "CD16+ NK"){
      gate_name <- "CD56_negative_myeloid"
    } else if (tail(test[[i]],1) == "CD11c+ CD11b+"){
      gate_name <- "CD11c+_CD11b+"
    } else if (tail(test[[i]],1) == "CD14+ monocytes"){
      gate_name <- "CD14+_monocyte"
    } else if (tail(test[[i]],1) == "CD14- DC"){
      gate_name <- "CD14+_DC"
    } else if (tail(test[[i]],1) == "Myeloid"){
      gate_name <- "Myeloid"
    }
  }
  
  #quick check: if the gate is still not identified, skip this entry
  if(is.na(gate_name)){
    next
  }
    
    #grab the number of cells and percentages for that population
    gate_count <- gates[i,"Cells"] #since the list and the dataframe are in the same order, we can just work with the list index for the data frame
    gate_count <- as.numeric(gate_count)
    gate_percentage <- gates[i,"Statistic"]
   gate_percentage <- as.numeric(gate_percentage)
    

  #make and export the data
  gate_row <- data.frame(Name = sample, population = gate_name, percentage = gate_percentage, count = gate_count)
  
  gate_info <- rbind(gate_info, gate_row)
}

#filter out some artefacts
gate_info$population[gate_info$population == "CD4+ CD4+ T cells"] <- "CD4+"
gate_info$population[gate_info$population == "CD8+ CD8+ T cells"] <- "CD8+"

#filter out any percentages that exceed 100 (this might be a FlowJo artefact?)
gate_info <- gate_info[gate_info$percentage <= 100,]

#convert to wide, separately for percentage and count
gate_percentages <- gate_info[,c("Name","population","percentage")]
gate_counts <- gate_info[,c("Name","population","count")]

gate_percentages <- spread(gate_percentages, key = population, value = percentage)
colnames(gate_percentages)[-1] <- paste0(colnames(gate_percentages)[-1], "_percentage")
gate_counts <- spread(gate_counts, key = population, value = count)
colnames(gate_counts)[-1] <- paste0(colnames(gate_counts)[-1], "_counts")

#get the condition and tumour out by strsplit
data <- strsplit(fragments_default$Name, " ")
type <- sapply(data,"[[",1)
type <- gsub("SLICE004_|SLICE005_|SLICE006_|SLICE007_", "", type)
condition <- sapply(data,"[[",2)
condition <- gsub("T0.fcs","unstim",condition)

matrix <- c()
for(i in seq_along(data)){
  if(any(grepl("-matrix",data[[i]]))){
    matrix_i <- "NO"
  } else{
    matrix_i <- "YES"
  }
  matrix <- c(matrix,matrix_i)
}

timepoint <- c()
for(i in seq_along(data)){
  if(any(grepl("24h",data[[i]]))){
    timepoint_i <- "24"
  } else if(any(grepl("48h",data[[i]]))){
    timepoint_i <- "48"
  } else if(any(grepl("72h",data[[i]]))){
    timepoint_i <- "72"
  } else{
    timepoint_i <- "0"
  }
  timepoint <- c(timepoint,timepoint_i)
}

#add the information back to the sample list
fragments_default$tumour <- NA
fragments_default$tumour[fragments_default$Experiment == "SLICE004"] <- "RE181"
fragments_default$tumour[fragments_default$Experiment == "SLICE005"] <- "RE186"
fragments_default$tumour[fragments_default$Experiment == "SLICE006"] <- "RE187"
fragments_default$tumour[fragments_default$Experiment == "SLICE007"] <- "RE187"

fragments_default$condition <- condition
fragments_default$timepoint <- timepoint
fragments_default$matrix <- matrix
fragments_default$type <- type

fragments_default$type[grepl("medium", fragments_default$Name)] <- "Slice medium"

fragments_default$cryopreservation <- "NO"
fragments_default$cryopreservation[fragments_default$Experiment == "SLICE007"] <- "YES"

#add gate information
gate_data_full <- full_join(gate_percentages,gate_counts, by = "Name")
fragments_default <- left_join(fragments_default, gate_data_full, by = "Name")

#add CD45% information
#The median of a marker is always below the population it refers to. Use this:
expression_file <- data.frame()

for(i in 1:nrow(flow_file)){
  if(grepl("Freq", flow_file$Name[[i]])){
    
    #grab median expression value
    row_stat <- as.numeric(flow_file$Statistic[[i]])
    
    #find the row it corresponds to
    for(j in 1:i){
      #look in all rows above the entry
      if(grepl("Fragment|Slice", flow_file$Name[[i-j]], ignore.case = F)){
        row <- flow_file$Name[[i-j]] #if there is a tumour sample found, couple the value to this sample and exit the loop
        break
      } else if (grepl("PBMC", flow_file$Name[[i-j]])) {
        row <- NA #if there is a PBMC sample found, do not enter a sample and exit the loop
        break
      } else {
        next #if no sample is found on this row, go to the row above it
      }
    } 
    
    if(is.na(row)){
      next #if the row is NA, this was a PBMC sample, so do not couple it
    }
    
    #make the data entry
    row_name <- strsplit(row, "/")
    sample_name <- sapply(row_name, "[[", 1)
    sample_name <- paste0(flow_file$Experiment[[i]], "_", sample_name)
    row_name <- sapply(row_name, "[[", length(row_name[[1]]))
    row_full <- data.frame(Name = sample_name, population = row_name, percent_of_CD45 = row_stat)
    
    #append to data.frame
    expression_file <- rbind(expression_file, row_full)
  }
}


CD45_percentages <- spread(expression_file, key = population, value = percent_of_CD45)
colnames(CD45_percentages)[-1] <- paste0(colnames(CD45_percentages)[-1], "_percentage_of_CD45")

fragments_default$Statistic <- NULL

fragments_default <- left_join(fragments_default, CD45_percentages, by = "Name")

#for the aCD3 fragment, calculate the percentage of CD3+ population by subtracting Myeloids from the T_cell_Myeloid population
fragments_default$`CD3+ T cells_percentage_of_CD45`[is.na(fragments_default$`CD3+ T cells_percentage_of_CD45`)] <- fragments_default$`T cells_Myeloids_percentage_of_CD45`[is.na(fragments_default$`CD3+ T cells_percentage_of_CD45`)] - fragments_default$Myeloid_percentage_of_CD45[is.na(fragments_default$`CD3+ T cells_percentage_of_CD45`)]

#other factor settings
fragments_default$timepoint <- factor(fragments_default$timepoint, levels = c("0", "24", "48","72"))
fragments_default$matrix <- as.factor(fragments_default$matrix)
fragments_default$cryopreservation <- as.factor(fragments_default$cryopreservation)
fragments_default$type <- as.factor(fragments_default$type)
fragments_default$condition <- factor(fragments_default$condition, level = c("unstim", "aPD1", "aCD3"))

#save the table
write.csv(fragments_default, "~/Documents/Analyses/SLICE/data_table.csv")


#---------------------------
#read in table
fragments_default <- read.csv("~/Documents/Analyses/SLICE/data_table.csv", check.names = F)
fragments_default$X <- NULL
fragments_default[,1] <- NULL


#plot the number of live cells in each sample
fragments_default.x <- fragments_default[fragments_default$matrix == "YES" & fragments_default$type != "Slice medium",]
fragments_default.x$category <- paste0(fragments_default.x$Experiment, "_", fragments_default.x$type, "_", fragments_default.x$condition)
cols <- pals::alphabet2(n = length(unique(fragments_default.x$category)))
names(cols) <- unique(fragments_default.x$category)

ggplot(fragments_default.x,aes(x=timepoint,y=live_of_ss_percentage,col=category)) +
  ylab("percentage live cells of single cells") + xlab("timepoint") + geom_point(aes(col=category), shape=16, size=3)  +
  geom_line(data = fragments_default.x, aes(group=category), na.rm=T) + theme_bw() +scale_colour_manual(values = cols)



#export to prism
fragments_default_frag <- fragments_default.x[fragments_default.x$type == "Fragment",]
fragments_default_frag$category <- paste0(fragments_default_frag$Experiment, "_", fragments_default_frag$condition)
fragments_default_frag <- fragments_default_frag[,c("timepoint","category","CD45_positive_live_percentage")]
colnames(fragments_default_frag) <- c("timepoint","category","value")
frag_wide <-reshape(fragments_default_frag, idvar = "timepoint", timevar = "category", direction = "wide")
colnames(frag_wide) <- c("timepoint",unique(fragments_default_frag$category))
timepoint <- c(0, 24, 48, 72)
frag_wide <- frag_wide[order(match(frag_wide$timepoint, timepoint)),]
frag_wide[1,grepl("SLICE004", colnames(frag_wide))] <- frag_wide[1,grepl("SLICE004_unstim", colnames(frag_wide))]
frag_wide[1,grepl("SLICE005", colnames(frag_wide))] <- frag_wide[1,grepl("SLICE005_unstim", colnames(frag_wide))]
frag_wide[1,grepl("SLICE006", colnames(frag_wide))] <- frag_wide[1,grepl("SLICE006_unstim", colnames(frag_wide))]
frag_wide[1,grepl("SLICE007", colnames(frag_wide))] <- frag_wide[1,grepl("SLICE007_unstim", colnames(frag_wide))]

write.csv(frag_wide,"~/Documents/Analyses/SLICE/live_CD45_fragment.csv")

fragments_default_slice <- fragments_default.x[fragments_default.x$type == "Slice",]
fragments_default_slice$category<- paste0(fragments_default_slice$Experiment, "_", fragments_default_slice$condition)
fragments_default_slice <- fragments_default_slice[,c("timepoint","category","CD45_positive_live_percentage")]
colnames(fragments_default_slice) <- c("timepoint","category","value")
slice_wide <-reshape(fragments_default_slice, idvar = "timepoint", timevar = "category", direction = "wide")
colnames(slice_wide) <- c("timepoint",unique(fragments_default_slice$category))
slice_wide <- slice_wide[order(match(slice_wide$timepoint, timepoint)),]

slice_wide[1,grepl("SLICE004", colnames(slice_wide))] <- slice_wide[1,grepl("SLICE004_unstim", colnames(slice_wide))]
slice_wide[1,grepl("SLICE005", colnames(slice_wide))] <- slice_wide[1,grepl("SLICE005_unstim", colnames(slice_wide))]
slice_wide[1,grepl("SLICE006", colnames(slice_wide))] <- slice_wide[1,grepl("SLICE006_unstim", colnames(slice_wide))]
slice_wide[1,grepl("SLICE007", colnames(slice_wide))] <- slice_wide[1,grepl("SLICE007_unstim", colnames(slice_wide))]

write.csv(slice_wide,"~/Documents/Analyses/SLICE/live_CD45_slice.csv")



#mass export to prism
fragments_default.x <- fragments_default.x[fragments_default.x$cryopreservation == "NO",]
fragments_default.x$type <- factor(fragments_default.x$type, levels = c("Fragment", "Slice"))
fragments_summarised <- split(fragments_default.x, fragments_default.x$type)

fragments_summarised_final <- list()
condition <- c("unstim", "aPD1", "aCD3")

for(i in seq_along(fragments_summarised)){
  fragments_summarized_data <- fragments_summarised[[i]] %>%
    group_by(condition, timepoint) %>% #group on condition (unstim/aPD1/aCD3) and timepoint (0/24/48/72)
    summarise(count = n(), across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) #take colMeans in all numeric columns in each group
  
  fragments_summarized_data[fragments_summarized_data == "NaN"] <- NA
  
  fragments_summarized_data <- fragments_summarized_data[order(match(fragments_summarized_data$condition, condition), fragments_summarized_data$timepoint),] #order based on condition and timepoint
  fragments_summarized_data$condition <- factor(fragments_summarized_data$condition, levels = condition)
  
  fragments_summarised_final[[length(fragments_summarised_final) + 1]] <- fragments_summarized_data #save data to list
  names(fragments_summarised_final)[[length(fragments_summarised_final)]] <- names(fragments_summarised)[[i]]
}


for(i in seq_along(fragments_summarised_final)){
  data <- as.data.frame(fragments_summarised_final[[i]])
  for (j in 1:ncol(fragments_summarised_final[[i]])){
    if(is.numeric(data[,j])){
      colname <- colnames(data)[[j]] #take the name of the statistic you are trying to convert
      long <- data[,c("condition","timepoint", colname)]
      colnames(long) <- c("condition","timepoint", "value")
      wide <- reshape(long, idvar = "timepoint", timevar = "condition", direction = "wide") #convert from long format to wide
      colnames(wide)[colnames(wide) == "value.unstim"] <- "unstim"
      colnames(wide)[colnames(wide) == "value.aPD1"] <- "aPD1"
      colnames(wide)[colnames(wide) == "value.aCD3"] <- "aCD3"

      file_location <- paste0("~/Documents/Analyses/SLICE/mass_data_export/",names(fragments_summarised_final)[[i]],"_",colname,".csv") #export to csv
      
      write.csv(wide,file_location)
    }
  }
}



#CD45+ cells
ggplot(fragments_default.x,aes(x=timepoint,y=CD45_positive_live_percentage,col=category, shape = type)) +
  ylab("Percentage of live CD45+ cells out of total CD45+") + xlab("timepoint") + geom_point(aes(col=category), size=3)  +
  geom_line(data = fragments_default.x, aes(group=category), na.rm=T) + theme_bw() +scale_colour_manual(values = cols)

#CD45- cells
ggplot(fragments_default.x,aes(x=timepoint,y=CD45_negative_live_percentage,col=category, shape = type)) +
  ylab("Percentage of live CD45- cells out of total CD45-") + xlab("timepoint") + geom_point(aes(col=category), size=3)  +
  geom_line(data = fragments_default.x, aes(group=category), na.rm=T) + theme_bw() +scale_colour_manual(values = cols)

#plot the composition of each sample
small <- fragments_default[,c("Experiment", "condition", "timepoint","matrix","type","B cells_percentage_of_CD45","CD4+ T cells_percentage_of_CD45",
                              "CD8+ T cells_percentage_of_CD45","CD56+ NK cells_percentage_of_CD45","CD56- myeloids_percentage_of_CD45")]

small$composition <- rowSums(small[,6:10])

small <- small[small$type != "Slice medium",]
colnames(small) <- c("Experiment", "condition", "timepoint","matrix","type","B cell","CD4+ T cell",
                     "CD8+ T cell","CD56+ NK cell","CD56- Myeloid", "composition")
small <- small[small$composition != 0,] #remove entries without CD45+ cells
small$Other <- 100-small$composition
small$Other[small$Other < 0] <- 0
small$composition <- NULL

small_long <- melt(setDT(small), id.vars = c("Experiment", "condition", "timepoint","matrix","type"), variable.name = "Population")
small_long <- as.data.frame(small_long)

small_long_004 <- small_long[small_long$Experiment == "SLICE004",]
small_long_005 <- small_long[small_long$Experiment == "SLICE005",]
small_long_005 <- small_long_005[small_long_005$matrix == "YES",]
small_long_006 <- small_long[small_long$Experiment == "SLICE006",]
small_long_007 <- small_long[small_long$Experiment == "SLICE007",]

library(pals)
palette <- pals::alphabet2(n = length(unique(small_long_004$Population)))
names(palette) <- unique(small_long_004$Population)

plot_004 <- ggplot(small_long_004,aes(x=timepoint,y=value,fill=Population)) +
  ylab("percentage of CD45") + xlab("timepoint") + geom_col(position = "stack") + facet_wrap(~type+condition) +
  scale_fill_manual(values= palette)

plot_005 <- ggplot(small_long_005,aes(x=timepoint,y=value,fill=Population)) +
  ylab("percentage of CD45") + xlab("timepoint") + geom_col(position = "stack") + facet_wrap(~type+condition) +
  scale_fill_manual(values= palette)

plot_006 <- ggplot(small_long_006,aes(x=timepoint,y=value,fill=Population)) +
  ylab("percentage of CD45") + xlab("timepoint") + geom_col(position = "stack") + facet_wrap(~type+condition) +
  scale_fill_manual(values= palette)

plot_007 <- ggplot(small_long_007,aes(x=timepoint,y=value,fill=Population)) +
  ylab("percentage of CD45") + xlab("timepoint") + geom_col(position = "stack") + facet_wrap(~type+condition) +
  scale_fill_manual(values= palette)

plot_grid(plot_004,plot_005,plot_006,plot_007,ncol=2,labels = c("SLICE004","SLICE005","SLICE006","SLICE007"))

#export to prism
small_prism <- small[small$matrix == "YES",]
small_prism$category <- paste0(small_prism$Experiment, "_", small_prism$condition, "_", small_prism$timepoint)
small_prism_frag <- small_prism[small_prism$type == "Fragment",]
small_prism_frag[,c("Experiment","condition", "timepoint","matrix","type")] <- NULL
rownames(small_prism_frag) <- small_prism_frag$category

write.csv(small_prism_frag, "~/Documents/Analyses/SLICE/Fragment_composition.csv")

small_prism_slice <- small_prism[small_prism$type == "Slice",]
small_prism_slice[,c("Experiment","condition", "timepoint","matrix","type")] <- NULL
rownames(small_prism_slice) <- small_prism_slice$category
write.csv(small_prism_slice, "~/Documents/Analyses/SLICE/Slice_composition.csv")

#plot aPD1 penetration
small <- fragments_default.x[,c("Experiment", "condition", "timepoint","matrix","type","CD4+ PD-1+_percentage","CD8+ PD-1+_percentage")]

small_long <- melt(setDT(small), id.vars = c("Experiment", "condition", "timepoint","matrix","type"), variable.name = "population")
small_long <- as.data.frame(small_long)

small_long$population <- gsub("_percentage", "", small_long$population)

small_long_004 <- small_long[small_long$Experiment == "SLICE004",]
small_long_005 <- small_long[small_long$Experiment == "SLICE005",]
small_long_006 <- small_long[small_long$Experiment == "SLICE006",]
small_long_007 <- small_long[small_long$Experiment == "SLICE007",]

pd1_004 <- ggplot(small_long_004,aes(x=timepoint,y=value,fill=population)) +
  ylab("percentage of parent") + xlab("timepoint") + geom_col(position = "dodge") + facet_wrap(~type+condition) +
  scale_fill_manual(values= c("darkgreen","darkcyan"))

pd1_005 <- ggplot(small_long_005,aes(x=timepoint,y=value,fill=population)) +
  ylab("percentage of parent") + xlab("timepoint") + geom_col(position = "dodge") + facet_wrap(~type+condition) +
  scale_fill_manual(values= c("darkgreen","darkcyan"))

pd1_006 <- ggplot(small_long_006,aes(x=timepoint,y=value,fill=population)) +
  ylab("percentage of parent") + xlab("timepoint") + geom_col(position = "dodge") + facet_wrap(~type+condition) +
  scale_fill_manual(values= c("darkgreen","darkcyan"))

pd1_007 <- ggplot(small_long_007,aes(x=timepoint,y=value,fill=population)) +
  ylab("percentage of parent") + xlab("timepoint") + geom_col(position = "dodge") + facet_wrap(~type+condition) +
  scale_fill_manual(values= c("darkgreen","darkcyan"))

plot_grid(pd1_004,pd1_005,pd1_006,pd1_007, labels = c("SLICE004","SLICE005","SLICE006","SLICE007"))

#investigate effects of matrigel
matrix <- fragments_default[fragments_default$Experiment == "SLICE005",]
matrix_medium <- matrix[matrix$type == "Slice medium",]
ggplot(matrix_medium,aes(x=matrix,y=Cells,fill=condition)) +
  ylab("Number of cells in sample") + xlab("Matrix") + geom_col(position = "dodge") +
  scale_fill_manual(values= c("lightpink", "purple", "turquoise"))

small <- matrix_medium[,c("Experiment", "condition", "timepoint","matrix","type","CD45_positive_total_counts","CD45_negative_total_counts")]

small_long <- melt(setDT(small), id.vars = c("Experiment", "condition", "timepoint","matrix","type"), variable.name = "population")
small_long <- as.data.frame(small_long)
small_long$population <- gsub("_counts", "", small_long$population)
small_long$population <- factor(small_long$population, levels=c("CD45_positive_total", "CD45_negative_total"), labels=c("CD45+", "CD45-"))

ggplot(small_long,aes(x=matrix,y=value,fill=population)) +
  ylab("Cell count") + xlab("Matrix") + geom_col(position = "stack") + facet_grid(~condition) +
  scale_fill_manual(values= c("orange","grey")) +theme_bw()

#export to prism
prism <- matrix_medium[,c("condition", "matrix", "CD45_positive_total_counts", "CD45_negative_total_counts")]
prism$matrix <- ifelse(prism$matrix == "YES", "+", "-")
prism$ID <- paste0(prism$condition, " ", prism$matrix)
prism <- prism[,3:5]
rownames(prism) <- prism$ID
prism$ID <- NULL

write.xlsx(prism, "~/Documents/Analyses/SLICE/matrix_data.xlsx", rowNames = T)

#different way to export to prism
fragments_48 <- fragments_default[fragments_default$timepoint == 48,]
fragments_48 <- fragments_48[fragments_48$cryopreservation == "NO",]
fragments_48$condition <- factor(fragments_48$condition, levels = c("unstim", "aPD1", "aCD3"))
fragments_48 <- fragments_48[order(fragments_48$type, fragments_48$condition, fragments_48$Experiment),]
fragments_48$ID <- paste0(fragments_48$type, "_", fragments_48$condition)

CD4_CD39_CD137 <- fragments_48[,c("ID","Experiment","CD4+ CD39+ CD137+_percentage")]
CD4_CD39_CD137 <- as.data.frame(t(CD4_CD39_CD137))

CD8_CD39_CD137 <- fragments_48[,c("ID","Experiment","CD8+ CD39+ CD137+_percentage")]
CD8_CD39_CD137 <- as.data.frame(t(CD8_CD39_CD137))

CD4_CD25 <- fragments_48[,c("ID","Experiment","CD4+ CD25+_percentage")]
CD4_CD25 <- as.data.frame(t(CD4_CD25))

CD8_CD25 <- fragments_48[,c("ID","Experiment","CD8+ CD25+_percentage")]
CD8_CD25 <- as.data.frame(t(CD8_CD25))

CD4_CD137 <- fragments_48[,c("ID","Experiment","CD4+ CD137+_percentage")]
CD4_CD137 <- as.data.frame(t(CD4_CD137))

CD8_CD137 <- fragments_48[,c("ID","Experiment","CD8+ CD137+_percentage")]
CD8_CD137 <- as.data.frame(t(CD8_CD137))

myeloid <- fragments_48[,c("ID","Experiment","Myeloid_percentage")]
myeloid <- as.data.frame(t(myeloid))

test <- fragments_default[fragments_default$Experiment == "SLICE007",]
test_small <- test[,c("Experiment", "condition", "timepoint","matrix","type","B cells_percentage_of_CD45","CD4+ T cells_percentage_of_CD45",
                      "CD8+ T cells_percentage_of_CD45","CD56+ NK cells_percentage_of_CD45","CD56- myeloids_percentage_of_CD45")]
test_small$condition <- factor(test_small$condition, levels = c("unstim", "aPD1", "aCD3"))
test_small <- test_small[order(test_small$type, test_small$condition, test_small$timepoint),]
test_small$composition <- rowSums(test_small[,6:10])

colnames(test_small) <- c("Experiment", "condition", "timepoint","matrix","type","B cell","CD4+ T cell",
                     "CD8+ T cell","CD56+ NK cell","CD56- Myeloid", "composition")
test_small <- test_small[test_small$composition != 0,] #remove entries without CD45+ cells
test_small$Other <- 100-test_small$composition
test_small$Other[test_small$Other < 0] <- 0
test_small$composition <- NULL

write.xlsx(test_small, "~/Documents/Analyses/SLICE/cryopreservation_cell_fractions.xlsx")
