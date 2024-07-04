##USEFUL SELF-DERIVED R FUNCTIONS
#Sidney van der Zande (last updated 12-06-2024)

#-----------
analyze_legendplex <- function(tumour_name, file_path, batch_mode = F, experiment = NA){
  suppressMessages(require(readxl))
  
  #read in the file
  suppressMessages(file <- as.data.frame(read_excel(file_path, sheet = 2)))
  
  #failsafe for incorrect column names/numbers
  if(ncol(file) != 28){
    stop("Incorrect excel template")
  } 
  
  
  
  #set column names
  colnames(file) <- c("well", "sample", "IL-5", "IL-13", "IL-2", "IL-6", "IL-9", "IL-10", "IFN-γ", "TNF-α",
                      "IL-17A", "IL-17F", "IL-4", "IL-22", "empty", "IL-8", "CXCL10", "CCL11", "CCL17", "CCL2", "CCL5",
                      "CCL3", "CXCL9", "CXCL5", "CCL20", "CXCL1", "CXCL11", "CCL4")
  
  #Formatting:
  file$empty <- NULL
  file <- apply(file, 2, function(y) gsub("NA", NA, y)) #set excluded samples to NA
  file <- as.data.frame(file)
  file <- file[!is.na(file$sample),]
  file[,-c(1:2)] <- apply(file[,-c(1:2)], 2, function(y) as.numeric(y))
  
  
  #failsafe if PDTF cytokines are missing
  if(!(all(c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20") %in% colnames(file)))){
    missing_col <- colnames(file)[!(c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20") %in% colnames(file))]
    stop(paste0("Error: missing the following cytokines: ", missing_col))
  }
  
  #Grab the samples that correspond to the tumour name
  if(any(grepl(tumour_name, file$sample))){
    sample <- file[grepl(tumour_name, file$sample),]
  } else{
    stop("That tumour was not found in this file.")
  }
  
  if(!is.na(experiment)){
    sample <- sample[grepl(experiment, sample$sample),]
  }
  
  if(batch_mode == F){
    print(paste0("Results for tumour ", tumour_name))
    print(paste0(nrow(sample), " samples were found for that tumour."))
  }
  
  #Split into unstim and aPD1
  unstim <- sample[grepl("unstim", sample$sample, ignore.case = T),]
  rownames(unstim) <- unstim$sample
  unstim$sample <- NULL
  unstim$well <- NULL
  
  PD1 <- sample[grepl("aPD1", sample$sample, ignore.case = T),]
  PD1 <- PD1[!grepl("IFN", PD1$sample, ignore.case = F),] #exclude the aPD1+aIFNyR1 condition
  if(nrow(PD1) != 0){
    rownames(PD1) <- PD1$sample
    PD1 <- PD1[,!(colnames(PD1) %in% c("sample", "well"))]
    PD1_present <- T
  } else{
    PD1_present <- F
  }
  
  PD1_IFNR <- sample[grepl("aPD1\\+aIFNyR1", sample$sample, ignore.case = F),]
  if(nrow(PD1_IFNR) != 0){
    rownames(PD1_IFNR) <- PD1_IFNR$sample
    PD1_IFNR <- PD1_IFNR[,!(colnames(PD1_IFNR) %in% c("sample", "well"))]
    PD1_IFNR_present <- T
  } else{
    PD1_IFNR_present <- F
  }
  
  hrIFNy <- sample[grepl("hrIFNy", sample$sample, ignore.case = F),]
  if(nrow(hrIFNy) != 0){
    rownames(hrIFNy) <- hrIFNy$sample
    hrIFNy <- hrIFNy[,!(colnames(hrIFNy) %in% c("sample", "well"))]
    hrIFNy_present <- T
  } else{
    hrIFNy_present <- F
  }
  
  extra_conditions <- F
  
  if(grepl("LP060", file_path)){ #in LP060, samples were included where RN440-15 and RN440-16 conditions exist
    RN440_15 <- sample[grepl("RN440-15", sample$sample, ignore.case = F),]
    RN440_16 <- sample[grepl("RN440-16", sample$sample, ignore.case = F),]
    
    rownames(RN440_15) <- RN440_15$sample
    RN440_15 <- RN440_15[,!(colnames(RN440_15) %in% c("sample", "well"))]
    
    rownames(RN440_16) <- RN440_16$sample
    RN440_16 <- RN440_16[,!(colnames(RN440_16) %in% c("sample", "well"))]
    extra_conditions <- T
  }
  
  #get the blanco values
  blank <- file[grepl("Standard 0pg/ml", file$sample, ignore.case = T),]
  blank <- blank[,!(colnames(blank) %in% c("sample", "well"))]
  
  #average the values
  blank_avg <- colMeans(blank)
  
  #save uncorrected values
  uncorrected <- rbind(unstim, PD1, PD1_IFNR, hrIFNy)
  if(extra_conditions == T){
    uncorrected <- rbind(uncorrected, RN440_15, RN440_16)
  }
  uncorrected <- rbind(uncorrected, blank_avg)
  rownames(uncorrected)[length(rownames(uncorrected))] <- "Blank average" #change the name of the blank
  
  #correct for blancos in raw values
  unstim <- sweep(unstim, 2, blank_avg, "-")
  if(PD1_present == T){
    suppressWarnings(PD1 <- sweep(PD1, 2, blank_avg, "-"))
  }
  if(PD1_IFNR_present == T){
    suppressWarnings(PD1_IFNR <- sweep(PD1_IFNR, 2, blank_avg, "-"))
  }
  if(hrIFNy_present == T){
    suppressWarnings(hrIFNy <- sweep(hrIFNy, 2, blank_avg, "-"))
  }
  
  if(extra_conditions == T){
    suppressWarnings(RN440_15 <- sweep(RN440_15, 2, blank_avg, "-"))
    suppressWarnings(RN440_16 <- sweep(RN440_16, 2, blank_avg, "-"))
    
    RN440_15[RN440_15 < 0] <- 0
    RN440_15_avg <- colMeans(RN440_15)
    
    RN440_16[RN440_16 < 0] <- 0
    RN440_16_avg <- colMeans(RN440_16)
    
  }
  #set all negative values to zero
  unstim[unstim < 0] <- 0
  PD1[PD1 < 0] <- 0
  PD1_IFNR[PD1_IFNR < 0] <- 0
  hrIFNy[hrIFNy < 0] <- 0
  
  #calculate averages for stim and unstim
  unstim_avg <- colMeans(unstim)
  PD1_avg <- colMeans(PD1)
  PD1_IFNR_avg <- colMeans(PD1_IFNR)
  hrIFNy_avg <- colMeans(hrIFNy)
  
  #calculate the PD1 delta (PD1 - unstim average)
  if(PD1_present == T){
    suppressWarnings(PD1_delta <-  sweep(PD1, 2, unstim_avg, "-"))
  }
  
  #calculate unstim deltas (unstim - unstim average)
  unstim_delta <- sweep(unstim, 2, unstim_avg, "-")
  
  #calculate the aPD1+aIFNyR1 deltas (aPD1+aIFNyR1 - unstim average)
  if(PD1_IFNR_present == T){
    suppressWarnings(PD1_IFNR_delta <- sweep(PD1_IFNR, 2, unstim_avg, "-"))
  }
  
  #calculate the hrIFNy deltas (hrIFNy - unstim average)
  if(hrIFNy_present == T){
    suppressWarnings(hrIFNy_delta <- sweep(hrIFNy, 2, unstim_avg, "-"))
  }
  
  #calculate the RN440-15 and RN440-16 (RN440-15/16 - unstim average)
  if(extra_conditions == T){
    suppressWarnings(RN440_15_delta <- sweep(RN440_15, 2, unstim_avg, "-"))
    suppressWarnings(RN440_16_delta <- sweep(RN440_16, 2, unstim_avg, "-"))
  }
  
  #calculate PDTF scores (only if PD1 is present, though)
  if(PD1_present == T){
    PDTF_matrix_PD1 <- PD1_delta[,c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")]
    PDTF_matrix_PD1[,1] <- ifelse(PDTF_matrix_PD1[,1] > 2.052, yes = 1, no = 0)
    PDTF_matrix_PD1[,2] <- ifelse(PDTF_matrix_PD1[,2] > 2.622, yes = 2, no = 0)
    PDTF_matrix_PD1[,3] <- ifelse(PDTF_matrix_PD1[,3] > 41.67, yes = 2, no = 0)
    PDTF_matrix_PD1[,4] <- ifelse(PDTF_matrix_PD1[,4] > 3.702, yes = 1, no = 0)
    PDTF_matrix_PD1[,5] <- ifelse(PDTF_matrix_PD1[,5] > 0.5597, yes = 2, no = 0)
    PDTF_matrix_PD1[,6] <- ifelse(PDTF_matrix_PD1[,6] > 395.3, yes = 1, no = 0)
    PDTF_matrix_PD1[,7] <- ifelse(PDTF_matrix_PD1[,7] > 482, yes = 2, no = 0)
    PDTF_matrix_PD1[,8] <- ifelse(PDTF_matrix_PD1[,8] > 0.3441, yes = 1, no = 0)
    PDTF_matrix_PD1[,9] <- ifelse(PDTF_matrix_PD1[,9] > 18.82, yes = 1, no = 0)
    PDTF_matrix_PD1[,10] <- ifelse(PDTF_matrix_PD1[,10] > 2.512, yes = 1, no = 0)
    PDTF_matrix_PD1[,11] <- ifelse(PDTF_matrix_PD1[,11] > 6.611, yes = 1, no = 0)
    
    #calculate the cumulative score, percentage and final responder status
    PDTF_matrix_PD1$total_score <- rowSums(PDTF_matrix_PD1)
    PDTF_matrix_PD1$percentage <- (PDTF_matrix_PD1$total_score/15)*100
    
    PDTF_matrix_PD1$PDTF_response <- ifelse(PDTF_matrix_PD1$percentage >= 40, yes = "PDTF-R", no = "PDTF-NR")
    PDTF_matrix_PD1$PDTF_response[PDTF_matrix_PD1$percentage >= 30 & PDTF_matrix_PD1$PDTF_response == "PDTF-NR"] <- "PDTF-BR"
    
    #calculate PDTF scores at baseline
    PDTF_matrix_unstim <- unstim_delta[,c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")]
    PDTF_matrix_unstim[,1] <- ifelse(PDTF_matrix_unstim[,1] > 2.052, yes = 1, no = 0)
    PDTF_matrix_unstim[,2] <- ifelse(PDTF_matrix_unstim[,2] > 2.622, yes = 2, no = 0)
    PDTF_matrix_unstim[,3] <- ifelse(PDTF_matrix_unstim[,3] > 41.67, yes = 2, no = 0)
    PDTF_matrix_unstim[,4] <- ifelse(PDTF_matrix_unstim[,4] > 3.702, yes = 1, no = 0)
    PDTF_matrix_unstim[,5] <- ifelse(PDTF_matrix_unstim[,5] > 0.5597, yes = 2, no = 0)
    PDTF_matrix_unstim[,6] <- ifelse(PDTF_matrix_unstim[,6] > 395.3, yes = 1, no = 0)
    PDTF_matrix_unstim[,7] <- ifelse(PDTF_matrix_unstim[,7] > 482, yes = 2, no = 0)
    PDTF_matrix_unstim[,8] <- ifelse(PDTF_matrix_unstim[,8] > 0.3441, yes = 1, no = 0)
    PDTF_matrix_unstim[,9] <- ifelse(PDTF_matrix_unstim[,9] > 18.82, yes = 1, no = 0)
    PDTF_matrix_unstim[,10] <- ifelse(PDTF_matrix_unstim[,10] > 2.512, yes = 1, no = 0)
    PDTF_matrix_unstim[,11] <- ifelse(PDTF_matrix_unstim[,11] > 6.611, yes = 1, no = 0)
    
    #calculate the cumulative score, percentage and final responder status at baseline
    PDTF_matrix_unstim$total_score <- rowSums(PDTF_matrix_unstim)
    PDTF_matrix_unstim$percentage <- (PDTF_matrix_unstim$total_score/15)*100
    
    PDTF_matrix_unstim$PDTF_response <- ifelse(PDTF_matrix_unstim$percentage >= 40, yes = "PDTF-R", no = "PDTF-NR")
    PDTF_matrix_unstim$PDTF_response[PDTF_matrix_unstim$percentage >= 30 & PDTF_matrix_unstim$PDTF_response == "PDTF-NR"] <- "PDTF-BR"
    
    #calculate overall deltas
    overall_delta_PD1<- PD1_avg - unstim_avg
    
    #only keep cytokines we need
    responder_delta <- overall_delta_PD1[names(overall_delta_PD1) %in% c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")]
    #correct order
    responder_delta <- responder_delta[order(factor(names(responder_delta), levels = c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")))]
    
    #calculate overall responder scores
    responder_delta[1] <- ifelse(responder_delta[1] > 2.052, yes = 1, no = 0)
    responder_delta[2] <- ifelse(responder_delta[2] > 2.622, yes = 2, no = 0)
    responder_delta[3] <- ifelse(responder_delta[3] > 41.67, yes = 2, no = 0)
    responder_delta[4] <- ifelse(responder_delta[4] > 3.702, yes = 1, no = 0)
    responder_delta[5] <- ifelse(responder_delta[5] > 0.5597, yes = 2, no = 0)
    responder_delta[6] <- ifelse(responder_delta[6] > 395.3, yes = 1, no = 0)
    responder_delta[7] <- ifelse(responder_delta[7] > 482, yes = 2, no = 0)
    responder_delta[8] <- ifelse(responder_delta[8] > 0.3441, yes = 1, no = 0)
    responder_delta[9] <- ifelse(responder_delta[9] > 18.82, yes = 1, no = 0)
    responder_delta[10] <- ifelse(responder_delta[10] > 2.512, yes = 1, no = 0)
    responder_delta[11] <- ifelse(responder_delta[11] > 6.611, yes = 1, no = 0)
    
    #Sum the scores, calculate percentages and give overall responder status
    responder <- ((sum(responder_delta))/15)*100
    
    if(responder >= 40){
      responder_output <- "PDTF-R"
    } else if(responder >= 30){
      responder_output <- "PDTF-BR"
    } else{
      responder_output <- "PDTF-NR"
    }
    
    responder <- round(responder, 1)
    
  }
  
  #calculate the other overall deltas
  if(PD1_IFNR_present == T){
    overall_delta_PD1_IFNR<- PD1_IFNR_avg - unstim_avg
  }
  if(hrIFNy_present == T){
    overall_delta_hrIFNy<- hrIFNy_avg - unstim_avg
  }
  
  if(extra_conditions == T){
    overall_delta_RN440_15 <- RN440_15_avg - unstim_avg
    overall_delta_RN440_16 <- RN440_16_avg - unstim_avg
  }
  
  #return table with corrected cytokine values
  suppressWarnings(corrected_values <- rbind(unstim,PD1,PD1_IFNR,hrIFNy,unstim_avg,PD1_avg,PD1_IFNR_avg,hrIFNy_avg))
  if(PD1_present == T & PD1_IFNR_present ==F & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,unstim_avg,PD1_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "aPD1 average"
    
    deltas <- rbind(unstim_delta, PD1_delta, overall_delta_PD1)
    rownames(deltas)[length(rownames(deltas))] <- "overall delta aPD1"
    
  } else if(PD1_present == T & PD1_IFNR_present ==T & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,PD1_IFNR,unstim_avg,PD1_avg,PD1_IFNR_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "aPD1_aIFNyR1 average"
    
    deltas <- rbind(unstim_delta, PD1_delta, PD1_IFNR_delta, overall_delta_PD1,overall_delta_PD1_IFNR)
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta aPD1_aIFNR1"
    
  } else if(PD1_present == T & PD1_IFNR_present ==F & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,hrIFNy,unstim_avg,PD1_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, PD1_delta, hrIFNy_delta, overall_delta_PD1,overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
    
  } else if(PD1_present == T & PD1_IFNR_present ==T & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,PD1_IFNR,hrIFNy,unstim_avg,PD1_avg,PD1_IFNR_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-3] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1_aIFNyR1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, PD1_delta, PD1_IFNR_delta, hrIFNy_delta, overall_delta_PD1,overall_delta_PD1_IFNR,overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))-2] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1_aIFNR1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
    
  } else if(PD1_present == F & PD1_IFNR_present ==F & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,unstim_avg)
    rownames(corrected_values)[length(rownames(corrected_values))] <- "unstim average"
    
    deltas <- rbind(unstim_delta)
    
  } else if(PD1_present == F & PD1_IFNR_present ==T & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1_IFNR,hrIFNy,unstim_avg,PD1_IFNR_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1_aIFNyR1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, PD1_IFNR_delta, hrIFNy_delta, overall_delta_PD1_IFNR,overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1_aIFNR1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
    
  } else if(PD1_present == F & PD1_IFNR_present ==T & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1_IFNR,unstim_avg,PD1_IFNR_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "aPD1_aIFNyR1 average"
    
    deltas <- rbind(unstim_delta, PD1_IFNR_delta, overall_delta_PD1_IFNR)
    rownames(deltas)[length(rownames(deltas))] <- "overall delta aPD1_aIFNR1"
    
  } else if(PD1_present == F & PD1_IFNR_present ==F & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,hrIFNy,unstim_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, hrIFNy_delta, overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
  }
  
  #add RN440-15 and RN440-16 if applicable
  if(extra_conditions == T){ #extra conditions only exist if other 4 conditions exist too
    corrected_values <- rbind(unstim,PD1,PD1_IFNR,hrIFNy,RN440_15,RN440_16,unstim_avg,PD1_avg,PD1_IFNR_avg,hrIFNy_avg,RN440_15_avg,RN440_16_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-5] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-4] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))-3] <- "aPD1_aIFNyR1 average"
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "hrIFNy average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "RN440-15 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "RN440-16 average"
    
    deltas <- rbind(unstim_delta, PD1_delta, PD1_IFNR_delta, hrIFNy_delta, RN440_15_delta, RN440_16_delta, overall_delta_PD1,overall_delta_PD1_IFNR,overall_delta_hrIFNy,overall_delta_RN440_15,overall_delta_RN440_16)
    rownames(deltas)[length(rownames(deltas))-4] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))-3] <- "overall delta aPD1_aIFNR1"
    rownames(deltas)[length(rownames(deltas))-2] <- "overall delta hrIFNy"
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta RN440-15"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta RN440-16"
    
  }
  
  
  #show responder status for individual fragments
  if(PD1_present == T){
    PDTF_matrix_full <- rbind(PDTF_matrix_unstim, PDTF_matrix_PD1)
    fragment_response <- data.frame(fragment = rownames(PDTF_matrix_full), response= PDTF_matrix_full$PDTF_response, percentage = PDTF_matrix_full$percentage)
    fragment_response$percentage <- round(fragment_response$percentage, 1)
    responder_total <- c(as.character(responder_delta), sum(responder_delta), as.character(responder), responder_output)
    PDTF_matrix_full <- rbind(PDTF_matrix_full, responder_total)
    rownames(PDTF_matrix_full)[length(rownames(PDTF_matrix_full))] <- "Overall"
    
    #produce output
    output <- list(uncorrected,corrected_values,deltas,PDTF_matrix_full)
    names(output) <- c("uncorrected values","corrected values","deltas","PDTF scores")
    
  } else{
    output <- list(uncorrected,corrected_values,deltas)
    names(output) <- c("uncorrected values","corrected values","deltas")
  }
  
  #Give output report:
  if(PD1_present == T){
    if(batch_mode == F){
      print("The fragments were classified as follows:")
      print(fragment_response)
    }
    
    print(paste0(tumour_name, " was overall classified as: ", responder_output, " (", responder, "%)"))
  }
  
  if(batch_mode == F){
    print("All values that were used for calculations can be found in the output of this function.")
  }
  
  return(output)
  
}

#--------

#Quality control for legendplexes based on fold-changes between conditions and unstim (concentrations)
legendplex_qc_concentration <- function(file_name){
  #read in data
  raw_mat <- read.xlsx(file_name, sheet = 2, colNames = T, rowNames = T) #sheet 2 contains corrected values
  
  #grab values for cytokines to be tested
  unstim <- raw_mat[grepl("unstim", rownames(raw_mat), ignore.case = T),]
  unstim <- unstim[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  unstim_average <- unlist(colMeans(unstim))
  unstim_average <- unlist(unstim_average)
  
  aPD1 <- raw_mat[grepl("aPD1 fragment", rownames(raw_mat), ignore.case = T),]
  aPD1 <- aPD1[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  aPD1_average <- unlist(colMeans(aPD1))
  aPD1_average <- unlist(aPD1_average)
  
  aPD1_aIFNR <- raw_mat[grepl("aPD1\\+aIFNyR1", rownames(raw_mat), ignore.case = T),]
  aPD1_aIFNR <- aPD1_aIFNR[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  aPD1_aIFNR_average <- unlist(colMeans(aPD1_aIFNR))
  aPD1_aIFNR_average <- unlist(aPD1_aIFNR_average)
  
  hrIFNy <- raw_mat[grepl("hrIFNy", rownames(raw_mat), ignore.case = T),]
  hrIFNy <- hrIFNy[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  hrIFNy_average <- unlist(colMeans(hrIFNy))
  hrIFNy_average <- unlist(hrIFNy_average)
  
  RN440_15 <- raw_mat[grepl("RN440-15", rownames(raw_mat), ignore.case = T),]
  RN440_15 <- RN440_15[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  RN440_15_average <- unlist(colMeans(RN440_15))
  RN440_15_average <- unlist(RN440_15_average)
  
  RN440_16 <- raw_mat[grepl("RN440-16", rownames(raw_mat), ignore.case = T),]
  RN440_16 <- RN440_16[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  RN440_16_average <- unlist(colMeans(RN440_16))
  RN440_16_average <- unlist(RN440_16_average)
  
  
  #calculate fold-changes
  FC_aPD1 <- aPD1_average/unstim_average
  names(FC_aPD1) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_aPD1_aIFNR <- aPD1_aIFNR_average/unstim_average
  names(FC_aPD1_aIFNR) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_hrIFNy <- hrIFNy_average/unstim_average
  names(FC_hrIFNy) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_RN440_15 <- RN440_15_average/unstim_average
  names(FC_RN440_15) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_RN440_16 <- RN440_16_average/unstim_average
  names(FC_RN440_16) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC<- list(FC_aPD1,FC_aPD1_aIFNR,FC_hrIFNy,FC_RN440_15,FC_RN440_16)
  names(FC) <- c("aPD1", "aPD1+aIFNyR1","hrIFNy","RN440-15","RN440-16")
  
  #for each fragment, extract the values for the five cytokines to be tested
  data_list <- list(unstim, aPD1, aPD1_aIFNR, hrIFNy,RN440_15,RN440_16)
  names(data_list) <- c("unstim", "aPD1", "aPD1+aIFNyR1", "hrIFNy","RN440-15","RN440-16")
  plot_table <- data.frame(matrix(nrow = 0, ncol=4))
  library(data.table)
  for(i in seq_along(data_list)){
    data_list[[i]]$name <- rownames(data_list[[i]])
    long <- melt(setDT(data_list[[i]]), id.vars = c("name"), variable.name = "factor")
    long$condition <- names(data_list)[[i]]
    plot_table <- rbind(plot_table, long)
  }
  
  #convert conditions to factors and split data frames based on the cytokine to be investigated
  plot_table$condition <- factor(plot_table$condition, levels = c("unstim", "aPD1", "aPD1+aIFNyR1", "hrIFNy","RN440-15","RN440-16"))
  
  plot_table_split <- split(plot_table, f = plot_table$factor)
  
  cols <- c("pink", "cyan", "darkblue", "gold", "red")
  names(cols) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  #plot values separately for each cytokine
  IL_6 <- ggplot(plot_table_split[["IL-6"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["IL-6"]]) + scale_color_manual(values = cols[["IL-6"]]) + ylab("concentration")
  IL_8 <- ggplot(plot_table_split[["IL-8"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["IL-8"]]) + scale_color_manual(values = cols[["IL-8"]]) + ylab("concentration")
  CCL2 <- ggplot(plot_table_split[["CCL2"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) + 
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["CCL2"]]) + scale_color_manual(values = cols[["CCL2"]]) + ylab("concentration")
  CXCL5 <- ggplot(plot_table_split[["CXCL5"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["CXCL5"]]) + scale_color_manual(values = cols[["CXCL5"]]) + ylab("concentration")
  CXCL1 <- ggplot(plot_table_split[["CXCL1"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["CXCL1"]]) + scale_color_manual(values = cols[["CXCL1"]]) + ylab("concentration")
  
  #combine into one plot
  require(cowplot)
  plots <- plot_grid(IL_6, IL_8, CCL2, CXCL5, CXCL1)
  plots
  
  quality_scores <- vector(mode = "list", length = length(FC))
  
  #calculate for each cytokine whether it is above/below the accepted fold-change threshold
  for (j in seq_along(FC)){
    for(r in FC[[j]]){
      if(r == "NaN"){
        r_score <- NA
      } else if (r <= 0.5 | r >= 2){
        r_score <- 1
      } else if (r <= 0.1 | r>= 10){
        r_score <- 2
      } else {
        r_score <- 0
      }
      quality_scores[[j]] <- c(quality_scores[[j]], r_score)
    }
  }
  
  names(quality_scores) <- names(FC)
  
  
  #for each condition, calculate whether the tumour passes QC or not
  #There are three rules for this:
  #1. None of the fold-changes can be higher than 10
  #2. From IL-6, IL-8 and CCL2, only two of them may have a fold-change higher than 2
  #3. If two out of three from IL-6, IL-8 and CCL2 have a fold-change higher than 2, only one out of CXCL5 and CXCL5 may have a fold-change higher than 2
  for(i in seq_along(quality_scores)){
    
    if(any(is.na(quality_scores[[i]]))){ #quality scores are NA if the condition does not exist
      next
    }
    names(quality_scores[[i]]) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
    stable_markers <- quality_scores[[i]][c("IL-6","IL-8","CCL2")]
    variable_markers <- quality_scores[[i]][c("CXCL5","CXCL1")]
    
    if(any(quality_scores[[i]] == 2)){ #if one marker has a 10-fold log change, fail quality check
      error <- quality_scores[quality_scores == 2]
      print(paste0("Quality check failed on ", names(quality_scores)[[i]], ": higher than 10-fold change for ", names(error)))
      
    } else if(sum(stable_markers) == 3){
      #if all 3 stable markers have a 2-fold change, fail quality check
      print(paste0("Quality check failed on ", names(quality_scores)[[i]], ": higher than 2-fold change for IL-6, IL-8 and CCL2"))
      
    } else if(sum(stable_markers) == 2 & sum(variable_markers) >= 1){ #if 2 stable markers have a 2-fold change and at least one variable marker has a 2-fold change, fail quality check
      error_stable <- stable_markers[stable_markers == 1]
      error_variable <- variable_markers[variable_markers == 1]
      names_stable <- names(error_stable)
      names_variable <- names(error_variable)
      
      print(paste0("Quality check failed on ", names(quality_scores)[[i]], ": higher than 2-fold change for ", names_stable, " and ", names_variable))
      
    } else{
      print(paste0("Quality check passed (", names(quality_scores)[[i]] ,")!"))
    }
    print(FC[[i]])
    
  }
  return(plots)
}
#---------------
#Function for performing quality control on Legendplexes based on ratio between unstim and all other conditions (MFI)
legendplex_qc_MFI <- function(MFI_file, experiment, tumour_name){
  raw_mat <- read.xlsx(MFI_file, sheet = 2, colNames = T, rowNames = F) 
  raw_mat <- raw_mat[!is.na(raw_mat$Sample),]
  
  
  colnames(raw_mat) <- c("well", "sample", "IL-5", "IL-13", "IL-2", "IL-6", "IL-9", "IL-10", "IFN-γ", "TNF-α",
                         "IL-17A", "IL-17F", "IL-4", "IL-22", "IL-8", "CXCL10", "CCL11", "CCL17", "CCL2", "CCL5",
                         "CCL3", "CXCL9", "CXCL5", "CCL20", "CXCL1", "CXCL11", "CCL4")
  
  #Formatting:
  raw_mat <- apply(raw_mat, 2, function(y) gsub("NA", NA, y)) #set excluded samples to NA
  raw_mat <- as.data.frame(raw_mat)
  raw_mat <- raw_mat[!is.na(raw_mat$sample),]
  raw_mat[,-c(1:2)] <- apply(raw_mat[,-c(1:2)], 2, function(y) as.numeric(y))
  
  raw_mat <- raw_mat[grepl(tumour_name,raw_mat$sample),]
  
  if(nrow(raw_mat) == 0){
    stop("Tumour not found in this Legendplex file")
  }
  
  if(tumour_name %in% c("LU035")){
    if(experiment == "TLS013"){
      raw_mat <- raw_mat[grepl("TLS013", raw_mat$sample),]
    } else if (experiment == "TLS015"){
      raw_mat <- raw_mat[grepl("TLS015", raw_mat$sample),]
    } #if the tumour was spread over multiple Legendplexes, this will not filter anything
  }
  
  if(tumour_name %in% c("LU052")){
    if(experiment == "TLS010"){
      raw_mat <- raw_mat[grepl("TLS010", raw_mat$sample),]
    } else if (experiment == "TLS014"){
      raw_mat <- raw_mat[grepl("TLS014", raw_mat$sample),]
    }
  }
  unstim <- raw_mat[grepl("unstim", raw_mat$sample, ignore.case = T),]
  unstim <- unstim[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  unstim_average <- unlist(colMeans(unstim))
  unstim_average <- unlist(unstim_average)
  
  aPD1 <- raw_mat[grepl("aPD1 fragment", raw_mat$sample, ignore.case = T),]
  aPD1 <- aPD1[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  aPD1_average <- unlist(colMeans(aPD1))
  aPD1_average <- unlist(aPD1_average)
  
  aPD1_aIFNR <- raw_mat[grepl("aPD1\\+aIFNyR1", raw_mat$sample, ignore.case = T),]
  aPD1_aIFNR <- aPD1_aIFNR[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  aPD1_aIFNR_average <- unlist(colMeans(aPD1_aIFNR))
  aPD1_aIFNR_average <- unlist(aPD1_aIFNR_average)
  
  hrIFNy <- raw_mat[grepl("hrIFNy", raw_mat$sample, ignore.case = T),]
  hrIFNy <- hrIFNy[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  hrIFNy_average <- unlist(colMeans(hrIFNy))
  hrIFNy_average <- unlist(hrIFNy_average)
  
  RN440_15 <- raw_mat[grepl("RN440-15", raw_mat$sample, ignore.case = T),]
  RN440_15 <- RN440_15[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  RN440_15_average <- unlist(colMeans(RN440_15))
  RN440_15_average <- unlist(RN440_15_average)
  
  RN440_16 <- raw_mat[grepl("RN440-16", raw_mat$sample, ignore.case = T),]
  RN440_16 <- RN440_16[,c("IL-6","IL-8","CCL2","CXCL5","CXCL1")]
  RN440_16_average <- unlist(colMeans(RN440_16))
  RN440_16_average <- unlist(RN440_16_average)
  
  FC_aPD1 <- aPD1_average/unstim_average
  names(FC_aPD1) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_aPD1_aIFNR <- aPD1_aIFNR_average/unstim_average
  names(FC_aPD1_aIFNR) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_hrIFNy <- hrIFNy_average/unstim_average
  names(FC_hrIFNy) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_RN440_15 <- RN440_15_average/unstim_average
  names(FC_RN440_15) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC_RN440_16 <- RN440_16_average/unstim_average
  names(FC_RN440_16) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  FC<- list(FC_aPD1,FC_aPD1_aIFNR,FC_hrIFNy,FC_RN440_15,FC_RN440_16)
  names(FC) <- c("aPD1", "aPD1+aIFNyR1","hrIFNy","RN440-15","RN440-16")
  
  data_list <- list(unstim, aPD1, aPD1_aIFNR, hrIFNy,RN440_15,RN440_16)
  names(data_list) <- c("unstim", "aPD1", "aPD1+aIFNyR1", "hrIFNy", "RN440-15","RN440-16")
  plot_table <- data.frame(matrix(nrow = 0, ncol=4))
  require(data.table)
  for(i in seq_along(data_list)){
    data_list[[i]]$name <- rownames(data_list[[i]])
    long <- melt(setDT(data_list[[i]]), id.vars = c("name"), variable.name = "factor")
    long$condition <- names(data_list)[[i]]
    plot_table <- rbind(plot_table, long)
  }
  
  plot_table$condition <- factor(plot_table$condition, levels = c("unstim", "aPD1", "aPD1+aIFNyR1", "hrIFNy","RN440-15","RN440-16"))
  
  plot_table_split <- split(plot_table, f = plot_table$factor)
  
  cols <- c("pink", "cyan", "darkblue", "gold", "red")
  names(cols) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
  
  IL_6 <- ggplot(plot_table_split[["IL-6"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["IL-6"]]) + scale_color_manual(values = cols[["IL-6"]]) + ylab("MFI")
  IL_8 <- ggplot(plot_table_split[["IL-8"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["IL-8"]]) + scale_color_manual(values = cols[["IL-8"]]) + ylab("MFI")
  CCL2 <- ggplot(plot_table_split[["CCL2"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["CCL2"]]) + scale_color_manual(values = cols[["CCL2"]]) + ylab("MFI")
  CXCL5 <- ggplot(plot_table_split[["CXCL5"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["CXCL5"]]) + scale_color_manual(values = cols[["CXCL5"]]) + ylab("MFI")
  CXCL1 <- ggplot(plot_table_split[["CXCL1"]], aes(x=condition, y = value, col = factor)) + geom_boxplot(position = "dodge", aes(col = factor)) +
    geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[["CXCL1"]]) + scale_color_manual(values = cols[["CXCL1"]]) + ylab("MFI")
  
  require(cowplot)
  plots <- plot_grid(IL_6, IL_8, CCL2, CXCL5, CXCL1)
  plots
  
  quality_scores <- vector(mode = "list", length = length(FC))
  
  for (j in seq_along(FC)){
    for(r in FC[[j]]){
      if(r == "NaN"){
        r_score <- NA
      } else if (r <= 0.5 | r >= 2){
        r_score <- 1
      } else if (r <= 0.1 | r>= 10){
        r_score <- 2
      } else {
        r_score <- 0
      }
      quality_scores[[j]] <- c(quality_scores[[j]], r_score)
    }
  }
  
  names(quality_scores) <- names(FC)
  
  
  for(i in seq_along(quality_scores)){
    
    if(any(is.na(quality_scores[[i]]))){
      next
    }
    names(quality_scores[[i]]) <- c("IL-6","IL-8","CCL2","CXCL5","CXCL1")
    stable_markers <- quality_scores[[i]][c("IL-6","IL-8","CCL2")]
    variable_markers <- quality_scores[[i]][c("CXCL5","CXCL1")]
    
    if(any(quality_scores[[i]] == 2)){ #if one marker has a 10-fold log change, fail quality check
      error <- quality_scores[quality_scores == 2]
      print(paste0(tumour_name, " Quality check failed on ", names(quality_scores)[[i]], ": higher than 10-fold change for ", names(error)))
      
    } else if(sum(stable_markers) == 3){
      #if all 3 stable markers have a 2-fold change, fail quality check
      print(paste0(tumour_name, " Quality check failed on ", names(quality_scores)[[i]], ": higher than 2-fold change for IL-6, IL-8 and CCL2"))
      
    } else if(sum(stable_markers) == 2 & sum(variable_markers) >= 1){ #if 2 stable markers have a 2-fold change and at least one variable marker has a 2-fold change, fail quality check
      error_stable <- stable_markers[stable_markers == 1]
      error_variable <- variable_markers[variable_markers == 1]
      names_stable <- names(error_stable)
      names_variable <- names(error_variable)
      
      print(paste0(tumour_name, " Quality check failed on ", names(quality_scores)[[i]], ": higher than 2-fold change for ", names_stable, " and ", names_variable))
      
    } else{
      print(paste0(tumour_name, " Quality check passed (", names(quality_scores)[[i]] ,")!"))
    }
    print(FC[[i]])
    
  }
  return(plots)
}
#----------------------
#v2 legendplex with normalisation
analyze_legendplex_v2 <- function(tumour_name, file_path, batch_mode = F, experiment = NA, normalisation = NA, LOD_file = NA){
  suppressMessages(require(readxl))
  
  #read in the file
  suppressMessages(file <- as.data.frame(read_excel(file_path, sheet = 2)))
  
  #failsafe for incorrect column names/numbers
  if(ncol(file) != 28){
    stop("Incorrect excel template")
  } 
  
  
  
  #set column names
  colnames(file) <- c("well", "sample", "IL-5", "IL-13", "IL-2", "IL-6", "IL-9", "IL-10", "IFN-γ", "TNF-α",
                      "IL-17A", "IL-17F", "IL-4", "IL-22", "empty", "IL-8", "CXCL10", "CCL11", "CCL17", "CCL2", "CCL5",
                      "CCL3", "CXCL9", "CXCL5", "CCL20", "CXCL1", "CXCL11", "CCL4")
  
  #Formatting:
  file$empty <- NULL
  file <- apply(file, 2, function(y) gsub("NA", NA, y)) #set excluded samples to NA
  file <- as.data.frame(file)
  file <- file[!is.na(file$sample),]
  file[,-c(1:2)] <- apply(file[,-c(1:2)], 2, function(y) as.numeric(y))
  
  
  #failsafe if PDTF cytokines are missing
  if(!(all(c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20") %in% colnames(file)))){
    missing_col <- colnames(file)[!(c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20") %in% colnames(file))]
    stop(paste0("Error: missing the following cytokines: ", missing_col))
  }
  
  file_raw <- file #save a copy of the non-normalised values
  
  #correct for blanks
  blank_average <- colSums(file[grepl("Standard 0pg/ml",file$sample),3:27]) #Blanks are in first two rows excluding the first two columns (sample information):take the average of the two
  file[,3:27] <- sweep(file[,3:27], 2, blank_average, "-") #subtract blanks from the data
  
  file_corrected <- file #save corrected but not-normalised data
  
  #setup for normalisation
  if(!is.na(normalisation)){
    if(is.na(LOD_file)){
      stop("Normalisation file not found.")
    }
    #These are the values of the highest standard, which is the highest detection limit
    upper <- c(9000,18000,5000,20000,16000,14000,15000,18000,6000,11000,9000,10000,16000,15000,14000,7000,14000,2000,31000,6000,6000,2000,5000,6000,3000)
    names(upper) <- c("IL-5", "IL-13", "IL-2", "IL-6", "IL-9", "IL-10", "IFN-γ", "TNF-α",
                      "IL-17A", "IL-17F", "IL-4", "IL-22", "IL-8", "CXCL10", "CCL11", "CCL17", "CCL2", "CCL5",
                      "CCL3", "CXCL9", "CXCL5", "CCL20", "CXCL1", "CXCL11", "CCL4")
    
    #The LODs are read in from a different document (values can be found on the first page of the individual Biolegend reports)
    suppressMessages(lower <- as.data.frame(read_excel(LOD_file, sheet = 1)))
    lower <- lower$LOD
    names(lower) <- names(upper)
    
    if(normalisation == "windsor"){
      for (i in 3:ncol(file)){ #skip the first two columns which contain sample information
        column <- file[,i]
        
        column_range <- column[column < upper[i-2] & column > lower[i-2]] #grab all samples that fall within the detection range (including reference files)
        #(i-2 because the columns start at index 3 and the thresholds start at index 1)
        column_range <- column_range[!is.na(column_range)] #preventing NA values from becoming the maximum/minimum

        maximum <- max(column_range)
        minimum <- min(column_range)
        
        column[column >= upper[i-2]] <- maximum #set everything above the detection limit to the highest value under the limit
        column[column <= lower[i-2]] <- minimum #set everything below the detection limit to the lowest value above the limit
        
        file[,i] <- column #replace column in original file with adjusted values
        
      }
    } else if (normalisation == "percentile"){
      for (i in 3:ncol(file)){ #skip the first two rows which contain sample information
        column <- file[,i]
        column_range <- column[column < upper[i-2] & column > lower[i-2]] #grab all samples that fall within the detection range (including reference files)
        #(i-2 because the columns start at index 3 and the thresholds start at index 1)
        column_range <- column_range[!is.na(column_range)] #preventing NA values from becoming the maximum/minimum
        
        blanks <- column[1:2] #save the blanks, because these will be subtracted later so I do not want them to be corrected
        
        percentile_99 <- quantile(column_range, c(.99))
        percentile_01 <- quantile(column_range, c(.01))
        
        column[column >= upper[i-2]] <- percentile_99 #set everything above the detection limit to the 99th percentile
        column[column <= lower[i-2]] <- percentile_01 #set everything below the detection limit to the first percentile
        column[1:2] <- blanks #restore the blanks to their original value
        
        file[,i] <- column #replace column in original file with adjusted values
        
      }
    }
    
  }
  
  #Grab the samples that correspond to the tumour name
  if(any(grepl(tumour_name, file$sample))){
    sample <- file[grepl(tumour_name, file$sample),]
    
    #also grab entries of that tumour in uncorrected and unnormalised data frames
    sample_raw <- file_raw[grepl(tumour_name, file_raw$sample),]
    sample_raw <- sample_raw[order(sample_raw$sample),]
    
    sample_corrected <- file_corrected[grepl(tumour_name, file_corrected$sample),]
    sample_corrected <- sample_corrected[order(sample_corrected$sample),]
    
  } else{
    stop("That tumour was not found in this file.")
  }
  
  if(!is.na(experiment)){
    sample <- sample[grepl(experiment, sample$sample),]
  }
  
  if(batch_mode == F){
    print(paste0("Results for tumour ", tumour_name))
    print(paste0(nrow(sample), " samples were found for that tumour."))
  }
  
  
  #Split into unstim and aPD1
  unstim <- sample[grepl("unstim", sample$sample, ignore.case = T),]
  rownames(unstim) <- unstim$sample
  unstim$sample <- NULL
  unstim$well <- NULL
  unstim[unstim < 0] <- 0 #set all negative values to zero
  
  
  PD1 <- sample[grepl("aPD1", sample$sample, ignore.case = T),]
  PD1 <- PD1[!grepl("IFN", PD1$sample, ignore.case = F),] #exclude the aPD1+aIFNyR1 condition
  if(nrow(PD1) != 0){
    rownames(PD1) <- PD1$sample
    PD1 <- PD1[,!(colnames(PD1) %in% c("sample", "well"))]
    PD1_present <- T
    PD1[PD1 < 0] <- 0 #set all negative values to zero
    
  } else{
    PD1_present <- F
  }
  
  PD1_IFNR <- sample[grepl("aPD1\\+aIFNyR1", sample$sample, ignore.case = F),]
  if(nrow(PD1_IFNR) != 0){
    rownames(PD1_IFNR) <- PD1_IFNR$sample
    PD1_IFNR <- PD1_IFNR[,!(colnames(PD1_IFNR) %in% c("sample", "well"))]
    PD1_IFNR_present <- T
    PD1_IFNR[PD1_IFNR < 0] <- 0
    
  } else{
    PD1_IFNR_present <- F
  }
  
  hrIFNy <- sample[grepl("hrIFNy", sample$sample, ignore.case = F),]
  if(nrow(hrIFNy) != 0){
    rownames(hrIFNy) <- hrIFNy$sample
    hrIFNy <- hrIFNy[,!(colnames(hrIFNy) %in% c("sample", "well"))]
    hrIFNy_present <- T
    hrIFNy[hrIFNy < 0] <- 0
    
  } else{
    hrIFNy_present <- F
  }
  
  extra_conditions <- F
  
  if(grepl("LP060", file_path)){ #in LP060, samples were included where RN440-15 and RN440-16 conditions exist
    RN440_15 <- sample[grepl("RN440-15", sample$sample, ignore.case = F),]
    RN440_16 <- sample[grepl("RN440-16", sample$sample, ignore.case = F),]
    
    rownames(RN440_15) <- RN440_15$sample
    RN440_15 <- RN440_15[,!(colnames(RN440_15) %in% c("sample", "well"))]
    RN440_15[RN440_15 < 0] <- 0
    
    rownames(RN440_16) <- RN440_16$sample
    RN440_16 <- RN440_16[,!(colnames(RN440_16) %in% c("sample", "well"))]
    RN440_16[RN440_16 < 0] <- 0
    
    extra_conditions <- T
  }
  
  
  #calculate averages for stim and unstim
  unstim_avg <- colMeans(unstim)
  PD1_avg <- colMeans(PD1)
  PD1_IFNR_avg <- colMeans(PD1_IFNR)
  hrIFNy_avg <- colMeans(hrIFNy)
  
  #calculate the PD1 delta (PD1 - unstim average)
  if(PD1_present == T){
    suppressWarnings(PD1_delta <-  sweep(PD1, 2, unstim_avg, "-"))
  }
  
  #calculate unstim deltas (unstim - unstim average)
  unstim_delta <- sweep(unstim, 2, unstim_avg, "-")
  
  #calculate the aPD1+aIFNyR1 deltas (aPD1+aIFNyR1 - unstim average)
  if(PD1_IFNR_present == T){
    suppressWarnings(PD1_IFNR_delta <- sweep(PD1_IFNR, 2, unstim_avg, "-"))
  }
  
  #calculate the hrIFNy deltas (hrIFNy - unstim average)
  if(hrIFNy_present == T){
    suppressWarnings(hrIFNy_delta <- sweep(hrIFNy, 2, unstim_avg, "-"))
  }
  
  #calculate the RN440-15 and RN440-16 (RN440-15/16 - unstim average)
  if(extra_conditions == T){
    RN440_15_avg <- colMeans(RN440_15)
    RN440_16_avg <- colMeans(RN440_16)
    
    suppressWarnings(RN440_15_delta <- sweep(RN440_15, 2, unstim_avg, "-"))
    suppressWarnings(RN440_16_delta <- sweep(RN440_16, 2, unstim_avg, "-"))
  }
  
  #calculate PDTF scores (only if PD1 is present, though)
  if(PD1_present == T){
    PDTF_matrix_PD1 <- PD1_delta[,c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")]
    PDTF_matrix_PD1[,1] <- ifelse(PDTF_matrix_PD1[,1] > 2.052, yes = 1, no = 0)
    PDTF_matrix_PD1[,2] <- ifelse(PDTF_matrix_PD1[,2] > 2.622, yes = 2, no = 0)
    PDTF_matrix_PD1[,3] <- ifelse(PDTF_matrix_PD1[,3] > 41.67, yes = 2, no = 0)
    PDTF_matrix_PD1[,4] <- ifelse(PDTF_matrix_PD1[,4] > 3.702, yes = 1, no = 0)
    PDTF_matrix_PD1[,5] <- ifelse(PDTF_matrix_PD1[,5] > 0.5597, yes = 2, no = 0)
    PDTF_matrix_PD1[,6] <- ifelse(PDTF_matrix_PD1[,6] > 395.3, yes = 1, no = 0)
    PDTF_matrix_PD1[,7] <- ifelse(PDTF_matrix_PD1[,7] > 482, yes = 2, no = 0)
    PDTF_matrix_PD1[,8] <- ifelse(PDTF_matrix_PD1[,8] > 0.3441, yes = 1, no = 0)
    PDTF_matrix_PD1[,9] <- ifelse(PDTF_matrix_PD1[,9] > 18.82, yes = 1, no = 0)
    PDTF_matrix_PD1[,10] <- ifelse(PDTF_matrix_PD1[,10] > 2.512, yes = 1, no = 0)
    PDTF_matrix_PD1[,11] <- ifelse(PDTF_matrix_PD1[,11] > 6.611, yes = 1, no = 0)
    
    #calculate the cumulative score, percentage and final responder status
    PDTF_matrix_PD1$total_score <- rowSums(PDTF_matrix_PD1)
    PDTF_matrix_PD1$percentage <- (PDTF_matrix_PD1$total_score/15)*100
    
    PDTF_matrix_PD1$PDTF_response <- ifelse(PDTF_matrix_PD1$percentage >= 40, yes = "PDTF-R", no = "PDTF-NR")
    PDTF_matrix_PD1$PDTF_response[PDTF_matrix_PD1$percentage >= 30 & PDTF_matrix_PD1$PDTF_response == "PDTF-NR"] <- "PDTF-BR"
    
    #calculate PDTF scores at baseline
    PDTF_matrix_unstim <- unstim_delta[,c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")]
    PDTF_matrix_unstim[,1] <- ifelse(PDTF_matrix_unstim[,1] > 2.052, yes = 1, no = 0)
    PDTF_matrix_unstim[,2] <- ifelse(PDTF_matrix_unstim[,2] > 2.622, yes = 2, no = 0)
    PDTF_matrix_unstim[,3] <- ifelse(PDTF_matrix_unstim[,3] > 41.67, yes = 2, no = 0)
    PDTF_matrix_unstim[,4] <- ifelse(PDTF_matrix_unstim[,4] > 3.702, yes = 1, no = 0)
    PDTF_matrix_unstim[,5] <- ifelse(PDTF_matrix_unstim[,5] > 0.5597, yes = 2, no = 0)
    PDTF_matrix_unstim[,6] <- ifelse(PDTF_matrix_unstim[,6] > 395.3, yes = 1, no = 0)
    PDTF_matrix_unstim[,7] <- ifelse(PDTF_matrix_unstim[,7] > 482, yes = 2, no = 0)
    PDTF_matrix_unstim[,8] <- ifelse(PDTF_matrix_unstim[,8] > 0.3441, yes = 1, no = 0)
    PDTF_matrix_unstim[,9] <- ifelse(PDTF_matrix_unstim[,9] > 18.82, yes = 1, no = 0)
    PDTF_matrix_unstim[,10] <- ifelse(PDTF_matrix_unstim[,10] > 2.512, yes = 1, no = 0)
    PDTF_matrix_unstim[,11] <- ifelse(PDTF_matrix_unstim[,11] > 6.611, yes = 1, no = 0)
    
    #calculate the cumulative score, percentage and final responder status at baseline
    PDTF_matrix_unstim$total_score <- rowSums(PDTF_matrix_unstim)
    PDTF_matrix_unstim$percentage <- (PDTF_matrix_unstim$total_score/15)*100
    
    PDTF_matrix_unstim$PDTF_response <- ifelse(PDTF_matrix_unstim$percentage >= 40, yes = "PDTF-R", no = "PDTF-NR")
    PDTF_matrix_unstim$PDTF_response[PDTF_matrix_unstim$percentage >= 30 & PDTF_matrix_unstim$PDTF_response == "PDTF-NR"] <- "PDTF-BR"
    
    #calculate overall deltas
    overall_delta_PD1<- PD1_avg - unstim_avg
    
    #only keep cytokines we need
    responder_delta <- overall_delta_PD1[names(overall_delta_PD1) %in% c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")]
    #correct order
    responder_delta <- responder_delta[order(factor(names(responder_delta), levels = c("IL-10", "IFN-γ", "CXCL10", "CCL5", "CCL17", "CXCL5", "CXCL1", "CXCL11", "CXCL9", "CCL4", "CCL20")))]
    
    #calculate overall responder scores
    responder_delta[1] <- ifelse(responder_delta[1] > 2.052, yes = 1, no = 0)
    responder_delta[2] <- ifelse(responder_delta[2] > 2.622, yes = 2, no = 0)
    responder_delta[3] <- ifelse(responder_delta[3] > 41.67, yes = 2, no = 0)
    responder_delta[4] <- ifelse(responder_delta[4] > 3.702, yes = 1, no = 0)
    responder_delta[5] <- ifelse(responder_delta[5] > 0.5597, yes = 2, no = 0)
    responder_delta[6] <- ifelse(responder_delta[6] > 395.3, yes = 1, no = 0)
    responder_delta[7] <- ifelse(responder_delta[7] > 482, yes = 2, no = 0)
    responder_delta[8] <- ifelse(responder_delta[8] > 0.3441, yes = 1, no = 0)
    responder_delta[9] <- ifelse(responder_delta[9] > 18.82, yes = 1, no = 0)
    responder_delta[10] <- ifelse(responder_delta[10] > 2.512, yes = 1, no = 0)
    responder_delta[11] <- ifelse(responder_delta[11] > 6.611, yes = 1, no = 0)
    
    #Sum the scores, calculate percentages and give overall responder status
    responder <- ((sum(responder_delta))/15)*100
    
    if(responder >= 40){
      responder_output <- "PDTF-R"
    } else if(responder >= 30){
      responder_output <- "PDTF-BR"
    } else{
      responder_output <- "PDTF-NR"
    }
    
    responder <- round(responder, 1)
    
  }
  
  #calculate the other overall deltas
  if(PD1_IFNR_present == T){
    overall_delta_PD1_IFNR<- PD1_IFNR_avg - unstim_avg
  }
  if(hrIFNy_present == T){
    overall_delta_hrIFNy<- hrIFNy_avg - unstim_avg
  }
  
  if(extra_conditions == T){
    overall_delta_RN440_15 <- RN440_15_avg - unstim_avg
    overall_delta_RN440_16 <- RN440_16_avg - unstim_avg
  }
  
  #return table with corrected cytokine values
  if(PD1_present == T & PD1_IFNR_present ==F & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,unstim_avg,PD1_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "aPD1 average"
    
    deltas <- rbind(unstim_delta, PD1_delta, overall_delta_PD1)
    rownames(deltas)[length(rownames(deltas))] <- "overall delta aPD1"
    
  } else if(PD1_present == T & PD1_IFNR_present ==T & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,PD1_IFNR,unstim_avg,PD1_avg,PD1_IFNR_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "aPD1_aIFNyR1 average"
    
    deltas <- rbind(unstim_delta, PD1_delta, PD1_IFNR_delta, overall_delta_PD1,overall_delta_PD1_IFNR)
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta aPD1_aIFNR1"
    
  } else if(PD1_present == T & PD1_IFNR_present ==F & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,hrIFNy,unstim_avg,PD1_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, PD1_delta, hrIFNy_delta, overall_delta_PD1,overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
    
  } else if(PD1_present == T & PD1_IFNR_present ==T & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1,PD1_IFNR,hrIFNy,unstim_avg,PD1_avg,PD1_IFNR_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-3] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1_aIFNyR1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, PD1_delta, PD1_IFNR_delta, hrIFNy_delta, overall_delta_PD1,overall_delta_PD1_IFNR,overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))-2] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1_aIFNR1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
    
  } else if(PD1_present == F & PD1_IFNR_present ==F & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,unstim_avg)
    rownames(corrected_values)[length(rownames(corrected_values))] <- "unstim average"
    
    deltas <- rbind(unstim_delta)
    
  } else if(PD1_present == F & PD1_IFNR_present ==T & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1_IFNR,hrIFNy,unstim_avg,PD1_IFNR_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "aPD1_aIFNyR1 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, PD1_IFNR_delta, hrIFNy_delta, overall_delta_PD1_IFNR,overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta aPD1_aIFNR1"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
    
  } else if(PD1_present == F & PD1_IFNR_present ==T & hrIFNy_present == F & extra_conditions == F){
    corrected_values <- rbind(unstim,PD1_IFNR,unstim_avg,PD1_IFNR_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "aPD1_aIFNyR1 average"
    
    deltas <- rbind(unstim_delta, PD1_IFNR_delta, overall_delta_PD1_IFNR)
    rownames(deltas)[length(rownames(deltas))] <- "overall delta aPD1_aIFNR1"
    
  } else if(PD1_present == F & PD1_IFNR_present ==F & hrIFNy_present == T & extra_conditions == F){
    corrected_values <- rbind(unstim,hrIFNy,unstim_avg,hrIFNy_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "hrIFNy average"
    
    deltas <- rbind(unstim_delta, hrIFNy_delta, overall_delta_hrIFNy)
    rownames(deltas)[length(rownames(deltas))] <- "overall delta hrIFNy"
  }
  
  #add RN440-15 and RN440-16 if applicable
  if(extra_conditions == T){ #extra conditions only exist if other 4 conditions exist too
    corrected_values <- rbind(unstim,PD1,PD1_IFNR,hrIFNy,RN440_15,RN440_16,unstim_avg,PD1_avg,PD1_IFNR_avg,hrIFNy_avg,RN440_15_avg,RN440_16_avg)
    rownames(corrected_values)[length(rownames(corrected_values))-5] <- "unstim average"
    rownames(corrected_values)[length(rownames(corrected_values))-4] <- "aPD1 average"
    rownames(corrected_values)[length(rownames(corrected_values))-3] <- "aPD1_aIFNyR1 average"
    rownames(corrected_values)[length(rownames(corrected_values))-2] <- "hrIFNy average"
    rownames(corrected_values)[length(rownames(corrected_values))-1] <- "RN440-15 average"
    rownames(corrected_values)[length(rownames(corrected_values))] <- "RN440-16 average"
    
    deltas <- rbind(unstim_delta, PD1_delta, PD1_IFNR_delta, hrIFNy_delta, RN440_15_delta, RN440_16_delta, overall_delta_PD1,overall_delta_PD1_IFNR,overall_delta_hrIFNy,overall_delta_RN440_15,overall_delta_RN440_16)
    rownames(deltas)[length(rownames(deltas))-4] <- "overall delta aPD1"
    rownames(deltas)[length(rownames(deltas))-3] <- "overall delta aPD1_aIFNR1"
    rownames(deltas)[length(rownames(deltas))-2] <- "overall delta hrIFNy"
    rownames(deltas)[length(rownames(deltas))-1] <- "overall delta RN440-15"
    rownames(deltas)[length(rownames(deltas))] <- "overall delta RN440-16"
    
  }
  
  
  #show responder status for individual fragments
  if(PD1_present == T){
    PDTF_matrix_full <- rbind(PDTF_matrix_unstim, PDTF_matrix_PD1)
    fragment_response <- data.frame(fragment = rownames(PDTF_matrix_full), response= PDTF_matrix_full$PDTF_response, percentage = PDTF_matrix_full$percentage)
    fragment_response$percentage <- round(fragment_response$percentage, 1)
    responder_total <- c(as.character(responder_delta), sum(responder_delta), as.character(responder), responder_output)
    PDTF_matrix_full <- rbind(PDTF_matrix_full, responder_total)
    rownames(PDTF_matrix_full)[length(rownames(PDTF_matrix_full))] <- "Overall"
    
    #produce output
    if(!is.na(normalisation)){
      output <- list(sample_raw,sample_corrected,corrected_values,deltas,PDTF_matrix_full)
      names(output) <- c("uncorrected values","corrected values","normalised values", "deltas","PDTF scores")
    } else{
      output <- list(sample_raw,sample_corrected,deltas,PDTF_matrix_full)
      names(output) <- c("uncorrected values","corrected values","deltas","PDTF scores")
    }

  } else{
    if(!is.na(normalisation)){
      output <- list(sample_raw,sample_corrected,corrected_values,deltas)
      names(output) <- c("uncorrected values","corrected values","normalised values","deltas")
    } else{
      output <- list(sample_raw,sample_corrected,deltas)
      names(output) <- c("uncorrected values","corrected values","deltas")
    }
  }
  
  #Give output report:
  if(PD1_present == T){
    if(batch_mode == F){
      print("The fragments were classified as follows:")
      print(fragment_response)
    }
    
    print(paste0(tumour_name, " was overall classified as: ", responder_output, " (", responder, "%)"))
  }
  
  if(batch_mode == F){
    print("All values that were used for calculations can be found in the output of this function.")
  }
  
  return(output)
  
}

#--------------------
#Function to visualise the ranges of cytokine/chemokine concentrations across different fragments/conditions
legendplex_fragment_control <- function(file_name, experiment){
  require(readxl)
  require(stringr)
  
  if(length(excel_sheets(file_name)) %in% c(4,3)){ #non-normalised values (with or without PDTF scores)
    raw_mat <- read.xlsx(file_name, sheet = 2, colNames = T, rowNames = T) #sheet 2 contains corrected values
  } else if(length(excel_sheets(file_name)) == 5){ #normalised values
    raw_mat <- read.xlsx(file_name, sheet = 3, colNames = T, rowNames = T) #sheet 3 contains corrected normalised values
  }
  
  #grab values for cytokines to be tested
  unstim <- raw_mat[grepl("unstim", rownames(raw_mat), ignore.case = T),]
  unstim_average <- unlist(colMeans(unstim))
  
  aPD1 <- raw_mat[grepl("aPD1 fragment", rownames(raw_mat), ignore.case = T),]
  aPD1_average <- unlist(colMeans(aPD1))
  
  aPD1_aIFNR <- raw_mat[grepl("aPD1\\+aIFNyR1", rownames(raw_mat), ignore.case = T),]
  aPD1_aIFNR_average <- unlist(colMeans(aPD1_aIFNR))
  
  hrIFNy <- raw_mat[grepl("hrIFNy", rownames(raw_mat), ignore.case = T),]
  hrIFNy_average <- unlist(colMeans(hrIFNy))
  hrIFNy_average <- unlist(hrIFNy_average)
  
  RN440_15 <- raw_mat[grepl("RN440-15", rownames(raw_mat), ignore.case = T),]
  RN440_15_average <- unlist(colMeans(RN440_15))
  
  RN440_16 <- raw_mat[grepl("RN440-16", rownames(raw_mat), ignore.case = T),]
  RN440_16_average <- unlist(colMeans(RN440_16))
  
  
  #calculate fold-changes
  FC_aPD1 <- aPD1_average/unstim_average
  names(FC_aPD1) <- colnames(aPD1)
  
  FC_aPD1_aIFNR <- aPD1_aIFNR_average/unstim_average
  names(FC_aPD1_aIFNR) <- colnames(aPD1_aIFNR)
  
  FC_hrIFNy <- hrIFNy_average/unstim_average
  names(FC_hrIFNy) <- colnames(hrIFNy)
  
  FC_RN440_15 <- RN440_15_average/unstim_average
  names(FC_RN440_15) <- colnames(RN440_15)
  
  FC_RN440_16 <- RN440_16_average/unstim_average
  names(FC_RN440_16) <- colnames(RN440_16)
  
  FC<- list(FC_aPD1,FC_aPD1_aIFNR,FC_hrIFNy,FC_RN440_15,FC_RN440_16)
  names(FC) <- c("aPD1", "aPD1+aIFNyR1","hrIFNy","RN440-15","RN440-16")
  
  #for each fragment, extract the values for the five cytokines to be tested
  data_list <- list(unstim, aPD1, aPD1_aIFNR, hrIFNy,RN440_15,RN440_16)
  names(data_list) <- c("unstim", "aPD1", "aPD1+aIFNyR1", "hrIFNy","RN440-15","RN440-16")
  plot_table <- data.frame(matrix(nrow = 0, ncol=4))
  library(data.table)
  for(i in seq_along(data_list)){
    data_list[[i]]$name <- rownames(data_list[[i]])
    long <- melt(setDT(data_list[[i]]), id.vars = c("name"), variable.name = "factor")
    long$condition <- names(data_list)[[i]]
    plot_table <- rbind(plot_table, long)
  }
  
  #convert conditions to factors and split data frames based on the cytokine to be investigated
  plot_table$condition <- factor(plot_table$condition, levels = c("unstim", "aPD1", "aPD1+aIFNyR1", "hrIFNy","RN440-15","RN440-16"))
  
  plot_table_split <- split(plot_table, f = plot_table$factor)
  
  cols <- c("#6A5F31","#922B3E","#20214F", "pink", "#4C2F27", "#26252D", "#82898F", "#C51D34", "#6C4675", "#646B63", "#308446", "#3B83BD",
            "#424632", "#F44611", "#2A6478", "#C1876B", "#025669", "#ED760E", "#C93C20", "#5B3A29", "#C6A664", "#354D73", "#6D3F5B", "#7E7B52", "#2F4538")
  names(cols) <- colnames(unstim)
  
  #plot values separately for each cytokine
  plot_list <- list()
  stats <- list()
  outlier_list <- list()
  
  for(i in names(cols)){
    plot <- ggplot(plot_table_split[[i]], aes(x=condition, y = value, col = cols[[i]])) + geom_boxplot(position = "dodge", aes(col = factor)) +
      geom_jitter(aes(col = factor)) + theme_minimal() + geom_hline(yintercept =  unstim_average[[i]]) + scale_color_manual(values = cols[[i]]) + ylab("concentration") +
      ggtitle(i) + guides(fill=guide_legend(title="mediator"))
    
    stats_i <- summary(plot_table_split[[i]]$value)
    
    test <- plot_table_split[[1]]$value
    z_scores <- (test - mean(test))/sd(test)
    
    outliers <- plot_table_split[[i]][z_scores > 3 | z_scores < -3,]
    
    plot_list[[(length(plot_list))+1]] <- plot
    names(plot_list)[[length(plot_list)]] <- i
    
    stats[[(length(stats))+1]] <- stats_i
    names(stats)[[length(stats)]] <- i
    
    outlier_list[[(length(outlier_list))+1]] <- outliers
    names(outlier_list)[[length(outlier_list)]] <- i
  }
  
  #get the tumour name from the file name
  tumour <- strsplit(file_name, "/")
  tumour <- tumour[[1]][[6]]
  tumour <- gsub(".xlsx", "", tumour)
  pdf(paste0("~/Documents/Analyses/TLS/QC plots/fragment_control_", tumour, "_", experiment, ".pdf"), width = 25,height = 25)
  require(ggpubr)
  require(gridExtra)
  suppressWarnings(grid.arrange(grobs = plot_list))
  dev.off()
  
  print(paste0("A file containing QC plots has been generated in the following location: ", paste0("~/Documents/Analyses/TLS/QC plots/fragment_control_", tumour, "_", experiment, ".pdf")))
  
  output <- list(stats,outlier_list)
  names(output) <- paste0(tumour,"_", experiment)
  
  return(output)
}


#-----------------
#Function to filter out outliers from individual data frames (based on outlier calculation over all fragments combined)
check_outliers_tumour <- function(data_frame, outliers, tumour, experiment){
  
  #check whether an entry with the same tumour/experiment name exists within the data frame containing outliers
  if(!is.na(experiment)){
    data_frame_outlier <- outliers[,grepl(paste0(tumour, "_", experiment), colnames(outliers))]
  } else {
    data_frame_outlier <- outliers[,grepl(tumour, colnames(outliers))]
  }
  
  #if outliers are found
  if(nrow(data_frame_outlier) != 0){
    names_mediators <- rownames(outliers)
    
    #check which conditions contain outliers
    unstim <- data_frame_outlier[,grepl("unstim", colnames(data_frame_outlier))]
    aPD1 <- data_frame_outlier[,grepl("aPD1 ", colnames(data_frame_outlier))]
    aPD1_aIFNR <- data_frame_outlier[,grepl("aPD1\\+aIFNyR1", colnames(data_frame_outlier))]
    hrIFNy <- data_frame_outlier[,grepl("hrIFNy", colnames(data_frame_outlier))]

    #starting data frames
    data_frame_filter <- data_frame
    data_frame_outlier <- as.data.frame(t(data_frame_outlier))
    
    if(length(unstim) != 0){ #if unstim contains an outlier
      names(unstim) <- names_mediators
      outliers_unstim <- names(unstim[is.na(unstim)]) #get the mediators for which there was an outlier
      
      data_frame_filter[grepl("unstim", rownames(data_frame_filter)),colnames(data_frame_filter) %in% outliers_unstim] <- NA #set the mediators to NA for all unstim fragments
      
    }
    if(length(aPD1) != 0){ #if aPD1 contains an outlier
      names(aPD1) <- names_mediators
      outliers_aPD1 <- names(aPD1[is.na(aPD1)]) #get the mediators for which there was an outlier
      
      data_frame_filter[grepl("aPD1 ", rownames(data_frame_filter)),colnames(data_frame_filter) %in% outliers_aPD1] <- NA #set the mediators to NA for all aPD1 fragments
      
    }
    
    if(length(aPD1_aIFNR) != 0){ #if aPD1 contains an outlier
      names(aPD1_aIFNR) <- names_mediators
      outliers_aPD1_aIFNR <- names(aPD1_aIFNR[is.na(aPD1_aIFNR)]) #get the mediators for which there was an outlier
      
      data_frame_filter[grepl("aPD1\\+aIFNyR1", rownames(data_frame_filter)),colnames(data_frame_filter) %in% outliers_aPD1_aIFNR] <- NA #set the mediators to NA for all aPD1+aIFNyR1 fragments
      
    }
    
    if(length(hrIFNy) != 0){ #if hrIFNy contains an outlier
      names(hrIFNy) <- names_mediators
      outliers_hrIFNy <- names(hrIFNy[is.na(hrIFNy)]) #get the mediators for which there was an outlier
      
      data_frame_filter[grepl("hrIFNy", rownames(data_frame_filter)),colnames(data_frame_filter) %in% outliers_hrIFNy] <- NA #set the mediators to NA for all hrIFNy fragments
      
    }
    
    
  } else {
    data_frame_filter <- data_frame #if no outliers were found, no correction will be done
  }
  
  return(data_frame_filter) #return the corrected object
}

#----------------
#Function to filter out outliers from individual data frames (based on outlier calculation over fragments individually)
check_outliers_fragment <- function(data_frame, outliers, tumour, experiment){
  
  #check whether an entry with the same tumour/experiment name exists within the data frame containing outliers
  if(!is.na(experiment)){
    data_frame_outlier <- outliers[,grepl(paste0(tumour, "_", experiment), colnames(outliers))]
    
    if(!is.data.frame(data_frame_outlier)){ #occurs if there is only one outlier for that tumour/experiment
      data_frame_outlier <- as.data.frame(data_frame_outlier)
      rownames(data_frame_outlier) <- colnames(data_frame)
      colnames(data_frame_outlier) <- colnames(outliers)[grepl(paste0(tumour, "_", experiment), colnames(outliers))]
    }
    
    colnames(data_frame_outlier) <- gsub(paste0("_", experiment), "", colnames(data_frame_outlier)) #remove experiment from name to match with original data
    
    #if the tumour is LU035 or LU052, these do have experiment names in their fragment names. Restore this:
    if(tumour %in% c("LU035","LU052")){
      colnames(data_frame_outlier) <- gsub(tumour, paste0(experiment, " ", tumour), colnames(data_frame_outlier)) #for example, "TLS013 LU035"
    }
    
  } else {
    data_frame_outlier <- outliers[,grepl(tumour, colnames(outliers))] #if no experiment name is in the fragment name, do not search for it
  }
  
  data_frame_filter <- data_frame
  
  #if outliers are found
  if(ncol(data_frame_outlier) != 0){
   for(i in 1:ncol(data_frame_outlier)){
     column <- data_frame_outlier[,i]
     names(column) <- rownames(data_frame_outlier)
     column_outlier <- names(column[is.na(column)]) #get the mediators for which there was an outlier in that fragment

     #filter the outlier values out
     name <- colnames(data_frame_outlier)[i]
     
     if(name %in% rownames(data_frame_filter)){ #only filter if the fragment is actually in the data frame
       data_frame_filter[name, column_outlier] <- NA
     }

   }
    
  } else {
    data_frame_filter <- data_frame #if no outliers were found, no correction will be done
  }
  
  return(data_frame_filter) #return the corrected object
}

#-------
#save pheatmaps (from https://gist.github.com/mathzero/a2070a24a6b418740c44a5c023f5c01e)
save_pheatmap <- function(x, filename, width=12, height=12){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  if(grepl(".png",filename)){
    png(filename, width=width, height=height, units = "in", res=300)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else if(grepl(".pdf",filename)){
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else{
    print("Filename did not contain '.png' or '.pdf'")
  }
}

#----
