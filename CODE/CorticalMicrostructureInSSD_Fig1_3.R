###############################################################################
# Title: Early deviations from normative brain morphology and cortical microstructure in schizophrenia spectrum disorders
# Author: Claudio Aleman Morillo . Universidad de Sevilla
# Date: 2025
# Purpose
#   This script reproduces the regression analyses across multiple neuroimaging
#   feature sets (cortical thickness, cortical WM metrics, subcortical volumes,
#   and global tracts) and exports:
#     (1) A set of SVG figures (brain maps or bar plots, depending on feature type)
#     (2) A multi-sheet Excel file with the regression results for each analysis
#
# Inputs
#   - A single merged dataset in CSV format (user-provided)
#
# Outputs
#   - SVG figures saved to: output_dir
#   - Excel workbook saved to: output_dir
#

# Reproducibility notes
#   - The script assumes the dataset contains all referenced variables/columns.
#
# How to run
#   1) Adjust 'data_dir', 'output_dir', and 'data_file'
###############################################################################

rm(list = ls())  # Clean environment for reproducibility

# ------------------------------- Libraries -----------------------------------

  library(dplyr)
  library(tidyr)
  library(lmerTest)
  library(ggplot2)
  library(ggseg)
  library(pscl)
  library(grid)
  library(openxlsx)
# --------------------- Functions -------------------------------------------

# Regression Function  
generalRegresionCalculus <- function(D,formulaText,varNames,listVariables,listNames,methodL) {
  
  
  res <- data.frame(matrix(0, ncol = length(listVariables), nrow = length(varNames)))
  names(res) <- listNames
  row.names(res) <- varNames
  
  resPValues <- data.frame(matrix(0, ncol = length(listVariables), nrow = length(varNames)))
  names(resPValues) <- listNames
  row.names(resPValues) <- varNames
  
  resRPValues <- data.frame(matrix(0, ncol = length(listVariables), nrow = length(varNames)))
  names(resRPValues) <- listNames
  row.names(resRPValues) <- varNames
  
  
  for (i in 1:length(listVariables)) {
    
    
    nameZone<-listVariables[[i]]
    
    formulaTextMod <- gsub("\\$", nameZone, formulaText)
   
    dataComplete <- model.frame(formulaTextMod, data = D, na.action = na.omit)
    
    if(methodL==0){
      linealModel <- glm(formulaTextMod, family = gaussian(), data = dataComplete)
      
    }else{
      linealModel <- glm(formulaTextMod, family = poisson(), data = dataComplete)
    }
    print(summary(linealModel))
    
    
    res[, i] <- as.numeric((coef(summary(linealModel)))[2:(length(varNames) + 1), 3])
    resPValues[, i] <- as.numeric((coef(summary(linealModel)))[2:(length(varNames) + 1), 4])
    resRPValues[, i] <- as.numeric((coef(summary(linealModel)))[2:(length(varNames) + 1), 4])
  }
  
  # Transpose results for readability
  res <- t(res)
  resPValues <- t(resPValues)
  resRPValues<-t(resRPValues)
  
  # Apply False Discovery Rate (FDR) adjustment to p-values
  for (z in 1:ncol(resPValues)) {
    
    resPValues[, z] <- p.adjust(resPValues[, z], method = "fdr")
    
  }
  
  
  
  # Combine results into a single data frame
  resD <- data.frame(res)
  respD <- data.frame(resPValues)
  resRPD<-data.frame(resRPValues)
  #browser()
  colnames(respD) <- paste0(colnames(respD), '_SigValues')
  colnames(resRPD) <- paste0(colnames(resRPD), '_PValues')
  resD$namesArea <- rownames(resD)
  respD$namesArea <- rownames(respD)
  resRPD$namesArea <- rownames(resRPD)
  
  # Merge all results into a final data frame
  A <- merge(resD, respD, by = 'namesArea')
  A <- merge(A, resRPD, by = 'namesArea')
  return(A)
}

# Function to create brain region maps with significance representation
regional_brainmap_representation_borders <- function(data, title, sup_lim, inf_lim, midd_p) {
  
  # Define cortical zones
  zones <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal', 'cuneus', 'entorhinal',
             'frontal pole', 'fusiform', 'inferior parietal', 'inferior temporal', 'insula',
             'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal', 'lingual',
             'medial orbitofrontal', 'middle temporal', 'paracentral', 'parahippocampal',
             'pars opercularis', 'pars orbitalis', 'pars triangularis', 'pericalcarine',
             'postcentral', 'posterior cingulate', 'precentral', 'precuneus',
             'rostral anterior cingulate', 'rostral middle frontal', 'superior frontal',
             'superior parietal', 'superior temporal', 'supramarginal', 'temporal pole',
             'transverse temporal')
  
  # Adjust data format for visualization
  rowOrder <- c(paste0('lh_', gsub(' ', '', zones)), paste0('rh_', gsub(' ', '', zones)))
  
  if (nrow(data) == 34) {
    someData <- tibble(
      region = zones,
      values = data$values,
      significant = data$sig,
      groups = c(rep(title, nrow(data)))
    )
  } else if (nrow(data) == 68) {
    data <- data[rowOrder, ]
    someData <- tibble(
      label = c(paste0('lh_', gsub(' ', '', zones)), paste0('rh_', gsub(' ', '', zones))),
      values = data$values,
      significant = data$sig,
      groups = c(rep(title, nrow(data)))
    )
  }
  
  # Define border colors based on significance
  if (all(as.logical(someData$significant))) {borders <- rgb(0, 0, 0, alpha = 1)}
  else if (all(!as.logical(someData$significant))){
    borders <- rgb(0, 0, 0, alpha = 0.1)}
  else {borders <- c(rgb(0, 0, 0, alpha = 0.1), "black")}
  
  # Plot brain regions
  someData %>%
    group_by(groups) %>%
    ggplot() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = midd_p, limits = c(inf_lim, sup_lim)) +
    geom_brain(mapping = aes(fill = values, col = significant, alpha = significant, size = significant),
               atlas = dk, position = position_brain(hemi ~ side)) +
    scale_colour_manual(values = borders) +
    scale_alpha_manual(values = c(0.6, 1)) +
    scale_size_manual(values = c(0.9, 1)) +
    facet_wrap(~groups) + theme_brain2(plot.background = "white")
}


# ----------------------------- User parameters ------------------------------
# Set generic, portable paths (edit these to your local setup)
data_dir   <- "path/to/data"      # e.g., "data"
output_dir <- "path/to/output"    # e.g., "results/figures"

data_file <- file.path(data_dir, "Datafile.csv")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------- Data load ----------------------------------
DataBase <- read.csv(data_file, sep = ",", dec = ".", stringsAsFactors = FALSE)

# --------------------------- Variable definitions ---------------------------

zones <- c(
  "bankssts", "caudal anterior cingulate", "caudal middle frontal", "cuneus", "entorhinal",
  "frontal pole", "fusiform", "inferior parietal", "inferior temporal", "insula",
  "isthmus cingulate", "lateral occipital", "lateral orbitofrontal", "lingual",
  "medial orbitofrontal", "middle temporal", "paracentral", "parahippocampal",
  "pars opercularis", "pars orbitalis", "pars triangularis", "pericalcarine",
  "postcentral", "posterior cingulate", "precentral", "precuneus",
  "rostral anterior cingulate", "rostral middle frontal", "superior frontal",
  "superior parietal", "superior temporal", "supramarginal", "temporal pole",
  "transverse temporal"
)

#### Definitions of names of the necessary columns. They must be in DataBase

common_columns <- c(
  "Global_Cognitive_Functioning_average", "Case", "Sex",
  "Age_inclusion", "SAPS_Total", "SANS_Total", "BPRS_Total","EstimatedTotalIntraCranialVol","Assessment","Age_MRI"
)


# Column-name templates used in the dataset
WM_cortical_columns <- c(
  paste0("lh_", gsub(" ", "", zones)),
  paste0("rh_", gsub(" ", "", zones))
)

thickness_cortical_columns <- c(
  paste0("lh.CT.", gsub(" ", "", zones)),
  paste0("rh.CT.", gsub(" ", "", zones))
)

tract_Names <- c(
  "ACR","ACR.L","ACR.R","ALIC","ALIC.L","ALIC.R","BCC","CC","CGC","CGC.L","CGC.R","CGH",
  "CGH.L","CGH.R","CR","CR.L","CR.R","CST","CST.L","CST.R","EC","EC.L","EC.R","FX",
  "FX.ST.L","FX.ST.R","FXST","GCC","IC","IC.L","IC.R","IFO","IFO.L","IFO.R","PCR","PCR.L",
  "PCR.R","PLIC","PLIC.L","PLIC.R","PTR","PTR.L","PTR.R","RLIC","RLIC.L","RLIC.R","SCC",
  "SCR","SCR.L","SCR.R","SFO","SFO.L","SFO.R","SLF","SLF.L","SLF.R","SS","SS.L","SS.R",
  "UNC","UNC.L","UNC.R"
)

global_tract_Names <- c(
  "ACR","ALIC","BCC","CC","CGC","CGH","CR","CST","EC","FX","FXST","GCC","IC","IFO","PCR",
  "PLIC","PTR","RLIC","SCC","SCR","SFO","SLF","SS","UNC"
)

Volumes_Names <- c(
  "Left.Cerebellum.Cortex",
  "Left.Thalamus",
  "Left.Caudate",
  "Left.Putamen",
  "Left.Pallidum",
  "Left.Hippocampus",
  "Left.Amygdala",
  "Left.Accumbens.area",
  "Right.Cerebellum.White.Matter",
  "Right.Cerebellum.Cortex",
  "Right.Thalamus",
  "Right.Caudate",
  "Right.Putamen",
  "Right.Pallidum",
  "Right.Hippocampus",
  "Right.Amygdala",
  "Right.Accumbens.area",
  "Left.Lateral.Ventricle",
  "Left.Inf.Lat.Vent",
  "X3rd.Ventricle",
  "X4th.Ventricle"
)

# Names used for representation 
cortical_Names_Representation_list <- WM_cortical_columns

# ------------------------------- Preprocessing (if required) ------------------------------
# Keep baseline assessment only
DataBase <- DataBase %>% filter(Assessment == 1)

# Center/scale covariates (as in the original script)
DataBase$EstimatedTotalIntraCranialVol <- scale(DataBase$EstimatedTotalIntraCranialVol)
DataBase$Age_MRI <- scale(DataBase$Age_MRI, center = TRUE, scale = TRUE)
DataBase$Sex <- as.factor(DataBase$Sex)

# Residualized eTIV (adjusted by Age_MRI and Sex)
DataBase$residual_etiv <- residuals(
  lm(EstimatedTotalIntraCranialVol ~ Age_MRI + Sex, data = DataBase)
)

# Rescale symptom totals (these values must be integers)
DataBase$BPRS_Total <- DataBase$BPRS_Total / 0.25
DataBase$SAPS_Total <- DataBase$SAPS_Total / 0.25
DataBase$SANS_Total <- DataBase$SANS_Total / 0.25

# -------------------------- Analysis configuration --------------------------
# Analysis codes:
# 00 diagnosis ~ meanCT
# 01 diagnosis ~ regional CT
# 02 diagnosis ~ subcortical GM volumes
# 03 diagnosis ~ cortical WM
# 04 diagnosis ~ global tracts
# 10 cognition ~ meanCT
# 11 cognition ~ regional CT
# 12 cognition ~ subcortical GM volumes
# 13 cognition ~ cortical WM
# 14 cognition ~ global tracts
# 20 SAPS ~ meanCT
# 21 SAPS ~ regional CT
# 22 SAPS ~ subcortical GM volumes
# 23 SAPS ~ cortical WM
# 24 SAPS ~ global tracts
# 30 SANS ~ meanCT
# 31 SANS ~ regional CT
# 32 SANS ~ subcortical GM volumes
# 33 SANS ~ cortical WM
# 34 SANS ~ global tracts
# 40 BPRS ~ meanCT
# 41 BPRS ~ regional CT
# 42 BPRS ~ subcortical GM volumes
# 43 BPRS ~ cortical WM
# 44 BPRS ~ global tracts

valuesComb <- 0:4
combinations <- apply(expand.grid(valuesComb, valuesComb), 1, paste0, collapse = "")

# Storage for Excel output
results_list <- list()

# ---------------------- Helper: analysis selector ---------------------------
get_analysis_spec <- function(optionCalculus) {
  
  # Returns a list with:
  # nameExcel, VarY, VarX, formula, vars, listColumns, listColumnsNames,
  # typeRegression, typeRepresentation
  
  if (optionCalculus == "00") {
    list(
      nameExcel = "diagnosis ~ mean CT",
      VarY = "diagnosis",
      VarX = "Mean CT",
      formula = "meanCT2 ~ 1 + Case + Age_inclusion + Sex + residual_etiv + mean_euler",
      vars = c("Diagnosis", "Age", "Sex", "residual", "mean_euler"),
      listColumns = c("meanCT2"),
      listColumnsNames = c("meanCT"),
      typeRegression = 0,
      typeRepresentation = 0
    )
    
  } else if (optionCalculus == "01") {
    list(
      nameExcel = "diagnosis ~ CT",
      VarY = "diagnosis",
      VarX = "CT",
      formula = "$ ~ 1 + Case + Age_inclusion + Sex + residual_etiv + mean_euler",
      vars = c("Diagnosis", "Age", "Sex", "residual", "mean_euler"),
      listColumns = thickness_cortical_columns,
      listColumnsNames = cortical_Names_Representation_list,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "02") {
    list(
      nameExcel = "diagnosis ~ GM Volumes",
      VarY = "diagnosis",
      VarX = "SubCortical GM Volumes",
      formula = "$ ~ 1 + Case + Age_inclusion + Sex + residual_etiv + mean_euler",
      vars = c("Diagnosis", "Age", "Sex", "residual", "mean_euler"),
      listColumns = Volumes_Names,
      listColumnsNames = Volumes_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "03") {
    list(
      nameExcel = "diagnosis ~ WM subCortical",
      VarY = "diagnosis",
      VarX = "WM subCortical",
      formula = "$ ~ 1 + Case + Age_inclusion + Sex",
      vars = c("Diagnosis", "Age", "Sex"),
      listColumns = WM_cortical_columns,
      listColumnsNames = WM_cortical_columns,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "04") {
    list(
      nameExcel = "diagnosis ~ tracts",
      VarY = "diagnosis",
      VarX = "tracts",
      formula = "$ ~ 1 + Case + Age_inclusion + Sex",
      vars = c("Diagnosis", "Age", "Sex"),
      listColumns = global_tract_Names,
      listColumnsNames = global_tract_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "10") {
    list(
      nameExcel = "Cognition ~ MeanCT",
      VarY = "Cognition",
      VarX = "MeanCT",
      formula = "Global_Cognitive_Functioning_average ~ 1 + Case + Age_inclusion + Sex + residual_etiv + mean_euler + $ + Case:$",
      vars = c("Diagnosis", "Age", "Sex", "etiv", "mean_euler", "variable", "ss FEP"),
      listColumns = c("meanCT2"),
      listColumnsNames = c("meanCT"),
      typeRegression = 0,
      typeRepresentation = 0
    )
    
  } else if (optionCalculus == "11") {
    list(
      nameExcel = "Cognition ~ CT",
      VarY = "Cognition",
      VarX = "CT",
      formula = "Global_Cognitive_Functioning_average ~ 1 + Case + Age_inclusion + Sex + residual_etiv + mean_euler + $ + Case:$",
      vars = c("Diagnosis", "Age", "Sex", "etiv", "mean_euler", "variable", "ss FEP"),
      listColumns = thickness_cortical_columns,
      listColumnsNames = cortical_Names_Representation_list,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "12") {
    list(
      nameExcel = "Cognition ~ GM Volumes",
      VarY = "Cognition",
      VarX = "GM Volumes",
      formula = "Global_Cognitive_Functioning_average ~ 1 + Case + Age_inclusion + Sex + residual_etiv + mean_euler + $ + Case:$",
      vars = c("Diagnosis", "Age", "Sex", "etiv", "mean_euler", "variable", "ss FEP"),
      listColumns = Volumes_Names,
      listColumnsNames = Volumes_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "13") {
    list(
      nameExcel = "Cognition ~ WM cortical",
      VarY = "Cognition",
      VarX = "WM cortical",
      formula = "Global_Cognitive_Functioning_average ~ 1 + Case + Age_inclusion + Sex + $ + Case:$",
      vars = c("Diagnosis", "Age", "Sex", "variable", "ss FEP"),
      listColumns = WM_cortical_columns,
      listColumnsNames = WM_cortical_columns,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "14") {
    list(
      nameExcel = "Cognition ~ Tracts",
      VarY = "Cognition",
      VarX = "Tracts",
      formula = "Global_Cognitive_Functioning_average ~ 1 + Case + Age_inclusion + Sex + $ + Case:$",
      vars = c("Diagnosis", "Age", "Sex", "variable", "ss FEP"),
      listColumns = global_tract_Names,
      listColumnsNames = global_tract_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "20") {
    list(
      nameExcel = "SAPS ~ MeanCT",
      VarY = "SAPS",
      VarX = "MeanCT",
      formula = "SAPS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = c("meanCT2"),
      listColumnsNames = c("meanCT"),
      typeRegression = 0,
      typeRepresentation = 0
    )
    
  } else if (optionCalculus == "21") {
    list(
      nameExcel = "SAPS ~ CT",
      VarY = "SAPS",
      VarX = "CT",
      formula = "SAPS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = thickness_cortical_columns,
      listColumnsNames = cortical_Names_Representation_list,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "22") {
    list(
      nameExcel = "SAPS ~ GM Volumes",
      VarY = "SAPS",
      VarX = "GM Volumes",
      formula = "SAPS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = Volumes_Names,
      listColumnsNames = Volumes_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "23") {
    list(
      nameExcel = "SAPS ~ WM cortical",
      VarY = "SAPS",
      VarX = "WM cortical",
      formula = "SAPS_Total ~ 1 + Age_inclusion + Sex + $",
      vars = c("Age", "Sex", "variable"),
      listColumns = WM_cortical_columns,
      listColumnsNames = WM_cortical_columns,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "24") {
    list(
      nameExcel = "SAPS ~ Tracts",
      VarY = "SAPS",
      VarX = "Tracts",
      formula = "SAPS_Total ~ 1 + Age_inclusion + Sex + $",
      vars = c("Age", "Sex", "variable"),
      listColumns = global_tract_Names,
      listColumnsNames = global_tract_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "30") {
    list(
      nameExcel = "SANS ~ MeanCT",
      VarY = "SANS",
      VarX = "MeanCT",
      formula = "SANS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = c("meanCT2"),
      listColumnsNames = c("meanCT"),
      typeRegression = 0,
      typeRepresentation = 0
    )
    
  } else if (optionCalculus == "31") {
    list(
      nameExcel = "SANS ~ CT",
      VarY = "SANS",
      VarX = "CT",
      formula = "SANS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = thickness_cortical_columns,
      listColumnsNames = cortical_Names_Representation_list,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "32") {
    list(
      nameExcel = "SANS ~ GM Volumes",
      VarY = "SANS",
      VarX = "GM Volumes",
      formula = "SANS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = Volumes_Names,
      listColumnsNames = Volumes_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "33") {
    list(
      nameExcel = "SANS ~ WM cortical",
      VarY = "SANS",
      VarX = "WM cortical",
      formula = "SANS_Total ~ 1 + Age_inclusion + Sex + $",
      vars = c("Age", "Sex", "variable"),
      listColumns = WM_cortical_columns,
      listColumnsNames = WM_cortical_columns,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "34") {
    list(
      nameExcel = "SANS ~ Tracts",
      VarY = "SANS",
      VarX = "Tracts",
      formula = "SANS_Total ~ 1 + Age_inclusion + Sex + $",
      vars = c("Age", "Sex", "variable"),
      listColumns = global_tract_Names,
      listColumnsNames = global_tract_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "40") {
    list(
      nameExcel = "BPRS ~ MeanCT",
      VarY = "BPRS",
      VarX = "MeanCT",
      formula = "BPRS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = c("meanCT2"),
      listColumnsNames = c("meanCT"),
      typeRegression = 0,
      typeRepresentation = 0
    )
    
  } else if (optionCalculus == "41") {
    list(
      nameExcel = "BPRS ~ CT",
      VarY = "BPRS",
      VarX = "CT",
      formula = "BPRS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = thickness_cortical_columns,
      listColumnsNames = cortical_Names_Representation_list,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "42") {
    list(
      nameExcel = "BPRS ~ GM Volumes",
      VarY = "BPRS",
      VarX = "GM Volumes",
      formula = "BPRS_Total ~ 1 + Age_inclusion + Sex + residual_etiv + mean_euler + $",
      vars = c("Age", "Sex", "etiv", "mean_euler", "variable"),
      listColumns = Volumes_Names,
      listColumnsNames = Volumes_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else if (optionCalculus == "43") {
    list(
      nameExcel = "BPRS ~ WM cortical",
      VarY = "BPRS",
      VarX = "WM cortical",
      formula = "BPRS_Total ~ 1 + Age_inclusion + Sex + $",
      vars = c("Age", "Sex", "variable"),
      listColumns = WM_cortical_columns,
      listColumnsNames = WM_cortical_columns,
      typeRegression = 0,
      typeRepresentation = 1
    )
    
  } else if (optionCalculus == "44") {
    list(
      nameExcel = "BPRS ~ Tracts",
      VarY = "BPRS",
      VarX = "Tracts",
      formula = "BPRS_Total ~ 1 + Age_inclusion + Sex + $",
      vars = c("Age", "Sex", "variable"),
      listColumns = global_tract_Names,
      listColumnsNames = global_tract_Names,
      typeRegression = 0,
      typeRepresentation = 2
    )
    
  } else {
    return(NULL)
  }
}

# ------------------------------ Main loop -----------------------------------
for (optionCalculus in combinations) {
  
  spec <- get_analysis_spec(optionCalculus)
  if (is.null(spec)) {
    message("Skipping invalid option: ", optionCalculus)
    next
  }
  
  message("Running: ", spec$nameExcel)
  
  # Run regression (external function)
  Results <- generalRegresionCalculus(
    DataBase,
    spec$formula,
    spec$vars,
    spec$listColumns,
    spec$listColumnsNames,
    spec$typeRegression
  )
  
  # Store for Excel
  results_list[[spec$nameExcel]] <- Results
  
  # ------------------------ Figure export (if needed) ------------------------
  if (spec$typeRepresentation %in% c(1, 2)) {
    
    for (s in seq_along(spec$vars)) {
      
      # The Results object is expected to have:
      #   - a column "namesArea"
      #   - t-values and p-value flags arranged in blocks
      t_block_size <- (ncol(Results) - 1) / 3
      i <- s
      
      colNameT_value <- colnames(Results)[i + 1]
      colNameP_value <- colnames(Results)[i + t_block_size + 1]
      
      DataF <- data.frame(
        values = as.numeric(Results[, colNameT_value]),
        sig    = (Results[, colNameP_value] == 1)
      )
      
      rownames(DataF) <- Results$namesArea
      repLimits <- max(abs(DataF$values), na.rm = TRUE)
      
      # File name (portable)
      filename <- file.path(
        output_dir,
        paste0("Image_", spec$VarY, "_", spec$VarX, "__", colNameT_value, ".svg")
      )
      
      # Close any open plotting devices (safe export)
      while (dev.cur() > 1) dev.off()
      
      svg(filename, width = 7, height = 5)
      on.exit({
        dev.off()
      }, add = TRUE)
      
      par(mfrow = c(1, 1))
      plot.new()
      grid.newpage()
      
      if (spec$typeRepresentation == 1) {
        # Cortical brain map representation (external function)
        print(
          regional_brainmap_representation_borders(
            DataF,
            paste0(gsub("\\$", spec$VarX, spec$formula), " ## ", colNameT_value),
            repLimits,
            -repLimits,
            0
          )
        )
        
      } else if (spec$typeRepresentation == 2) {
        # Bar plot representation
        print(
          ggplot(DataF, aes(y = rownames(DataF), x = values, fill = sig)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
            scale_x_continuous(limits = c(-repLimits, repLimits)) +
            labs(
              title = gsub("\\$", spec$VarX, spec$formula),
              x = colNameT_value,
              y = spec$VarX
            ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 0),
              axis.text.y = element_text(size = 12, face = "bold")
            )
        )
      }
      
      message("Saved figure: ", filename)
    }
  } else {
    print(Results)
  }
}

# ------------------------------ Excel export --------------------------------
message("Creating Excel workbook...")

wb <- createWorkbook()
for (sheet_name in names(results_list)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, results_list[[sheet_name]])
}

excel_file <- file.path(output_dir, "Results_Figures_1_3.xlsx")
saveWorkbook(wb, excel_file, overwrite = TRUE)

message("Excel created: ", excel_file)
###############################################################################
