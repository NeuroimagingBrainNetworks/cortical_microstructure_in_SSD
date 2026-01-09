###############################################################################
# Title: Early deviations from normative brain morphology and cortical microstructure in schizophrenia spectrum disorders
# Author: Claudio Aleman Morillo . Universidad de Sevilla
# Date: 2025
#
# Purpose
#   This script runs three regression analyses used in Figure 2 and exports:
#     (1) An Excel workbook with three sheets of results
#     (2) Optional visualizations (brain map + bar plot) saved as SVG
#
# Analyses
#   A) AverageFA ~ Subcortical GM volumes
#   B) sGMV     ~ WM tracts
#   C) CT       ~ SWM (cortical thickness vs superficial white matter)
#
# Inputs
#   - A merged dataset in CSV format containing all required predictors/outcomes
#
# Outputs
#   - Excel workbook: paperAnalisisResults_Figure2.xlsx
#   - Optional bar plot: barPlotTracts_Figure2.svg
#
#
# How to run
#   1) Adjust 'data_dir' , 'output_dir' and 'DataFile.csv'

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

###############################################################

#----------------------------#
# 1. MAIN REGRESSION FUNCTIONS
#----------------------------#

# Function for general multiple regressions across a list of variables
generalRegresionCalculus <- function(D, formulaText, varNames, listVariables, listNames, methodL) {
  res <- resPValues <- resRPValues <- matrix(0, ncol = length(listVariables), nrow = length(varNames)) %>% as.data.frame()
  colnames(res) <- colnames(resPValues) <- colnames(resRPValues) <- listNames
  rownames(res) <- rownames(resPValues) <- rownames(resRPValues) <- varNames
  
  for (i in seq_along(listVariables)) {
    formulaTextMod <- gsub("\\$", listVariables[i], formulaText)
    dataComplete <- model.frame(formulaTextMod, data = D, na.action = na.omit)
    
    model <- if (methodL == 0) glm(formulaTextMod, family = gaussian(), data = dataComplete)
    else glm(formulaTextMod, family = poisson(), data = dataComplete)
    print(summary(model))
    coef_summary <- summary(model)$coefficients
    res[, i] <- round(coef_summary[2:(length(varNames) + 1), 3], 2)    # T-values
    resPValues[, i] <- coef_summary[2:(length(varNames) + 1), 4]       # Raw p-values
    resRPValues[, i] <- coef_summary[2:(length(varNames) + 1), 4]      # To adjust later
  }
  
  # Transpose results
  res <- t(res)
  resPValues <- t(resPValues)
  resRPValues <- t(resRPValues)
  
  for (z in seq_len(ncol(resPValues))) {
    resPValues[, z] <- round(resPValues[, z], 2)
    resRPValues[, z] <- round(p.adjust(resRPValues[, z], method = "fdr"), 3)
  }
  
  # Combine into final result
  resD <- data.frame(res, namesArea = rownames(res))
  respD <- data.frame(resPValues, namesArea = rownames(resPValues))
  resRPD <- data.frame(resRPValues, namesArea = rownames(resRPValues))
  colnames(respD)[-ncol(respD)] <- paste0(colnames(respD)[-ncol(respD)], "_Raw_pValues")
  colnames(resRPD)[-ncol(resRPD)] <- paste0(colnames(resRPD)[-ncol(resRPD)], "_Corrected_FDR_pValues")
  
  final <- merge(resD, respD, by = 'namesArea') %>%
    merge(resRPD, by = 'namesArea')
  return(final)
}

# Regression function for cortical thickness vs. superficial white matter
generalRegresionCalculus_WM_CT <- function(D, zones, varNames, methodL) {
  WM <- c(paste0('lh_', gsub(' ', '', zones)), paste0('rh_', gsub(' ', '', zones)))
  CT <- c(paste0('lh.CT.', gsub(' ', '', zones)), paste0('rh.CT.', gsub(' ', '', zones)))
  
  res <- resPValues <- resRPValues <- matrix(0, ncol = length(WM), nrow = length(varNames)) %>% as.data.frame()
  colnames(res) <- colnames(resPValues) <- colnames(resRPValues) <- WM
  rownames(res) <- rownames(resPValues) <- rownames(resRPValues) <- varNames
  
  for (i in seq_along(WM)) {
    formula <- '$Y ~ 1 + Case + Sex + Age_inclusion + residual_etiv + mean_euler+$X + $X:Case'
    formula_mod <- gsub("\\$Y", CT[i], formula)
    formula_mod <- gsub("\\$X", WM[i], formula_mod)
    dataComplete <- model.frame(formula_mod, data = D, na.action = na.omit)
    
    model <- if (methodL == 0) glm(formula_mod, family = gaussian(), data = dataComplete)
    else glm(formula_mod, family = poisson(), data = dataComplete)
    
    coef_summary <- summary(model)$coefficients
    res[, i] <- round(coef_summary[2:(length(varNames) + 1), 3], 2)
    resPValues[, i] <- coef_summary[2:(length(varNames) + 1), 4]
    resRPValues[, i] <- coef_summary[2:(length(varNames) + 1), 4]
  }
  
  # Transpose and correct p-values
  res <- t(res); resPValues <- t(resPValues); resRPValues <- t(resRPValues)
  for (z in seq_len(ncol(resPValues))) {
    resPValues[, z] <- round(resPValues[, z], 2)
    resRPValues[, z] <- round(p.adjust(resRPValues[, z], method = "fdr"), 3)
  }
  
  resD <- data.frame(res, namesArea = rownames(res))
  respD <- data.frame(resPValues, namesArea = rownames(resPValues))
  resRPD <- data.frame(resRPValues, namesArea = rownames(resRPValues))
  colnames(respD)[-ncol(respD)] <- paste0(colnames(respD)[-ncol(respD)], "_Raw_pValues")
  colnames(resRPD)[-ncol(resRPD)] <- paste0(colnames(resRPD)[-ncol(resRPD)], "_Corrected_FDR_pValues")
  
  final <- merge(resD, respD, by = 'namesArea') %>%
    merge(resRPD, by = 'namesArea')
  return(final)
}

#----------------------------#
# 1.2. REPRESENTATION FUNCTIONS
#----------------------------#

# Function: Plot brain regions with statistical significance borders using ggseg
regional_brainmap_representation_borders <- function(data, title, sup_lim, inf_lim, midd_p) {
  
  # Define Desikan-Killiany atlas regions
  zones <- c('bankssts', 'caudal anterior cingulate', 'caudal middle frontal', 'cuneus', 'entorhinal',
             'frontal pole', 'fusiform', 'inferior parietal', 'inferior temporal', 'insula',
             'isthmus cingulate', 'lateral occipital', 'lateral orbitofrontal', 'lingual',
             'medial orbitofrontal', 'middle temporal', 'paracentral', 'parahippocampal',
             'pars opercularis', 'pars orbitalis', 'pars triangularis', 'pericalcarine',
             'postcentral', 'posterior cingulate', 'precentral', 'precuneus',
             'rostral anterior cingulate', 'rostral middle frontal', 'superior frontal',
             'superior parietal', 'superior temporal', 'supramarginal', 'temporal pole',
             'transverse temporal')
  
  # Reorder if both hemispheres are included
  rowOrder <- c(paste0('lh_', gsub(' ', '', zones)), paste0('rh_', gsub(' ', '', zones)))
  
  if (nrow(data) == 34) {
    # One hemisphere only
    someData <- tibble(
      region = zones,
      values = data$values,
      significant = data$sig,
      groups = rep(title, nrow(data))
    )
  } else if (nrow(data) == 68) {
    # Both hemispheres
    data <- data[rowOrder, ]
    someData <- tibble(
      label = rowOrder,
      values = data$values,
      significant = data$sig,
      groups = rep(title, 68)
    )
  }
  
  # Define border transparency based on significance
  if (all(someData$significant)) {
    borders <- rgb(0, 0, 0, alpha = 1)
  } else if (all(!someData$significant)) {
    borders <- rgb(0, 0, 0, alpha = 0.1)
  } else {
    borders <- c(rgb(0, 0, 0, alpha = 0.1), "black")
  }
  
  # Plot using ggseg
  someData %>%
    group_by(groups) %>%
    ggplot() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = midd_p, limits = c(inf_lim, sup_lim)) +
    geom_brain(aes(fill = values, col = significant, alpha = significant, size = significant),
               atlas = dk, position = position_brain(hemi ~ side)) +
    scale_colour_manual(values = borders) +
    scale_alpha_manual(values = c(0.6, 1)) +
    scale_size_manual(values = c(0.9, 1)) +
    facet_wrap(~groups) +
    theme_brain2(plot.background = "white")
}


##################################################################
# ----------------------------- User parameters ------------------------------
# Generic, portable paths (edit to your setup)
data_dir   <- "path/to/data"      # e.g., "data"
output_dir <- "path/to/output"    # e.g., "results/figure2"

data_file <- file.path(data_dir, "DataFile.csv")

# Output files
excel_outfile <- file.path(output_dir, "paperAnalisisResults_Figure2.xlsx")
barplot_outfile <- file.path(output_dir, "barPlotTracts_Figure2.svg")

# Create output directory if missing
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ------------------------------- Data load ----------------------------------
DataBase <- read.csv(data_file, sep = ",", dec = ".", stringsAsFactors = FALSE)

# ----------------------- Required columns --------------

# required_commons_cols <- c("Assessment", "EstimatedTotalIntraCranialVol", "Age_MRI", "Sex",
#                    "Age_inclusion", "mean_euler", "Case","sGMV","EstimatedTotalIntraCranialVol")

# ------------------------------- Preprocessing (if required)  ------------------------------


# Baseline assessment only
DataBase <- DataBase %>% filter(Assessment == 1)

# Standardize covariates (as in original)
DataBase$EstimatedTotalIntraCranialVol <- scale(DataBase$EstimatedTotalIntraCranialVol)
DataBase$Age_MRI <- scale(DataBase$Age_MRI, center = TRUE, scale = TRUE)
DataBase$Sex <- as.factor(DataBase$Sex)

# Residualized eTIV (adjusted by Age_MRI and Sex)
DataBase$residual_etiv <- scale(residuals(
  lm(EstimatedTotalIntraCranialVol ~ Age_MRI + Sex, data = DataBase)
))



# ---------------------------- Variable/columns definitions ---------------------------
zones <- c(
  "bankssts", "caudal anterior cingulate", "caudal middle frontal", "cuneus", "entorhinal",
  "frontal pole", "fusiform", "inferior parietal", "inferior temporal", "insula", "isthmus cingulate",
  "lateral occipital", "lateral orbitofrontal", "lingual", "medial orbitofrontal", "middle temporal",
  "paracentral", "parahippocampal", "pars opercularis", "pars orbitalis", "pars triangularis",
  "pericalcarine", "postcentral", "posterior cingulate", "precentral", "precuneus",
  "rostral anterior cingulate", "rostral middle frontal", "superior frontal", "superior parietal",
  "superior temporal", "supramarginal", "temporal pole", "transverse temporal"
)

global_tract_Names <- c(
  "ACR","ALIC","BCC","CC","CGC","CGH","CR","CST","EC","FX","FXST","GCC","IC","IFO",
  "PCR","PLIC","PTR","RLIC","SCC","SCR","SFO","SLF","SS","UNC"
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


# ---------------------------- Regression analyses ---------------------------

# A) Average FA ~ Subcortical GM volumes
# Model term labels 
vars <- c("Diagnosis", "Sex", "Age", "etiv", "mean_euler", "varible", "variableFep")

formula1 <- "AverageFA ~ 1 + Case + Sex + Age_inclusion + residual_etiv + mean_euler + $ + $:Case"
ResultsC <- generalRegresionCalculus(
  DataBase,
  formula1,
  vars,
  Volumes_Names,
  Volumes_Names,
  typeRegression = 0
)

# B) sGMV ~ WM tracts
# Model term labels 
vars <- c("Diagnosis", "Sex", "Age", "etiv", "mean_euler", "varible", "variableFep")
formula2 <- "sGMV ~ 1 + Case + Sex + Age_inclusion + residual_etiv + mean_euler + $ + $:Case"
ResultsD <- generalRegresionCalculus(
  DataBase,
  formula2,
  vars,
  global_tract_Names,
  global_tract_Names,
  typeRegression = 0
)

# C) CT ~ SWM (Cortical thickness vs superficial white matter)
ResultsB <- generalRegresionCalculus_WM_CT(
  DataBase,
  zones,
  vars,
  typeRegression = 0
)

# ------------------------------- Export results ------------------------------
results_list <- list(
  "Average FA ~ GM volumes" = ResultsC,
  "sGMV ~ WM tracts"        = ResultsD,
  "CT ~ SWM"                = ResultsB
)

message("Writing Excel workbook: ", excel_outfile)

wb <- createWorkbook()
for (sheet in names(results_list)) {
  # Excel sheet names are limited to 31 characters
  safe_sheet <- substr(sheet, 1, 31)
  addWorksheet(wb, safe_sheet)
  writeData(wb, safe_sheet, results_list[[sheet]])
}
saveWorkbook(wb, excel_outfile, overwrite = TRUE)

# ------------------------------ Optional plots ------------------------------
# Toggle plotting (set to FALSE if you only want the Excel outputs)
make_plots <- TRUE

if (make_plots) {
  
  # ------------------------- Brain map representation -----------------------
  # Choose predictor index for visualization
  numberVar <- 7 # Select a column to representation (number of column)
  varName <- vars[numberVar]
  pVarName <- paste0(varName, "_Corrected_FDR_pValues")
  
  # Select result table (example: CT ~ SWM)
  Results <- ResultsB
  
  region_labels <- Results$namesArea
  is_significant <- Results[[pVarName]] <= 0.05
  effect_values <- Results[[varName]]
  
  dataRepresentation <- data.frame(values = effect_values, sig = is_significant)
  rownames(dataRepresentation) <- region_labels
  repLimits <- max(abs(dataRepresentation$values), na.rm = TRUE)
  
  # This produces a plot in the current graphics device
  regional_brainmap_representation_borders(
    data    = dataRepresentation,
    title   = paste0("CT ~ SWM: ", varName),
    sup_lim = repLimits,
    inf_lim = -repLimits,
    midd_p  = 0
  )
  
  # ------------------------------ Bar plot (SVG) ----------------------------
  # Example: sGMV ~ WM tracts table (ResultsD)
  Results <- ResultsD
  
  colNameT_value <- vars[7]  # variable of interest (as in original)
  colNameP_value <- paste0(colNameT_value, "_Corrected_FDR_pValues")
  
  DataF <- data.frame(
    values = Results[[colNameT_value]],
    sig    = Results[[colNameP_value]] < 0.05
  )
  rownames(DataF) <- Results$namesArea
  
  message("Saving bar plot: ", barplot_outfile)
  
  # Close any open devices except the null device
  while (dev.cur() > 1) dev.off()
  
  svg(barplot_outfile, width = 7, height = 5)
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  
  print(
    ggplot(DataF, aes(y = rownames(DataF), x = values, fill = sig)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
      scale_x_continuous(limits = c(-2.5, 2.5)) +
      labs(
        title = paste0("Effect of ", colNameT_value, " on WM tracts"),
        x = "T-value",
        y = "Region"
      ) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 16, face = "bold"))
  )
  
  dev.off()
}
###############################################################################
