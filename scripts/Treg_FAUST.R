library(tidyverse)
library(flowWorkspace)
library(flowCore)
library(here)
library(faust)
library(ggplot2)

# Load gating set 
gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet"))

# Set FAUST directory
if(!dir.exists(here::here("out/FAUST"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/FAUST")))
  dir.create(here::here("out/FAUST"), recursive = T)
}

faust_path <- here::here("out/FAUST")

# Set markers for FAUST analysis
faust_markers <- c("CTLA4", "CD4","OX40", "CD154", "IL10", "CD39", "FOXP3", 
                   "CCR7", "CD137", "CD8a", "CD25", "CD73")

# Run FAUST (this took ~4 hours)
system.time(
  faust(
    gatingSet = gs,
    startingCellPop = "CD3+ Lymphocytes",
    activeChannels = faust_markers,
    projectPath = faust_path,
    annotationsApproved = TRUE,
    seedValue = 123,
    threadNum = parallelly::availableCores(omit = 1),
    debugFlag = TRUE
  )
)

# Examine output (annotated count matrix)
count_df <- as.data.frame(readRDS(file.path(faust_path,"faustData","faustCountMatrix.rds")))
count_df <- merge(pData(gs), count_df, by = 0)
