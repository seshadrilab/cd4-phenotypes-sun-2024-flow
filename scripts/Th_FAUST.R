library(tidyverse)
library(flowWorkspace)
library(flowCore)
library(here)
library(faust)
library(ggplot2)
library(ggdendro)

# Load gating set 
gs <- load_gs(here::here("out/GatingSets/RSTR_Th_GatingSet"))

# Set FAUST directory
if(!dir.exists(here::here("out/FAUST_Th"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/FAUST_Th")))
  dir.create(here::here("out/FAUST_Th"), recursive = T)
}

faust_path <- here::here("out/FAUST_Th")

# Set markers for FAUST analysis
faust_markers <- c("CTLA4", "CD4","OX40", "CD154", "CXCR3", "TBET", "CCR6", "IL17a", 
                   "RORyT", "CCR7", "CD137", "CD8a", "IFNg", "CD45RA")

# Run FAUST (this took ~5 hours)
# First, run FAUST with annotationApproved set to FALSE to generate channel depth scores.
# Then, run FAUST with annotationsApproved set to TRUE and the revised depth score thresholds.
# TODO: Try modifying depthScoreThreshold and selectionQuantile 
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

# # Examine output (annotated count matrix)
# count_df <- as.data.frame(readRDS(file.path(faust_path,"faustData","faustCountMatrix.rds")))
# count_df <- merge(pData(gs), count_df, by = 0)

# Extract annotated count matrix
countMatrix <- readRDS(file.path(faust_path,"faustData","faustCountMatrix.rds"))
countMatrix <- merge(pData(gs), countMatrix, by = 0)

# Make suitable for ggplot
# TODO: Add pData info 
count_long <- countMatrix %>%
  as.data.frame %>%
  mutate(sample = rownames(.)) %>%
  
  gather(key = population, value = count, 1:(ncol(.) - 1) )

# Run clustering
count_dendro <- as.dendrogram(hclust(d = dist(x = t(countMatrix))))

# Create dendrogram plot
dendro.plot <- ggdendrogram(data = count_dendro, rotate = TRUE) +
  theme(axis.text.y = element_text(size = 6))

# Extract the order
count_order <- order.dendrogram(count_dendro)

# Order the levels 
count_long$population <- factor(x = count_long$population,
                                levels = rownames(t(countMatrix))[count_order],
                                ordered = TRUE)

# Create heatmap plot
ggplot(data = count_long, aes(x = sample, y = population)) +
  geom_tile(aes(fill = count)) +
  scale_fill_distiller(trans ="log10", palette = 2, type = "div") +
  theme(legend.position = "top") +
  ggtitle("Cell counts of populations by samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Generate bivariate dotplots showing the gating strategie to gate out discovered cell pops
# TODO: change the rownames of countMatrix so the file names are more informative
pops <- names(which(colSums(countMatrix) > 1000)) # at least 1000 cells across all samples
pops <- setdiff(pops,"0_0_0_0_0") #We don't want plots for the not-annotated cells.
for (col in pops) {
  for (r in rownames(countMatrix)) {
    faust:::plotFaustGates(col, r, faust_path)
  }
}
