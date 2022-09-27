library(tidyverse)
library(flowWorkspace)
library(flowCore)
library(here)
library(faust)
library(ggplot2)
library(ggdendro)
library(ComplexHeatmap)

# Load gating set 
gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet"))

# Set FAUST directory
if(!dir.exists(here::here("out/FAUST_Treg"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/FAUST_Treg")))
  dir.create(here::here("out/FAUST_Treg"), recursive = T)
}

faust_path <- here::here("out/FAUST_Treg")

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

# # Examine output (annotated count matrix)
# count_df <- as.data.frame(readRDS(file.path(faust_path,"faustData","faustCountMatrix.rds")))
# count_df <- merge(pData(gs), count_df, by = 0)

# Extract annotated count matrix
countMatrix <- readRDS(file.path(faust_path,"faustData","faustCountMatrix.rds"))
pdata <- pData(gs) %>%
  rownames_to_column("sample")

# Make suitable for ggplot
count_long <- countMatrix %>%
  as.data.frame %>%
  mutate(sample = rownames(.)) %>%
  left_join(pdata, by = "sample") %>%
  gather(key = population, value = count, 1:(ncol(.) - 9))

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

# Create complex heatmap
# TODO: Change color scale
countMatrix_log10 <- log10(countMatrix)
countMatrix_log10[countMatrix_log10 == -Inf] <- 0
col_fun <- circlize::colorRamp2(c(0, 80000, 160100), c("blue", "white", "red"))
heatmap <- Heatmap(t(countMatrix), name = "Cell Counts by Sample",
                   show_row_dend = FALSE, show_column_dend = FALSE,
                   row_names_side = "left", col = col_fun)
heatmap

# Generate bivariate dotplots showing the gating strategies to gate out discovered cell pops
count_df <- as.data.frame(countMatrix) %>% 
  mutate(sample = rownames(.)) %>%
  left_join(pdata, by = "sample")
count_df$sample_desc <- paste(count_df$`SAMPLE ID`, "_", count_df$Status, "_", count_df$Stim)
count_df <- count_df %>%
  select(-c(sample, Status, `WELL ID`, name, `EXPERIMENT NAME`, `SAMPLE ID`, Stim, `PLATE NAME`, `$DATE`)) %>%
  column_to_rownames("sample_desc")

pops <- names(which(colSums(count_df) > 1000)) # At least 1000 cells across all samples
pops <- setdiff(pops,"0_0_0_0_0") # We don't want plots for the not-annotated cells.

for (col in pops) {
  for (r in rownames(count_df)) {
    faust:::plotFaustGates(col, r, faust_path)
  }
}
