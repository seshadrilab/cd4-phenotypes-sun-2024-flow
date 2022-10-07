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

# Perform marker scoring, selection, and threshold standardization 
# Review the diagnostic plots before proceeding to the discovery and annotation stage
system.time(
  generateAnnotationThresholds(
    gatingSet = gs,
    startingCellPop = "CD3+ Lymphocytes",
    activeChannels = faust_markers,
    projectPath = faust_path,
    seedValue = 123,
    threadNum = parallelly::availableCores(omit = 1)
  )
)

# Run FAUST (this took ~5 hours)
# First, run FAUST with annotationApproved set to FALSE to generate channel depth scores.
# Then, run FAUST with annotationsApproved set to TRUE and the revised depth score thresholds.
# TODO: Try modifying depthScoreThreshold and selectionQuantile 
# NOTE: phenotype discovery occurs in `discoverPhenotypes`. `faust` is a wrapper around `generateAnnotationThresholds` and `discoverPhenotypes`
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

# Extract annotated count matrix
countMatrix <- readRDS(file.path(faust_path,"faustData","faustCountMatrix.rds"))
pdata <- pData(gs) %>%
  rownames_to_column("sample")

# Set heatmap color scale
col_fun_prop <- circlize::colorRamp2(c(1, 1000, 170000), c("blue", "white", "red"))

# Set 0 values to NA
countMatrix_for_heatmap <- countMatrix
countMatrix_for_heatmap[countMatrix_for_heatmap == 0] <- NA

# Set column order
pdata <- pdata %>%
  arrange(Status, Stim)
col_order <- pdata$sample
countMatrix_for_heatmap <- t(countMatrix_for_heatmap)
countMatrix_for_heatmap <- countMatrix_for_heatmap[, col_order]

# Set column annotation
ann <- pdata %>%
  select(Status, Stim)
colors <- list("Status" = c("Pneg" = "#984EA3", "TST+" = "#4DAF4A"),
               "Stim" = c("DMSO" = "#bebada", "PP1" = "#8dd3c7", "TB WCL" = "#ffffb3"))
colAnn <- HeatmapAnnotation (df = ann,
                             which = "col",
                             col = colors,
                             annotation_width = unit(c(1,4), "cm"),
                             gap = unit(1, "mm"))

# Create complex heatmap plot
hmap <- Heatmap(countMatrix_for_heatmap,
                show_row_dend = FALSE, 
                show_column_dend = FALSE, 
                show_column_names = FALSE,
                cluster_columns = FALSE, 
                cluster_rows = TRUE,
                column_order = col_order,
                row_names_side = "left", 
                col = col_fun_prop, 
                na_col = "black", 
                top_annotation = colAnn,
                row_names_max_width = unit(15, "cm"),
                heatmap_legend_param = list(col_fun = col_fun_prop, title = "Count", 
                                            at = c(1, 100, 1000, 10000, 170000), break_dist = 1))

draw(hmap, heatmap_legend_side="right", annotation_legend_side="right")

# Save heatmap
cairo_pdf(file=here::here("out/FAUST_Th/FAUST_Th_count_heatmap.pdf"), width=15, height=7,
          onefile = TRUE, bg = "transparent", family = "Arial") 
print(draw(hmap, gap = unit(0.1, "in"), heatmap_legend_side="right", annotation_legend_side="right"))
dev.off()

# Generate bivariate dotplots showing the gating strategies to gate out discovered cell pops
# count_df <- as.data.frame(countMatrix) %>% 
#   mutate(sample = rownames(.)) %>%
#   left_join(pdata, by = "sample")
# count_df$sample_desc <- paste(count_df$`SAMPLE ID`, "_", count_df$Status, "_", count_df$Stim)
# count_df <- count_df %>%
#   select(-c(sample, Status, `WELL ID`, name, `EXPERIMENT NAME`, `SAMPLE ID`, Stim, `PLATE NAME`, `$DATE`)) %>%
#   column_to_rownames("sample_desc")

pops <- names(which(colSums(countMatrix) > 1000)) # At least 1000 cells across all samples
pops <- setdiff(pops,"0_0_0_0_0") # We don't want plots for the not-annotated cells.

for (col in pops) {
  for (r in rownames(countMatrix)) {
    faust:::plotFaustGates(col, r, faust_path)
  }
}