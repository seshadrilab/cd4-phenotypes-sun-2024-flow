library(tidyverse)
library(flowWorkspace)
library(flowCore)
library(here)
library(faust)
library(ggplot2)
library(ggdendro)
library(ComplexHeatmap)
library(stringr)

# Load gating set 
gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet"))

# Set FAUST directory
if(!dir.exists(here::here("out/FAUST_Treg_V2"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/FAUST_Treg_V2")))
  dir.create(here::here("out/FAUST_Treg_V2"), recursive = T)
}

faust_path <- here::here("out/FAUST_Treg_V2")

# Set variables for FAUST analysis
faust_markers <- c("CTLA4", "CD4","OX40", "CD154", "IL10", "CD39", "FOXP3", 
                   "CCR7", "CD137", "CD8a", "CD25", "CD73")
eu <- "SAMPLE ID"
start_node <- "CD3+ Lymphocytes"

# Perform marker scoring, selection, and threshold standardization (for me, this took ~1.5 hours) 
# Review the diagnostic plots before proceeding to the discovery and annotation stage
system.time(
  faust(
    gatingSet = gs,
    experimentalUnit = eu,
    startingCellPop = start_node,
    activeChannels = faust_markers,
    projectPath = faust_path,
    depthScoreThreshold = 0.01, # default
    selectionQuantile = 0.5, # default
    annotationsApproved = FALSE, # set to FALSE before we inspect the score line plots
    seedValue = 123,
    threadNum = parallelly::availableCores(omit = 1),
    debugFlag = TRUE
  )
)

# # Link scoreLines plot to the experimental unit data
# all_eus <- list.files(file.path(faust_path,"faustData","expUnitData"))
# firstEu <- TRUE
# for (eu in all_eus) {
#   eu_data <- as.data.frame(readRDS(file.path(faust_path,"faustData","expUnitData",eu,"expUnitExprs.rds")))
#   eu_data$`eu_id` <- eu
#   if (firstEu) {
#     all_eu_data <- eu_data
#     firstEu <- FALSE
#   }
#   else {
#     all_eu_data <- rbind(all_eu_data,eu_data)
#   }
# }
# depth_score_demo_df <- gather(all_eu_data,key="Marker",value="Expression",-eu_id)
# p1 <- ggplot(depth_score_demo_df, aes(x = Expression, y = eu_id,fill=Marker))+
#   facet_wrap(~Marker,ncol=2)+
#   geom_density_ridges(rel_min_height=0.01)+
#   theme_ridges(center_axis_label=TRUE)+
#   ylab("Simulated experimental unit")+
#   xlab("Simulated expression")+
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.position="bottom")
# p1

# Set new FAUST directory
if(!dir.exists(here::here("out/FAUST_Treg_V3"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/FAUST_Treg_V3")))
  dir.create(here::here("out/FAUST_Treg_V3"), recursive = T)
}

faust_path <- here::here("out/FAUST_Treg_V3")

# Perform cell discovery with the the revised depth score thresholds (this took ~5 hours)
# Running FAUST with annotationsApproved set to TRUE signals that tuning is complete
system.time(
  faust(
    gatingSet = gs,
    experimentalUnit = eu,
    startingCellPop = start_node,
    activeChannels = faust_markers,
    projectPath = faust_path,
    depthScoreThreshold = 0.005,
    selectionQuantile = 1,
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
col_fun_prop <- circlize::colorRamp2(c(1, 1000, 130000), c("blue", "white", "red"))

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
                row_names_max_width = unit(20, "cm"),
                heatmap_legend_param = list(col_fun = col_fun_prop, title = "Count", 
                                            at = c(1, 100, 1000, 10000, 130000), break_dist = 1))

draw(hmap, heatmap_legend_side="right", annotation_legend_side="right")

# Save heatmap
cairo_pdf(file=here::here("out/FAUST_Treg_V3/FAUST_Treg_count_heatmap.pdf"), width=15, height=10,
          onefile = TRUE, bg = "transparent", family = "Arial") 
print(draw(hmap, gap = unit(0.1, "in"), heatmap_legend_side="right", annotation_legend_side="right"))
dev.off()

# Determine how many subpopulations are CD4+ and/or CD8a+
subpops <- colnames(countMatrix)
sum(str_count(subpops, "CD4\\+.*CD8a\\-"))
sum(str_count(subpops, "CD4\\-.*CD8a\\+"))
sum(str_count(subpops, "CD4\\+.*CD8a\\+"))
sum(str_count(subpops, "CD4\\-.*CD8a\\-"))

# Find which subpopulations are Tregs (CD4+FOXP3+CD25+)
str_count(subpops, "CD4\\+.*CD8a\\-.*FOXP3\\+.*CD25\\+")
which(grepl("CD4\\+.*CD8a\\-.*FOXP3\\+.*CD25\\+", subpops))
subpops[39]
subpops[45]
subpops[49]

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
