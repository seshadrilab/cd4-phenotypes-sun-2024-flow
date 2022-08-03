library(tidyverse)
library(flowWorkspace)
library(here)
library(ggplot2)
library(ggbeeswarm)
library(ggh4x)
source(here::here("scripts/Helper_Functions.R")) # for make_mag_plots()

# Plot background-corrected frequencies of phenotypic and functional markers among CD4 COMPASS-selected T cells

if(!dir.exists(here::here("out/Treg_CD4_COMPASS_Subsets_Analysis"))) {
  cat(sprintf("Creating folder %s\n", here::here("out/Treg_CD4_COMPASS_Subsets_Analysis")))
  dir.create(here::here("out/Treg_CD4_COMPASS_Subsets_Analysis"), recursive = T)
}

# Load GatingSet
gs <- load_gs(here::here("out/GatingSets/RSTR_Treg_GatingSet_with_COMPASS_Subsets/"))

# Arial font setup. Downloaded afms from https://github.com/microsoft/microsoft-r-open/tree/ec3fd89e5fb5794bd8149905c134ad801bb61800
Arial <- Type1Font(family = "Arial",
                   metrics = c(here::here("data/Arial_afm/ArialMT.afm"),
                               here::here("data/Arial_afm/ArialMT-Bold.afm"),
                               here::here("data/Arial_afm/ArialMT-Italic.afm"),
                               here::here("data/Arial_afm/ArialMT-BoldItalic.afm")))
windowsFonts("Arial" = windowsFont("Arial"))
pdfFonts(Arial = Arial)

# Extract frequencies
# cd4_compass_treg_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+CD25+"
# cd4_compass_foxp3_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+"
# cd4_compass_cd25_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD25+"
# cd4_compass_cd39cd73_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+CD73+"
# cd4_compass_cd39_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+"
# cd4_compass_cd73_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD73+"
# cd4_compass_ccr7_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CCR7+"
# cd4_compass_il10_parent <- "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/IL10+"
# 
# compass_parent_paths <- c(cd4_compass_treg_parent, cd4_compass_foxp3_parent, cd4_compass_cd25_parent, cd4_compass_cd39cd73_parent,
#                           cd4_compass_cd39_parent, cd4_compass_cd73_parent, cd4_compass_ccr7_parent, cd4_compass_il10_parent)

nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+CD25+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD25+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/IL10+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+CD73+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD73+",
           "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CCR7+")
nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/CD4\\+\\/CD4_COMPASS_Subsets\\/", "")

count_list <- list()

for(pop in nodes) {
  current_counts <- pData(gs) %>%
    tibble::rownames_to_column("rowname") %>%
    left_join(gs_pop_get_count_fast(gs, subpopulations = pop), by = c("rowname" = "name")) %>%
    dplyr::rename(Subpop = Count) %>%
    dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Status, Subpop, ParentCount) %>%
    mutate(Freq = (Subpop/ParentCount)*100)
  
  current_counts$Status <- factor(current_counts$Status, levels = c("Pneg", "TST+"))
  count_list[[pop]] <- current_counts
}

names(count_list) <- nodes_short

# treg_counts <- pData(gs) %>%
#   tibble::rownames_to_column("rowname") %>%
#   left_join(gs_pop_get_count_fast(gs, subpopulations = cd4_compass_treg_parent), by = c("rowname" = "name")) %>%
#   dplyr::rename(Subpop = Count) %>%
#   dplyr::select(rowname, "SAMPLE ID", "EXPERIMENT NAME", Stim, Status, Subpop, ParentCount) %>%
#   mutate(Freq = (Subpop/ParentCount)*100)
# treg_counts$Status <- factor(treg_counts$Status, levels = c("Pneg", "TST+"))


# # Define function
# # Argument "pop" is the node of interest
# plot_pop <- function(pop, freq, bg_freq, stim) {     
#   parent <- sub("(.*)\\/.*", "\\1", pop)
#   stim_dat <- freq %>%
#     mutate(stim_prop = !!as.name(pop) / ParentCount) %>%
#     select(`SAMPLE ID`, Status, stim_prop) %>%
#     drop_na()
#   bg_dat <- bg_freq %>%
#     mutate(bg_prop = !!as.name(pop) / ParentCount) %>%
#     select(`SAMPLE ID`, Status, bg_prop) %>%
#     drop_na()
#   bg_corr_dat <- stim_dat %>%
#     left_join(bg_dat, c("SAMPLE ID", "Status")) %>%
#     mutate(prop = stim_prop - bg_prop)
#   p <- wilcox.test(prop ~ Status, data = bg_corr_dat, paired = FALSE)$p.value
#   p.unadj.text <- if_else(p < 0.001, "p<0.001", paste0("p=", formatC(round(p, 3), format='f', digits=3)))
#   
#   current_plot <- ggplot(bg_corr_dat, aes(Status, prop)) +
#     geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total")) +
#     theme_bw(base_size = 22) +
#     geom_quasirandom(size=3, shape = 16, width = 0.3, aes(color=!!as.symbol("Status"))) +
#     labs(title = stim,
#          subtitle = pop,
#          y = "% CD4 T cells") +
#     theme(text = element_text(family="Arial"),
#           axis.title.x = element_blank(),
#           axis.title.y = element_text(size = 21),
#           axis.text.y = element_text(color="black", size=17),
#           axis.text.x = element_text(color="black", size=20),
#           plot.title = element_text(hjust = 0.5, size=21),
#           plot.subtitle = element_text(hjust = 0.5, size=13),
#           panel.grid.major.x = element_blank(),
#           legend.position = "none",
#           plot.margin = margin(0.3, 0.2, 0.1, 0.2, "cm")) +
#     scale_color_manual(values = c("Pneg" = "#984EA3", "TST+" = "#4DAF4A")) +
#     scale_y_continuous(labels = function(x) paste0(x*100)) +
#     force_panelsizes(rows = unit(3.5, "in"),
#                      cols = unit(3, "in"))
#   
#   plot_ylims <- ggplot_build(current_plot)$layout$panel_params[[1]]$y.range
#     
#   current_plot <- current_plot +
#     annotate("text", x = 1.5, y = plot_ylims[2] + 0.01*diff(plot_ylims),
#              label = p.unadj.text, size=5.5) +
#     coord_cartesian(ylim = c(plot_ylims[[1]], plot_ylims[[2]] + 0.09*diff(plot_ylims)))
# }

# Get nodes
# treg_node <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+",
#                "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/FOXP3+CD25+")
# treg_node_short <- str_replace(treg_node, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/", "")

# nodes <- c("/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+CD25+",
#            "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/FOXP3+",
#            "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD25+",
#            "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/IL10+",
#            "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+CD73+",
#            "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD39+",
#            "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CD73+",
#            "/Time/Cells/CD3+CD14-CD19-/Singlets/Live/CD3+ Lymphocytes/CD4+/CD4_COMPASS_Subsets/CCR7+")
# nodes_short <- str_replace(nodes, "\\/Time\\/Cells\\/CD3\\+\\CD14\\-\\CD19\\-\\/Singlets\\/Live\\/CD3\\+\\ Lymphocytes\\/CD4\\+\\/CD4_COMPASS_Subsets\\/", "")

# # Get Treg counts
# dmso_treg_freq <- subset(gs, Stim == "DMSO") %>%
#   gs_pop_get_count_with_meta(subpopulations = treg_node) %>%
#   pivot_wider(names_from = Population, values_from = Count) %>%
#   rename_at(vars(all_of(treg_node)), ~ treg_node_short) 
# 
# pp1_treg_freq <- subset(gs, Stim == "PP1") %>%
#   gs_pop_get_count_with_meta(subpopulations = treg_node) %>%
#   pivot_wider(names_from = Population, values_from = Count) %>%
#   rename_at(vars(all_of(treg_node)), ~ treg_node_short) 
# 
# tbwcl_treg_freq <- subset(gs, Stim == "TB WCL") %>%
#   gs_pop_get_count_with_meta(subpopulations = treg_node) %>%
#   pivot_wider(names_from = Population, values_from = Count) %>%
#   rename_at(vars(all_of(treg_node)), ~ treg_node_short) 

# Get phenotypic and functional marker counts
# dmso_freq <- subset(gs, Stim == "DMSO") %>%
#   gs_pop_get_count_with_meta(subpopulations = nodes) %>%
#   pivot_wider(names_from = Population, values_from = Count) %>%
#   rename_at(vars(all_of(nodes)), ~ nodes_short) 
# 
# pp1_freq <- subset(gs, Stim == "PP1") %>%
#   gs_pop_get_count_with_meta(subpopulations = nodes) %>%
#   pivot_wider(names_from = Population, values_from = Count) %>%
#   rename_at(vars(all_of(nodes)), ~ nodes_short) 
# 
# tbwcl_freq <- subset(gs, Stim == "TB WCL") %>%
#   gs_pop_get_count_with_meta(subpopulations = nodes) %>%
#   pivot_wider(names_from = Population, values_from = Count) %>%
#   rename_at(vars(all_of(nodes)), ~ nodes_short) 

# Plot frequencies and perform Wilcoxon Rank Sum test between status groups
# stim <- "PP1"
# png(file=here::here(sprintf("out/Treg_CD4_COMPASS_Subsets_Analysis/%s_%s_vs_Status.png", stim, 
#                             sub("\\/", "_", treg_node_short[2]))), width=300, height=400, units = "px")
# print(plot_pop(pop = treg_node_short[2], freq = pp1_treg_freq, stim = "PP1"))
# dev.off()
# 
# stim <- "TB WCL"
# png(file=here::here(sprintf("out/Treg_CD4_COMPASS_Subsets_Analysis/%s_%s_vs_Status.png", stim, 
#                             sub("\\/", "_", treg_node_short[2]))), width=300, height=400, units = "px")
# print(plot_pop(pop = treg_node_short[2], freq = tbwcl_treg_freq, stim = "TB WCL"))
# dev.off()

# Set color scheme and stims
fill_colors <- c("Pneg" = "#984EA3", "TST+" = "#4DAF4A")
stims <- c("DMSO", "PP1", "TB WCL")

# Get plots
foxp3_cd25_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[1]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[1],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 5))
                              })
names(foxp3_cd25_plots) <- stims

foxp3_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[2]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[2],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 50))
                                })
names(foxp3_plots) <- stims

cd25_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[3]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[3],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 70))
                                })
names(cd25_plots) <- stims

il10_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[4]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[4],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 0.3))
                                })
names(il10_plots) <- stims

cd39_cd73_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[5]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[5],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 30))
                                })
names(cd39_cd73_plots) <- stims

cd39_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[6]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[6],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 30))
                                })
names(cd39_plots) <- stims

cd73_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[7]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[7],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(0, 35))
                                })
names(cd73_plots) <- stims

ccr7_plots <- purrr::pmap(.l = list(stims),
                                .f = function(n) {
                                  make_mag_plots(count_list[[nodes_short[8]]], current_stim = n, num_comparisons = length(stims), groups_to_compare = c("Pneg", "TST+"),
                                                 paired = FALSE, adjust_p = FALSE, fill_colors = fill_colors, group_by_colname = "Status", subtitle = nodes_short[8],
                                                 y_axis_text = "% COMPASS-selected CD4 T Cells", y_axis_size = 15,   ylim = c(70, 100))
                                })
names(ccr7_plots) <- stims

# Save plots
cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_FOXP3_CD25_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(foxp3_cd25_plots[["DMSO"]])
print(foxp3_cd25_plots[["PP1"]])
print(foxp3_cd25_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_FOXP3_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(foxp3_plots[["DMSO"]])
print(foxp3_plots[["PP1"]])
print(foxp3_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD25_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd25_plots[["DMSO"]])
print(cd25_plots[["PP1"]])
print(cd25_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_IL10_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(il10_plots[["DMSO"]])
print(il10_plots[["PP1"]])
print(il10_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD39_CD73_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd39_cd73_plots[["DMSO"]])
print(cd39_cd73_plots[["PP1"]])
print(cd39_cd73_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD39_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd39_plots[["DMSO"]])
print(cd39_plots[["PP1"]])
print(cd39_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CD73_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(cd73_plots[["DMSO"]])
print(cd73_plots[["PP1"]])
print(cd73_plots[["TB WCL"]])
dev.off()

cairo_pdf(file=here::here("out/Treg_CD4_COMPASS_Subsets_Analysis/Treg_CD4_COMPASS_CCR7_plots.pdf"), width=4, height=5,
          onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica.
print(ccr7_plots[["DMSO"]])
print(ccr7_plots[["PP1"]])
print(ccr7_plots[["TB WCL"]])
dev.off()

# Plot and save figures
# stim <- "PP1"
# cairo_pdf(file=here::here(sprintf("out/Treg_CD4_COMPASS_Subsets_Analysis/%s_Freq_vs_Status.pdf", stim)), 
#           width=4, height=5, onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica
# for(pop in nodes_short[2:length(nodes_short)]) {
#   print(plot_pop(pop, freq = pp1_freq, bg_freq = dmso_freq, stim = stim))
#   }
# dev.off()
# 
# stim <- "TB WCL"
# cairo_pdf(file=here::here(sprintf("out/Treg_CD4_COMPASS_Subsets_Analysis/%s_Freq_vs_Status.pdf", stim)), 
#           width=4, height=5, onefile = TRUE, bg = "transparent", family = "Arial") # default unit is inches, default font is Helvetica
# for(pop in nodes_short[2:length(nodes_short)]) {
#   print(plot_pop(pop, freq = tbwcl_freq, bg_freq = dmso_freq, stim = stim))
# }
# dev.off()
