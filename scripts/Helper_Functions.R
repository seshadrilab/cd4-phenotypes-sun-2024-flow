# t-SNE Helper functions

# Distribute totalEvents as evenly as possible among a vector of different sub-group sizes
# Recursive function
#
# Example:
# cd1c_panel_counts %>%
#   dplyr::select(NHP, TimepointTissue, Batch, rCD1c_GMM_clean) %>% 
#   group_by(TimepointTissue, Batch) %>% 
#   nest() %>% 
#   ungroup() %>% 
#   mutate(nsamp = map2(376, data, function(totalEvents, df) {
#     distributeEvents(totalEvents, df$rCD1c_GMM_clean)
#   })) %>% 
#   unnest(cols = c(data, nsamp))
distributeEvents <- function(totalEvents, subGroupSizes) {
  #print(totalEvents)
  #print(subGroupSizes)
  if(totalEvents == 0) {
    rep(0, length(subGroupSizes))
  } else if(sum(subGroupSizes) < totalEvents) {
    stop("Not enough events to achieve requested sample size")
  } else if(any(subGroupSizes == 0)) {
    # This branch exists in case the first call to this function includes some subGroups of size 0. 
    output <- rep(0, length(subGroupSizes))
    nonZeroIndices <- which(subGroupSizes > 0)
    output[nonZeroIndices] <- distributeEvents(totalEvents, subGroupSizes[nonZeroIndices])
  } else {
    numSubGroups <- length(subGroupSizes)
    total_div_numSubGroups <- totalEvents %/% numSubGroups # ideally, how many events would be sampled from each sub-group?
    remainder <- totalEvents %% numSubGroups
    minGroupSize <- min(subGroupSizes)
    
    if(total_div_numSubGroups == 0) {
      # totalEvents is non-zero but total_div_numSubGroups is zero, meaning there is a non-zero remainder
      stopifnot(remainder != 0)
      # If the function got to this point, we know that all subGroups are non-zero (see previous branch)
      stopifnot(all(subGroupSizes > 0))
      
      # Ok, now time to distribute the remainder among subGroupSizes events
      # We know that remainder is smaller than the number of subGroups
      stopifnot(remainder < numSubGroups)
      # So now we have to evenly and randomly distribute the remainder among the subGroups
      # This is the only random portion of this function
      output <- rep(0, numSubGroups)
      output[sample.int(numSubGroups, remainder)] <- 1
      # This is a base case, so we return output. No more recursive calls.
      output
    } else {
      howMuchToAssign <- min(total_div_numSubGroups, minGroupSize)
      output <- rep(howMuchToAssign, numSubGroups)
      
      #print(output)
      
      remaining_subGroupSizes <- subGroupSizes - howMuchToAssign
      nonZeroIndices <- which(remaining_subGroupSizes > 0)
      if(length(nonZeroIndices)) {
        toAdd <- rep(0, numSubGroups)
        toAdd[nonZeroIndices] <- distributeEvents(totalEvents - sum(output), remaining_subGroupSizes[nonZeroIndices])
        output + toAdd
      } else {
        output
      }
    }
  }
}

############################################################################################################################

# A function to sample n events from a GatingHierarchy containing a single sample
sampleGatingHierarchy <- function(gh, parentGate, n = NULL, otherGates = NULL) {
  library(openCyto)
  library(CytoML) # 1.12.0
  library(flowCore) # required for description()
  library(flowWorkspace) # required for gh_pop_get_data()
  
  stopifnot(length(gh) == 1)
  allMarkerNames <- pData(parameters(gh_pop_get_data(gh)))[,c(1,2)] # First column is flow channel, second is marker name
  if (any(is.na(allMarkerNames[,2])) | length(unique(allMarkerNames[,2])) < length(allMarkerNames[,2]))
    stop ("all marker names (even FSC-A and Time) must be assigned and be unique")
  # May want to loosen above requirement, using channel names where marker names are unavailable.
  
  if(!is.null(n)) {
    # First take length n sample from all events in parentGate
    availableEvents <- gh_pop_get_count(gh, parentGate)
    nSampled <- sample.int(availableEvents, size = n)
  }
  
  # Then extract boolean gate data and mfi data, and cbind them.
  # For now, in the interest of saving memory, the data is subset to the sampled events as soon as possible.
  parentGateIndices <- gh_pop_get_indices(gh, parentGate) # relative to all data in gh, i.e. root node. a TRUE/FALSE vector. Used for subsetting boolean data
  gates2Extract <- unique(c(parentGate, if(is.null(otherGates)) {gh_get_pop_paths(gh)} else {otherGates}))
  if(!is.null(n)) {
    perCellGateMembership <- data.frame(lapply(gates2Extract, function(currentGate) {
      as.integer(gh_pop_get_indices_mat(gh, currentGate)[parentGateIndices,][nSampled]) }))
  } else {
    perCellGateMembership <- data.frame(lapply(gates2Extract, function(currentGate) {
      as.integer(gh_pop_get_indices_mat(gh, currentGate)[parentGateIndices,]) }))
  }
  colnames(perCellGateMembership) <- gates2Extract
  
  if(!is.null(n)) {
    perCellMFIData <- exprs(gh_pop_get_data(gh, parentGate))[nSampled,]
  } else {
    perCellMFIData <- exprs(gh_pop_get_data(gh, parentGate))
  }
  colnames(perCellMFIData) <- allMarkerNames[match(colnames(perCellMFIData), allMarkerNames[,1]), 2]
  
  # Combine Gate membership data with MFI expression data
  stopifnot(nrow(perCellGateMembership) == nrow(perCellMFIData))
  gateAndMFIData <- cbind(perCellGateMembership, perCellMFIData)
  
  # Add metadata and return 
  cbind(pData(gh), gateAndMFIData, row.names = NULL)
}

############################################################################################################################
############################################################################################################################

# COMPASS Helper functions

# TODO: load libraries. COMPASS, grid, ?
runCompassOnce <- function(gs,
                           seed=NULL,
                           outDir,
                           parentNode,
                           nodeMarkerMap,
                           uniqueIdentifier,
                           treatmentCol="trt",
                           currentTreatment="Treatment",
                           currentControl="Control",
                           stratifyBy=NULL,
                           iter=40000,
                           eventCountFilterThreshold=0,
                           textForRunOutputId=NULL) {
  require(COMPASS)
  require(grid)
  
  currentRunTextForConsole <- paste0(parentNode, " ", currentTreatment)
  currentRunTextForConsole <- if(is.null(textForRunOutputId)) {
    currentRunTextForConsole
  } else {
    sprintf("%s (%s)", currentRunTextForConsole, textForRunOutputId)
  }
  message(sprintf("Running COMPASS for %s", currentRunTextForConsole))
  
  # Set the seed
  if (!is.null(seed)) {
    rngKind <- "L'Ecuyer-CMRG" # This seems to be the recommended kind of RNG for reproducible parallel processing
    RNGkind(kind = rngKind)
    message(sprintf("Setting %s seed to %s", rngKind, seed))
    set.seed(seed)
  }
  
  # Create a COMPASSContainer from the GatingSet or GatingSetList.
  CC <- COMPASS::COMPASSContainerFromGatingSet(gs, node=parentNode, individual_id=uniqueIdentifier,
                                               mp=nodeMarkerMap, countFilterThreshold=eventCountFilterThreshold)
  
  # Run COMPASS
  # The treatment and control arguments for COMPASS::COMPASS() accept expressions, which makes it difficult to use COMPASS::COMPASS() programmatically.
  # This work-around uses bquote() to build up an expression consisting of a call to COMPASS::COMPASS(). eval() then evaluates this expression.
  # bquote() returns an expression consisting of the code you provided, but it first evaluates everything in .() 
  # It is not an ideal solution, but it works for now.
  fit <- eval(bquote(COMPASS::COMPASS( CC,
                                       treatment= .(as.name(treatmentCol)) == .(eval(currentTreatment)), 
                                       control= .(as.name(treatmentCol)) == .(eval(currentControl)),
                                       iterations=iter
  )))
  
  message("COMPASS complete, now saving output")
  
  if(is.null(textForRunOutputId)) {
    textForRunOutputId <- paste0(parentNode, "_", currentTreatment)
  }
  
  # Save the COMPASS run output as an RDS file for safekeeping
  saveRDS(fit, file.path(outDir, sprintf("COMPASSResult_%s.rds", textForRunOutputId)))
  
  # Save the Functionality and Polyfunctionality Scores to a tsv
  # TODO: rewrite this using stack(COMPASS::FunctionalityScore(fit)) instead of data.frame()
  FS <- COMPASS::FunctionalityScore(fit)
  PFS <- COMPASS::PolyfunctionalityScore(fit)
  FS_df <- data.frame(tmp = names(FS), FS = FS)
  colnames(FS_df) <- c(uniqueIdentifier, "FS")
  PFS_df <- data.frame(tmp = names(PFS), PFS = PFS)
  colnames(PFS_df) <- c(uniqueIdentifier, "PFS")
  FS_PFS_df <- merge(FS_df, PFS_df, by = uniqueIdentifier)
  write.table(FS_PFS_df,
              file = file.path(outDir, sprintf("FS_PFS_%s.tsv", textForRunOutputId)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Plot a heatmap of the mean probability of response
  plotTitleSuffix <- paste(c(",\n", run, ", ", parentNode, " Cells"), collapse="")
  cytokine_annotation_colors <- c("black", "black", "black", "black", "black", "black", "black")
  
  png(filename=file.path(outDir, sprintf("HeatmapMeanProbResponse_%s.png", textForRunOutputId)),
      width=800, height=650)
  # COMPASS::plot.COMPASSResult implements the generic S3 method graphics::plot(), which uses the object's class (in this case a COMPASSResult) to dispatch to the COMPASS plotting function
  # The COMPASS package doesn't allow me to call COMPASS::plot.COMPASSResult directly
  try(grid::grid.draw(print(graphics::plot(fit, stratifyBy, show_rownames = TRUE,
                                           main = sprintf("Heatmap of Mean Probability of Response %s", textForRunOutputId),
                                           fontsize=14, fontsize_row=13, fontsize_col=11,
                                           cytokine_annotation_colors=cytokine_annotation_colors))))
  dev.off()
  
  message(sprintf("Done with run %s", currentRunTextForConsole))
}

############################################################################################################################

# Like plot.COMPASSResult, but using the ComplexHeatmap package
# The argument "subsets_keep" specifies which samples to keep based on Run_SAMPLE_ID
plot.COMPASSResult.ComplexHeatmap <- function(cr,
                                              subset_keep = NULL,
                                              row_annotation = NULL,
                                              cytokine_order_for_annotation = NULL,
                                              cytokine_row_name_text = NULL,
                                              dichotomize_by_cytokine = NULL,
                                              dichotomize_by_cytokine_color = NULL,
                                              row_annotation_colors = NULL,
                                              staircase_cytokine_annotation = TRUE,
                                              row_gap = unit(0, "in")) {
  library(tidyverse)
  library(ComplexHeatmap)
  library(COMPASS)
  library(RColorBrewer)
  
  if(!is.null(row_annotation) & is.null(row_annotation_colors)) { stop("row_annotation_colors must not be NULL if row_annotation is not NULL")}
  
  mean_gamma <- cr$fit$mean_gamma
  cats <- as.data.frame(cr$fit$categories[, -ncol(cr$fit$categories),drop=FALSE]) # drop the "Counts" column
  numMarkers <- ncol(cats)
  rownames(cats) <- colnames(mean_gamma)
  meta <- cr$data$meta
  
  # Subset individuals/samples to keep based on rowname (Run_SAMPLE_ID)
  if(!is.null(subset_keep)) {
    mean_gamma <- mean_gamma[grep(subset_keep, rownames(mean_gamma)), ]
    meta <- meta[grep(subset_keep, rownames(meta)), ]
  }
  
  # Filter the subsets to those where the average mean_gamma is greater than the threshold (default in heatmap and this function is 0.01)
  compassSubsetsFiltered <- names(which(apply(mean_gamma, 2, function(x) { mean(x, na.rm = TRUE) }) > 0.01))
  # And remove the subset with 0 positive markers
  compassSubsetsFiltered <- compassSubsetsFiltered[lengths(regmatches(compassSubsetsFiltered, gregexpr("!", compassSubsetsFiltered))) != numMarkers]
  
  # Subset the cats rows to compassSubsetsFiltered, and put the columns in the order of cytokine_order_for_annotation
  cats <- cats[compassSubsetsFiltered,]
  if(!is.null(cytokine_order_for_annotation)) {
    cats <- cats[, cytokine_order_for_annotation]
  }
  
  # Before plotting, put the categories data frame rows in the desired order (columns were already re-ordered above)
  # This is essentially the same code I added to the pheatmap function
  cats <- cats[rev(do.call(order, cats)),,drop=FALSE]
  # And order the cats df rows by degrees (number of cytokines in subset)
  ckr<-apply(cats,1,function(x)sum(as.numeric(as.character(x))))
  cats = cats[order(ckr),]
  if(!is.null(dichotomize_by_cytokine)) {
    # And then dichotomize the cats df rows so that all subsets containing the cytokine in dichotomize_by_cytokine (e.g. "IFNg") appear last
    cats <- cats[order(cats[,dichotomize_by_cytokine]),]
  }
  
  # Now re-order the columns of mean_gamma to match the rows of cats
  mean_gamma <- mean_gamma[,rownames(cats)]
  
  # If dichotomizing and coloring by a cytokine, specify which cells in the cats matrix should be colored differently
  if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
    current_cats_rownames <- rownames(cats)
    cats <- cats %>% 
      mutate_all(~ ifelse(. == 1 & !!as.symbol(dichotomize_by_cytokine) > 0, 2, .))
    rownames(cats) <- current_cats_rownames
  }
  
  # Now re-order the rows of mean_gamma by FunctionalityScore
  FS_order <- order(FunctionalityScore(mean_gamma, n = numMarkers), decreasing=T)
  mean_gamma <- mean_gamma[FS_order,]
  # And make sure to update the metadata rows
  meta <- meta[FS_order,]
  stopifnot(all.equal(rownames(mean_gamma), rownames(meta)))
  
  # Now additionally order the rows of meta and then mean_gamma based on row_annotation
  suppressWarnings(if(!is.null(row_annotation)) {
    if(length(row_annotation) == 1) {
      meta <- meta %>% arrange(!! rlang::sym(row_annotation))
    } else if(length(row_annotation) > 1) {
      meta <- meta %>% arrange(!!! rlang::syms(row_annotation))
    }
    mean_gamma <- mean_gamma[match(rownames(meta), rownames(mean_gamma)),]
  })
  
  ht_opt$simple_anno_size = unit(2.5, "mm")
  heatmap_main <- Heatmap(mean_gamma,
                          cluster_rows = FALSE, 
                          show_row_dend = FALSE, 
                          cluster_columns = FALSE, cluster_column_slices = FALSE,
                          row_split = suppressWarnings(if(!is.null(row_annotation)) { meta %>% dplyr::select(all_of(row_annotation)) } else { NULL }),
                          border = "black",
                          show_column_names = FALSE,
                          show_row_names = FALSE, row_title = NULL,
                          right_annotation = suppressWarnings(if(!is.null(row_annotation)) { 
                            HeatmapAnnotation(df = meta %>% dplyr::select(all_of(row_annotation)),
                                              col = row_annotation_colors,
                                              which = "row", border = F,
                                              show_annotation_name = F,
                                              annotation_legend_param = list(border=F,
                                                                             title_gp = gpar(fontsize = 13, fontface = "bold"),
                                                                             labels_gp = gpar(fontsize = 12)))
                          } else {
                            NULL
                          }),
                          heatmap_legend_param = list(border=F, title_gp = gpar(fontsize = 13, fontface = "bold"),
                                                      labels_gp = gpar(fontsize = 12), legend_direction = "horizontal", title_position = "topcenter"),
                          name="Response", use_raster = F,
                          col = colorRampPalette(brewer.pal(9,"Purples"))(20),
                          height = unit(5, "in"),
                          row_gap = row_gap)
  # draw(heatmap_main)
  if(!is.null(cytokine_order_for_annotation)) {
    colnames(cats) <- cytokine_row_name_text
  }
  heatmap_cats <- Heatmap(cats %>% dplyr::select(rev(everything())) %>% t(),
                          cluster_rows = FALSE, 
                          show_row_dend = FALSE, 
                          cluster_columns = FALSE,
                          border = T,
                          show_column_names = FALSE,
                          show_row_names = TRUE,
                          row_names_gp = gpar(fontsize=14),
                          row_names_side = "left",
                          row_title = NULL,
                          col = if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
                            c("white", "black", dichotomize_by_cytokine_color) } else { c("white", "black") },
                          show_heatmap_legend = FALSE,
                          rect_gp = gpar(col = "white", lwd = 1),
                          height = unit(1.2, "in"))
  # draw(heatmap_cats)
  
  # draw(heatmap_main %v% heatmap_cats, gap = unit(0.1, "in"), merge_legends = T) # heatmap_legend_side = "left"
  # ht_opt(RESET = TRUE)
  
  heatmap_main %v% heatmap_cats
}

############################################################################################################################

# Boxplot of background-corrected subset percents out of the parent population
# cr is the COMPASSResult object
make_boxplot_for_COMPASS_run <- function(cr, subset_keep = NULL, run_name, output_folder=NA, current_ylim=NULL, add_legend=FALSE, legend_position=c(0.02, 1),
                                         paired=FALSE, important_label=NULL, include_label=FALSE, color_important=FALSE, save_test_results=TRUE, p_text_size=5, include_0_line=FALSE, zeroed_BgCorr_stats = FALSE, zeroed_BgCorr_plot = FALSE, 
                                         plot_width=7, plot_height=6, dichotomize_by_cytokine=NULL, dichotomize_by_cytokine_color=NULL, group_by_colname="Group",
                                         group_by_order=c("Hospitalized", "Non-hospitalized"),
                                         group_by_colors=c("Hospitalized" = "#757bbcb2", "Non-hospitalized" = "#b0d2c8bf"),
                                         group_by_labels = c("Conv Hosp", "Conv Non-Hosp"),
                                         parentSubset="CD4+", cytokine_order_for_annotation=NULL, cytokine_row_name_text=NULL,
                                         cats_heatmap_left_padding = 3, main_title = NULL, onlyShowPBelowAlpha = TRUE, adj_pval = FALSE) {
  # zeroed_BgCorr_stats If TRUE, the returned magnitude values and statistics are calculated using zeroed values
  # zeroed_BgCorr_plot If TRUE, the points on the plot are zeroed (regardless of zeroed_BgCorr_stats)
  library(ComplexHeatmap)
  library(grid)
  library(gridExtra)
  library(gtable)
  
  mean_gamma <- cr$fit$mean_gamma
  cats <- as.data.frame(cr$fit$categories[, -ncol(cr$fit$categories),drop=FALSE]) # drop the "Counts" column
  numMarkers <- ncol(cats)
  rownames(cats) <- colnames(mean_gamma)
  
  # Subset individuals/samples to keep based on "SAMPLE ID"
  if(!is.null(subset_keep)) {
    mean_gamma <- mean_gamma[grep(subset_keep, rownames(mean_gamma)), ]
  }
  
  # Filter the subsets to those where the average mean_gamma is greater than the threshold (default in heatmap and this function is 0.01)
  compassSubsetsFiltered <- names(which(apply(mean_gamma, 2, function(x) { mean(x, na.rm = TRUE) }) > 0.01))
  # Note that we don't need mean_gamma anymore for this task
  # And remove the subset with 0 positive markers
  compassSubsetsFiltered <- compassSubsetsFiltered[lengths(regmatches(compassSubsetsFiltered, gregexpr("!", compassSubsetsFiltered))) != numMarkers]
  
  # Subset the cats rows to compassSubsetsFiltered, and put the columns in the order of cytokine_order_for_annotation
  cats <- cats[compassSubsetsFiltered,]
  if(!is.null(cytokine_order_for_annotation)) {
    cats <- cats[, cytokine_order_for_annotation]
  }
  
  stim_counts <- as.data.frame(cr$data$n_s) %>%
    mutate(Individual = rownames(cr$data$n_s)) %>%
    dplyr::select(c("Individual", all_of(compassSubsetsFiltered))) %>% 
    dplyr::left_join(cr$data$counts_s %>%
                       stack() %>%
                       rename("ParentCount" = "values", "Individual" = "ind"),
                     by = "Individual") %>%
    dplyr::left_join(cr$data$meta %>% 
                       dplyr::select(!!as.symbol(cr$data$individual_id), !!as.symbol(group_by_colname)),
                     by = c("Individual"=cr$data$individual_id)) %>% 
    mutate(Stim = "Dummy_Stim_Name")
  
  bg_counts <- as.data.frame(cr$data$n_u) %>%
    mutate(Individual = rownames(cr$data$n_u)) %>%
    dplyr::select(c("Individual", all_of(compassSubsetsFiltered))) %>% 
    dplyr::left_join(cr$data$counts_u %>%
                       stack() %>%
                       rename("ParentCount" = "values", "Individual" = "ind"),
                     by = "Individual") %>%
    dplyr::left_join(cr$data$meta %>% 
                       dplyr::select(!!as.symbol(cr$data$individual_id), !!as.symbol(group_by_colname)),
                     by = c("Individual"=cr$data$individual_id)) %>% 
    mutate(Stim = "DMSO")
  
  dat_bgCorr_long <- bind_rows(bg_counts, stim_counts) %>% 
    mutate_at(.vars = compassSubsetsFiltered, `/`, quote(ParentCount)) %>%  # convert counts to proportions
    dplyr::select(-ParentCount) %>% 
    tidyr::pivot_longer(cols = -c("Individual", "Stim", !!as.symbol(group_by_colname)),
                        names_to = "BooleanSubset",
                        values_to = "Proportion") %>% 
    tidyr::pivot_wider(id_cols = c("Individual", !!as.symbol(group_by_colname), "BooleanSubset"),
                       names_from = Stim,
                       values_from = Proportion) %>% 
    drop_na() %>% # Filter out rows with NA
    mutate(BgCorr = if(zeroed_BgCorr_stats) {pmax(0, Dummy_Stim_Name - DMSO)} else {Dummy_Stim_Name - DMSO},
           BgCorr_plot = if(zeroed_BgCorr_plot) {pmax(0, Dummy_Stim_Name - DMSO)} else {Dummy_Stim_Name - DMSO},
           !!as.symbol(group_by_colname) := factor(!!as.symbol(group_by_colname), levels = group_by_order)) %>% 
    dplyr::select(-c(DMSO, Dummy_Stim_Name))
  
  dat_bgCorr_wide <- dat_bgCorr_long %>% 
    tidyr::pivot_wider(id_cols = c("Individual", !!as.symbol(group_by_colname)),
                       names_from = BooleanSubset,
                       values_from = BgCorr)
  
  # Subset individuals/samples to keep based on "SAMPLE ID"
  if(!is.null(subset_keep)) {
    dat_bgCorr_wide <- dat_bgCorr_wide[grep(subset_keep, dat_bgCorr_wide$Individual), ]
    dat_bgCorr_long <- dat_bgCorr_long[grep(subset_keep, dat_bgCorr_long$Individual), ]
  }
  
  # Non-parametric statistical hypothesis tests
  if(paired) {
    # Need both time points for a single PTID to do Wilcoxon signed rank test
    # Drop any PTIDs without both time points
    dat_bgCorr_wide <- dat_bgCorr_wide %>%
      mutate(PTID = str_replace(.$Individual, c("-1|-2"), "")) %>%
      group_by(PTID) %>%
      filter(n()>1) %>% 
      ungroup() %>%
      select(-PTID)
    
    dat_bgCorr_long <- dat_bgCorr_long %>%
      mutate(PTID = str_replace(.$Individual, c("-1|-2"), "")) %>%
      group_by(PTID, BooleanSubset) %>%
      filter(n()>1) %>% 
      ungroup() %>%
      select(-PTID)
    
    tests <- lapply(compassSubsetsFiltered, function(boolSubset) {
      wilcox.test(as.formula(sprintf("`%s` ~ %s", boolSubset, group_by_colname)), data=dat_bgCorr_wide, paired=T)
    }) 
  } else {
    # Wilcoxon rank sum test works with unequal sample sizes 
    tests <- lapply(compassSubsetsFiltered, function(boolSubset) {
      wilcox.test(as.formula(sprintf("`%s` ~ %s", boolSubset, group_by_colname)), data=dat_bgCorr_wide)
    }) 
  }
  
  if(adj_pval) {
    pvals_df <- data.frame(BooleanSubset = compassSubsetsFiltered,
                           p = unlist(lapply(tests, function(x) {x$p.value}))) %>% 
      mutate(p.adj = p.adjust(p, method = "bonferroni")) %>% # Strict
      mutate(p.adj.text = if_else(p.adj < 0.001, "p<.001", paste0("p=", sub("0.", ".", round(p.adj, 3))))) 
  } else {
    pvals_df <- data.frame(BooleanSubset = compassSubsetsFiltered,
                           p = unlist(lapply(tests, function(x) {x$p.value}))) %>% 
      mutate(p.text = if_else(p < 0.05, "*", paste0("p=", sub("0.", ".", round(p, 3)))))
  }
  
  # Before plotting, put the categories data frame rows in the desired order (columns were already re-ordered above)
  # This is essentially the same code I added to the pheatmap function
  cats <- cats[rev(do.call(order, cats)),,drop=FALSE]
  # And order the cats df rows by degrees (number of cytokines in subset)
  ckr<-apply(cats,1,function(x)sum(as.numeric(as.character(x))))
  cats = cats[order(ckr),]
  if(!is.null(dichotomize_by_cytokine)) {
    # And then dichotomize the cats df rows so that all subsets containing the cytokine in dichotomize_by_cytokine (e.g. "IFNg") appear last
    cats <- cats[order(cats[,dichotomize_by_cytokine]),]
  }
  
  # Use the cats df row order to order the boolean subsets in dat_bgCorr_long and pvals_df
  dat_bgCorr_long$BooleanSubset <- factor(dat_bgCorr_long$BooleanSubset, levels = rownames(cats))
  pvals_df$BooleanSubset <- factor(pvals_df$BooleanSubset, levels = rownames(cats))
  
  # If dichotomizing and coloring by a cytokine, specify which cells in the cats matrix should be colored differently
  if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
    current_cats_rownames <- rownames(cats)
    cats <- cats %>% 
      mutate_all(~ ifelse(. == 1 & !!as.symbol(dichotomize_by_cytokine) > 0, 2, .))
    rownames(cats) <- current_cats_rownames
  }
  
  # Calculate medians of each group for each subset
  dat_bgCorr_medians <- dat_bgCorr_long %>%
    dplyr::group_by(!!as.symbol(group_by_colname), BooleanSubset) %>%
    dplyr::summarise(BgCorr = median(BgCorr),
                     BgCorr_plot = median(BgCorr_plot))
  
  if(save_test_results) {
    # Save some form of the test results to disk so it doesn't only exist in the plot as adjusted p-values
    pvals_df_for_file <- cats %>% rownames_to_column("BooleanSubset") %>% 
      dplyr::left_join(pvals_df) %>% 
      dplyr::left_join(dat_bgCorr_medians %>%
                         pivot_wider(id_cols = "BooleanSubset",
                                     names_from = !!as.symbol(group_by_colname),
                                     values_from = BgCorr,
                                     names_prefix = "med_")) 
    
    if(adj_pval) {
      pvals_df_for_file <- arrange(pvals_df_for_file, p.adj)
    } else {
      pvals_df_for_file <- arrange(pvals_df_for_file, p)
    }
    
    if(!is.na(output_folder)) {
      test_results_file_path <- file.path(output_folder,
                                          sprintf("%s_%s_%s_BooleanSubsets_BgCorrProps_MannWhitney.tsv",
                                                  run_name,
                                                  if(subset_keep == "H-") {
                                                    "N"
                                                    } else if(subset_keep == "C-") {
                                                      "C"
                                                    } else if(subset_keep == "H-1|C-1") {
                                                      "PRE"
                                                    } else {"POST"},
                                                  if(zeroed_BgCorr_stats) {"Zeroed"} else {"NotZeroed"}))
      write.table(pvals_df_for_file,
                  file = test_results_file_path,
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
  # Add a new variable to dat_bgCorr_long to indicate the important samples to label
  dat_bgCorr_long$important <- dplyr::case_when(
    grepl(paste(important_label, collapse = "|"), dat_bgCorr_long$Individual) ~ TRUE,
    TRUE ~ FALSE
  )
  
  # Draw the boxplot
  p_boxplot <- ggplot(dat_bgCorr_long, aes(x = !!as.symbol(group_by_colname), y = BgCorr_plot,
                                           #fill = !!as.symbol(group_by_colname),
                                           group = !!as.symbol(group_by_colname)))
  if(include_0_line) {
    p_boxplot <- p_boxplot + geom_hline(yintercept = 0, linetype="dashed", alpha = 0.5)
  }
  p_boxplot <- p_boxplot +
    geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total")) +
    geom_quasirandom(size=1.25, shape = 16, width = 0.2, aes(color=!!as.symbol(group_by_colname))) +
    facet_grid(. ~ BooleanSubset) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text = element_text(color="black", size=18),
          axis.title = element_text(size=25),
          text = element_text(family="Arial"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + #,
    # plot.margin = add_margin(b=0, unit="cm")) + 
    labs(y=sprintf("%% Responding %s T cells", sub("+", "", parentSubset, fixed=T)))
  
  if(include_label) {
    p_boxplot <- p_boxplot +
    geom_text(aes(label = ifelse(important, Individual, ""), group = Individual),
              position = position_dodge(width = 0.9), size = 2.5)
  }
  
  if(color_important) {
    p_boxplot <- p_boxplot +
      geom_point(data=dat_bgCorr_long[dat_bgCorr_long$important == TRUE,], 
                 color = "red", size = 0.5)
  }
  
  if(!is.null(main_title)) {
    p_boxplot <- p_boxplot +
      theme(plot.title = element_text(hjust=0.5, size=26)) +
      labs(title = main_title)
  }
  if(!is.null(group_by_colors)) {
    p_boxplot <- p_boxplot + if(is.null(group_by_labels)) {
      scale_color_manual(values=group_by_colors)
    } else {
      scale_color_manual(values=group_by_colors,
                        labels=group_by_labels)
    }
  }
  
  # Adjust ylim here manually if specified
  if(is.null(current_ylim)) {
    p_boxplot <- p_boxplot +
      scale_y_continuous(labels = function(x) paste0(x*100))
  } else {
    p_boxplot <- p_boxplot +
      scale_y_continuous(labels = function(x) paste0(x*100), limits=current_ylim)
  }
  
  if(add_legend) {
    p_boxplot <- p_boxplot +
      theme(legend.justification = c(0,1),
            legend.position = legend_position,
            legend.text=element_text(size=16),
            legend.title = element_blank())
  } else {
    p_boxplot <- p_boxplot +
      theme(legend.position = "none")
  }
  
  showSignificanceBracket <- TRUE
  p_alpha <- 0.05
  #onlyShowPBelowAlpha <- TRUE
  if(showSignificanceBracket) {
    
    get_y_pos <- function(boolSubsets) {
      sapply(boolSubsets, function(boolSubset) {
        boolSubset <- as.character(boolSubset)
        
        find_whisker_max <- function(x) {
          # From geom_boxplot: The upper whisker extends from the hinge to the largest value no further
          # than 1.5 * IQR from the hinge (where IQR is the inter-quartile range, or distance between the first and third quartiles).
          upper_whisker_limit <- quantile(x, probs = c(0.75)) + 1.5*IQR(x);
          max(x[x <= upper_whisker_limit])
        }
        whisker_max <- dat_bgCorr_wide %>%
          group_by(!!as.symbol(group_by_colname)) %>%
          summarise(whisker_max = find_whisker_max(!!as.symbol(boolSubset))) %>% 
          dplyr::pull(whisker_max) %>% 
          max()
        whisker_max + if(is.null(current_ylim)) {max(dat_bgCorr_long$BgCorr)/20} else {current_ylim[[2]]/20}
      })
      
    }
    
    annotation_df <- pvals_df %>% 
      mutate(start = group_by_order[[1]],
             end = group_by_order[[2]],
             y_pos = get_y_pos(BooleanSubset))
    
    if(onlyShowPBelowAlpha) {
      if(adj_pval) {
        annotation_df <- annotation_df %>% dplyr::filter(p.adj < p_alpha) 
      } else {
        annotation_df <- annotation_df %>% dplyr::filter(p < p_alpha)
      }
    }
    
    # If I don't use the full path for ggsignif::geom_signif, it may try to use a global environment variable GeomSignif and ignore manual = T. Odd.
    if(adj_pval) {
      p_boxplot <- p_boxplot +
        ggsignif::geom_signif(inherit.aes=F,data=annotation_df,
                              aes_string(xmin="start", xmax="end", annotations="p.adj.text", y_position="y_pos"), # , family="Arial"
                              tip_length = c(0, 0),
                              fontface = "bold",
                              textsize=p_text_size,
                              size = 1,
                              manual = TRUE)
    } else {
      p_boxplot <- p_boxplot +
        ggsignif::geom_signif(inherit.aes=F,data=annotation_df,
                              aes_string(xmin="start", xmax="end", annotations="p.text", y_position="y_pos"), # , family="Arial"
                              tip_length = c(0, 0),
                              fontface = "bold",
                              textsize=p_text_size,
                              size = 0.75,
                              manual = TRUE)
    }
  }
  
  # Now make the categories legend
  if(!is.null(cytokine_order_for_annotation)) {
    colnames(cats) <- cytokine_row_name_text
  }
  heatmap_cats <- Heatmap(cats %>% dplyr::select(rev(everything())) %>% t(),
                          cluster_rows = FALSE, 
                          show_row_dend = FALSE, 
                          cluster_columns = FALSE,
                          border = T,
                          show_column_names = FALSE,
                          show_row_names = TRUE,
                          row_names_gp = gpar(fontsize=14),
                          row_names_side = "left",
                          row_title = NULL,
                          col = if(!is.null(dichotomize_by_cytokine) & !is.null(dichotomize_by_cytokine_color)) {
                            c("white", "black", dichotomize_by_cytokine_color) } else { c("white", "black") },
                          show_heatmap_legend = FALSE,
                          rect_gp = gpar(col = "white", lwd = 1),
                          height = unit(1.2, "in"))
  
  # alignment help from https://support.bioconductor.org/p/103113/
  heatmap_cats_grob <- grid.grabExpr(draw(heatmap_cats, padding = unit(c(0,cats_heatmap_left_padding,0,2), "mm")))
  p_boxplot_grob <- ggplotGrob(p_boxplot)
  combined_plot <- arrangeGrob(p_boxplot_grob, heatmap_cats_grob, nrow = 2, heights = c(1, 0.28))
  
  to_return <- list("Plot" = combined_plot,
                    "BgCorrMagnitudes" = dat_bgCorr_wide)
  if(save_test_results) {
    to_return$Test_Results <- pvals_df_for_file
  }
  to_return
}

############################################################################################################################

# Plots of functionality scores and polyfunctionality scores for all stims for a given timepoint and infection status group
# Use Wilcoxon signed-rank test to compare across stim within a donor.
# Show the medians in each group.
# Note: all p-values unadjusted

fs_pfs_plot <- function(df, FS_or_PFS = "FS", cd4_or_cd8 = "CD4", group = "Naive PRE", fill_colors = NULL, group_by_colname = NULL) {
  d <- df %>% dplyr::filter(parent == cd4_or_cd8 & Group == group)
  
  d_wide <- d %>%  
    pivot_wider(id_cols = c(`PATIENT ID`, Group), names_from = Stim, values_from = !!as.name(FS_or_PFS)) %>%
    na.omit() # can't perform signed-rank test with missing data 
  
  # Signed rank test
  n_s1_signed_rank_result <- wilcox.test(d_wide$NCAP, d_wide$S1, paired = T)
  s1_s2_signed_rank_result <- wilcox.test(d_wide$S1, d_wide$S2, paired = T)
  n_s2_signed_rank_result <- wilcox.test(d_wide$NCAP, d_wide$S2, paired = T)
  
  stim_x_order <- 1:3
  d$Stim <- factor(d$Stim, levels = c("NCAP", "S1", "S2"))
  names(stim_x_order) <- levels(d$Stim)
  test_df <- data.frame(group1 = c("NCAP", "S1", "NCAP"),
                                group2 = c("S1", "S2", "S2"),
                                p.val = c(n_s1_signed_rank_result$p.value,
                                          s1_s2_signed_rank_result$p.value,
                                          n_s2_signed_rank_result$p.value)) %>% 
    mutate(group1.xloc = unname(stim_x_order[group1]),
           group2.xloc = unname(stim_x_order[group2]),
           p_val_text = sapply(p.val, function(p) {
             if(p < 0.001) {
               "p<0.001"
             } else {
               paste0("p=", round(p, 3))
             }
           }),
           y_pos = c(d %>% dplyr::filter(Stim %in% c("NCAP", "S1")) %>%
                       dplyr::pull(!!FS_or_PFS) %>% max() + 0.1*diff(range(d %>% dplyr::pull(!!FS_or_PFS))),
                     d %>% dplyr::filter(Stim %in% c("S1", "S2")) %>%
                       dplyr::pull(!!FS_or_PFS) %>% max() + 0.1*diff(range(d %>% dplyr::pull(!!FS_or_PFS))),
                     d %>% dplyr::filter(Stim %in% c("NCAP", "S1", "S2")) %>%
                       dplyr::pull(!!FS_or_PFS) %>% max() + 0.25*diff(range(d %>% dplyr::pull(!!FS_or_PFS)))),
           geom_signif_group = paste0(group1, group2))
  
  medians <- d %>% group_by(Group, Infection_Status, Timepoint, Stim) %>%
    summarise(med = median(!!as.name(FS_or_PFS))) %>% ungroup() %>% 
    mutate(x.min.segment = 1:3 - 0.1,
           x.end.segment = 1:3 + 0.1)
  
  ggplot(d, aes(Stim, !!as.name(FS_or_PFS))) +
    geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total")) +
    geom_quasirandom(size=3, shape = 16, width = 0.3, aes(color=!!as.symbol(group_by_colname))) +
    labs(y = if(FS_or_PFS == "FS") {sprintf("%s Functionality Score", cd4_or_cd8)} else if(FS_or_PFS == "PFS") {sprintf("%s Polyfunctionality Score", cd4_or_cd8)},
         title = group) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(color="black", size=17),
          axis.text.x = element_text(color="black", size=20),
          plot.title = element_text(hjust = 0.5, size=21),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(0.3, 0.2, 0.1, 0.2, "cm")) +
    scale_x_discrete(expand = c(0.15, 0.15),
                     labels = c("NCAP", "S1", "S2")) +
    scale_color_manual(values = fill_colors[[group]]) +
    ggsignif::geom_signif(inherit.aes=F,data=test_df,
                          aes_string(xmin="group1.xloc", xmax="group2.xloc",
                                     annotations="p_val_text", y_position="y_pos",
                                     group="geom_signif_group"), # , family="Arial"
                          tip_length = c(0, 0),
                          textsize=5,
                          size = 0.75,
                          manual = TRUE) +
    coord_cartesian(ylim = c(NA, max(test_df$y_pos) + 0.1*diff(range(d %>% dplyr::pull(!!FS_or_PFS)))))
}

############################################################################################################################

split_fs_pfs_plot <- function(df, FS_or_PFS = "FS", current_stim = "S1",
                              cd4_or_cd8 = "CD4", compare_naive = FALSE,
                              compare_conv = FALSE, compare_POST = FALSE,
                              compare_PRE = FALSE, compare_intra = FALSE,
                              group_by_colname = NULL, groups_to_compare = NULL, 
                              fill_colors = NULL, ylim = NULL) {
  d <- df %>% 
    dplyr::filter(Stim == current_stim & parent == cd4_or_cd8)
  d_pre <- d %>% 
    dplyr::filter(Timepoint == "PRE")
  d_post <- d %>%
    dplyr::filter(Timepoint == "POST")
  
  d_wide <- d %>%  
    pivot_wider(id_cols = c(`PATIENT ID`, Infection_Status), names_from = Timepoint, values_from = !!as.name(FS_or_PFS)) %>%
    na.omit() # can't perform signed-rank test with missing data
  d_wide_naive <- d_wide %>% dplyr::filter(Infection_Status == "Naive")
  d_wide_conv <- d_wide %>% dplyr::filter(Infection_Status == "Conv")
  
  # Signed rank test for Naive, pre-to-post-vac
  naive_signed_rank_result <- wilcox.test(d_wide_naive$PRE, d_wide_naive$POST, paired = T)
  # Signed rank test for Convalescent, pre-to-post-vac
  conv_signed_rank_result <- wilcox.test(d_wide_conv$PRE, d_wide_conv$POST, paired = T)
  
  # mann-whitney/wilcoxon rank-sum test for Pre-vac, across infection status group
  pre_mann_whitney_result <- wilcox.test(d_pre %>% dplyr::pull(!!FS_or_PFS) ~ as.factor(d_pre %>% dplyr::pull("Infection_Status")))
  # mann-whitney/wilcoxon rank-sum test for Post-vac, across infection status group
  post_mann_whitney_result <- wilcox.test(d_post %>% dplyr::pull(!!FS_or_PFS) ~ as.factor(d_post %>% dplyr::pull("Infection_Status")))
  
  timepoint_infection_status_x_order <- 1:4
  names(timepoint_infection_status_x_order) <- levels(d$Group)
  test_df <- data.frame(group1 = c("Naive PRE", "Conv PRE", "Naive PRE", "Naive POST"),
                                group2 = c("Naive POST", "Conv POST", "Conv PRE", "Conv POST"),
                                p.val = c(naive_signed_rank_result$p.value,
                                          conv_signed_rank_result$p.value,
                                          pre_mann_whitney_result$p.value,
                                          post_mann_whitney_result$p.value)) %>% 
    mutate(group1.xloc = unname(timepoint_infection_status_x_order[group1]),
           group2.xloc = unname(timepoint_infection_status_x_order[group2]),
           p_val_text = sapply(p.val, function(p) {
             if(p < 0.001) {
               "p<0.001"
             } else {
               paste0("p=", formatC(round(p, 3), format='f', digits=3))
             }
           }),
           # y_pos = c(d %>% dplyr::filter(Group %in% c("Naive PRE", "Naive POST")) %>%
           #             dplyr::pull(!!FS_or_PFS) %>% max() + 0.1*diff(range(d %>% dplyr::pull(!!FS_or_PFS))),
           #           d %>% dplyr::filter(Group %in% c("Conv PRE", "Conv POST")) %>%
           #             dplyr::pull(!!FS_or_PFS) %>% max() + 0.1*diff(range(d %>% dplyr::pull(!!FS_or_PFS))),
           #           d %>% dplyr::filter(Group %in% c("Naive PRE", "Naive POST", "Conv PRE")) %>%
           #             dplyr::pull(!!FS_or_PFS) %>% max() + 0.25*diff(range(d %>% dplyr::pull(!!FS_or_PFS))),
           #           d %>% dplyr::filter(Group %in% c("Naive PRE", "Naive POST", "Conv PRE", "Conv POST")) %>%
           #             dplyr::pull(!!FS_or_PFS) %>% max() + 0.4*diff(range(d %>% dplyr::pull(!!FS_or_PFS)))),
           geom_signif_group = paste0(group1, group2))
  
  medians <- d %>% group_by(Group, Infection_Status, Timepoint) %>%
    summarise(med = median(!!as.name(FS_or_PFS))) %>% ungroup() %>% 
    mutate(x.min.segment = 1:4 - 0.1,
           x.end.segment = 1:4 + 0.1)
  
  group_lab <- c("Naive PRE" = "Naive\nPRE", 
                 "Naive POST" = "Naive\nPOST",
                 "Conv PRE" = "Conv\nPRE",
                 "Conv POST" = "Conv\nPOST")
  
  if(compare_naive){
    d <- d %>%
      filter(Infection_Status %in% "Naive")
    
    medians <- medians %>%
      filter(Infection_Status %in% "Naive") %>%
      select(-c(x.min.segment, x.end.segment)) %>%
      add_column(x.min.segment = c(0.9,1.9), x.end.segment = c(1.1,2.1))
    
    group_lab <- c("Naive PRE" = "Naive\nPRE",
                   "Naive POST" = "Naive\nPOST")
    
    test_df <- test_df %>%
      filter(geom_signif_group %in% "Naive PRENaive POST") %>%
      mutate(group1.xloc = replace(group1.xloc, group1.xloc != 1, 1)) %>%
      mutate(group2.xloc = replace(group2.xloc, group2.xloc != 2, 2))
  }
  
  if(compare_conv){
    d <- d %>%
      filter(Infection_Status %in% "Conv")
    
    medians <- medians %>%
      filter(Infection_Status %in% "Conv") %>%
      select(-c(x.min.segment, x.end.segment)) %>%
      add_column(x.min.segment = c(0.9,1.9), x.end.segment = c(1.1,2.1))
    
    group_lab <- c("Conv PRE" = "Conv\nPRE",
                   "Conv POST" = "Conv\nPOST")
    
    test_df <- test_df %>%
      filter(geom_signif_group %in% "Conv PREConv POST") %>%
      mutate(group1.xloc = replace(group1.xloc, group1.xloc != 1, 1)) %>%
      mutate(group2.xloc = replace(group2.xloc, group2.xloc != 2, 2))
  }
  
  if(compare_POST){
    d <- d %>%
      filter(Timepoint %in% "POST")
    
    medians <- medians %>%
      filter(Timepoint %in% "POST") %>%
      select(-c(x.min.segment, x.end.segment)) %>%
      add_column(x.min.segment = c(0.9,1.9), x.end.segment = c(1.1,2.1))
    
    group_lab <- c("Naive POST" = "Naive\nPOST",
                   "Conv POST" = "Conv\nPOST")
    
    test_df <- test_df %>%
      filter(geom_signif_group %in% "Naive POSTConv POST") %>%
      mutate(group1.xloc = replace(group1.xloc, group1.xloc != 1, 1)) %>%
      mutate(group2.xloc = replace(group2.xloc, group2.xloc != 2, 2))
  }
  
  if(compare_PRE){
    d <- d %>%
      filter(Timepoint %in% "PRE")
    
    medians <- medians %>%
      filter(Timepoint %in% "PRE") %>%
      select(-c(x.min.segment, x.end.segment)) %>%
      add_column(x.min.segment = c(0.9,1.9), x.end.segment = c(1.1,2.1))
    
    group_lab <- c("Naive PRE" = "Naive\nPRE",
                   "Conv PRE" = "Conv\nPRE") 
    
    test_df <- test_df %>%
      filter(geom_signif_group %in% "Naive PREConv PRE") %>%
      mutate(group1.xloc = replace(group1.xloc, group1.xloc != 1, 1)) %>%
      mutate(group2.xloc = replace(group2.xloc, group2.xloc != 2, 2))
  }
  
  if(compare_intra) {
    plot <- ggplot(d, aes(Group, !!as.name(FS_or_PFS))) +
      geom_line(aes(group = `PATIENT ID`), color = fill_colors[[groups_to_compare]][1], size = 0.6, alpha = 0.5) +
      geom_segment(data = medians,
                   aes(x=x.min.segment, xend=x.end.segment, y=med, yend=med),
                   inherit.aes=FALSE, size = 1) +
      geom_point(size = 3, shape = 21, stroke = 1.5, aes(fill=!!as.symbol(group_by_colname), color=!!as.symbol(group_by_colname))) +
      scale_color_manual(values = fill_colors[[groups_to_compare]]) +
      scale_fill_manual(values = ggplot2::alpha(fill_colors[[groups_to_compare]], 0.3)) +
      scale_x_discrete(expand = c(0.1, 0.1),
                       labels = group_lab) 
  } else {
    plot <- ggplot(d, aes(Group, !!as.name(FS_or_PFS))) +
      geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total")) +
      geom_quasirandom(size=3, shape = 16, width = 0.3, aes(color=!!as.symbol(group_by_colname))) +
      scale_color_manual(values = fill_colors[[groups_to_compare]]) +
      scale_x_discrete(expand = c(0.3, 0.3),
                       labels = group_lab) 
  }
  
  plot <- plot +
    labs(y = if(FS_or_PFS == "FS") {sprintf("%s Functionality Score", cd4_or_cd8)} else if(FS_or_PFS == "PFS") {sprintf("%s Polyfunctionality Score", cd4_or_cd8)},
         title = current_stim) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(color="black", size=17),
          axis.text.x = element_text(color="black", size=20),
          plot.title = element_text(hjust = 0.5, size=21),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(0.3, 0.2, 0.1, 0.2, "cm")) +
    force_panelsizes(rows = unit(3.5, "in"),
                     cols = unit(3, "in"))
  if(!is.null(ylim)) {
    plot <- plot + 
      coord_cartesian(ylim = ylim) +
      annotate("text", x = 1.5, y = ylim[2] - 0.09*diff(ylim),
               label = test_df$p_val_text, size=5.5)
  } else {
    plot_ylims <- ggplot_build(plot)$layout$panel_params[[1]]$y.range
    plot <- plot + 
      annotate("text", x = 1.5, y = plot_ylims[2] + 0.01*diff(plot_ylims),
               label = test_df$p_val_text, size=5.5) +
      coord_cartesian(ylim = c(plot_ylims[[1]], plot_ylims[[2]] + 0.09*diff(plot_ylims)))
  } 
}

############################################################################################################################

# Individual magnitude/frequency plots

make_mag_plots <- function(counts, counts_no_outlier = NULL, compare_time, keep, current_stim,
                           groups_to_compare, paired, fill_colors,
                           group_by_colname, subtitle, ylim = NULL, y_axis_text,
                           y_axis_size, axis_break = NULL) {
  if(compare_time) {
    counts <- counts %>%
      dplyr::filter(grepl(keep, Group) & Stim == current_stim)
  } else {
    counts <- counts %>%
      dplyr::filter(Timepoint == keep & Stim == current_stim)
  }
  
  if(!is.null(counts_no_outlier)) {
    if(compare_time) {
      counts_no_outlier <- counts_no_outlier %>%
        dplyr::filter(grepl(keep, Group) & Stim == current_stim)
    } else {
      counts_no_outlier <- counts_no_outlier %>%
        dplyr::filter(Timepoint == keep & Stim == current_stim)
    }
    counts_to_plot <- counts_no_outlier
  } else {
    counts_to_plot <- counts
  }
  
  # P-values are unadjusted
  if(paired) {
    # Signed-rank test
    test <- wilcox.test(Freq ~ Group, data = counts, paired = TRUE)
  } else {
    # Rank sum test
    test <- wilcox.test(Freq ~ Group, data = counts, paired = FALSE)
  }
  
  if(compare_time) {
    timepoint_infection_status_x_order <- 1:2
    names(timepoint_infection_status_x_order) <- groups_to_compare
    test_df <- data.frame(group1 = groups_to_compare[1],
                          group2 = groups_to_compare[2],
                          p.val = test$p.value) %>% 
      mutate(group1.xloc = unname(timepoint_infection_status_x_order[group1]),
             group2.xloc = unname(timepoint_infection_status_x_order[group2]),
             p_val_text = sapply(p.val, function(p) {
               if(p < 0.001) {
                 "p<0.001"
               } else {
                 paste0("p=", formatC(round(p, 3), format='f', digits=3))
               }
             }),
             geom_signif_group = paste0(group1, group2))
    
    medians <- counts %>% group_by(Group) %>%
      summarise(med = median(Freq)) %>% ungroup() %>% 
      mutate(x.min.segment = 1:2 - 0.1,
             x.end.segment = 1:2 + 0.1)
  } else {
    test_df <- data.frame(p = as.numeric(unlist(test)["p.value"])) %>%
      mutate(p_val_text = if_else(p < 0.001, "p<0.001", paste0("p=", formatC(round(p, 3), format='f', digits=3))))
  }
  
  if(compare_time) {
    current_plot <- ggplot(counts_to_plot, aes(x = Group, y = Freq)) +
      geom_line(aes(group = `PATIENT ID`), color = fill_colors[[keep]][1], size = 0.6, alpha = 0.5) +
      geom_segment(data = medians,
                   aes(x=x.min.segment, xend=x.end.segment, y=med, yend=med),
                   inherit.aes=FALSE, size = 1) +
      geom_point(size = 3, shape = 21, stroke = 1.5, aes(fill=!!as.symbol(group_by_colname), color=!!as.symbol(group_by_colname))) +
      scale_color_manual(values = fill_colors[[keep]]) +
      scale_fill_manual(values = ggplot2::alpha(fill_colors[[keep]], 0.3)) +
      scale_x_discrete(expand = c(0.1, 0.1), labels = stringr::str_replace_all(groups_to_compare, " ", "\n"))
  } else {
    current_plot <- ggplot(counts_to_plot, aes(x = Group, y = Freq)) +
      geom_boxplot(outlier.shape=NA, position = position_dodge2(preserve = "total")) +
      geom_quasirandom(size=3, shape = 16, width = 0.3, aes(color=!!as.symbol(group_by_colname))) +
      scale_color_manual(values = fill_colors[[keep]]) +
      scale_x_discrete(labels = stringr::str_replace_all(groups_to_compare, " ", "\n"), expand = c(0.3, 0.3))
  }
  
  current_plot <- current_plot +
    theme_bw() +
    theme(text = element_text(family="Arial"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = y_axis_size),
          axis.text.y = element_text(color="black", size=17),
          axis.text.x = element_text(color="black", size=20),
          plot.title = element_text(hjust = 0.5, size=21),
          plot.subtitle = element_text(hjust = 0.5, size=13),
          panel.grid.major.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(0.3, 0.2, 0.1, 0.2, "cm")) +
    labs(title = as.character(current_stim),
         subtitle = subtitle,
         y = y_axis_text) +
    force_panelsizes(rows = unit(3.5, "in"),
                     cols = unit(3, "in"))
    
  if(!is.null(ylim)) {
    current_plot <- current_plot + 
      coord_cartesian(ylim = ylim) +
      annotate("text", x = 1.5, y = ylim[2] - 0.09*diff(ylim),
               label = test_df$p_val_text, size=5.5)
  } else {
    plot_ylims <- ggplot_build(current_plot)$layout$panel_params[[1]]$y.range
    current_plot <- current_plot + 
      annotate("text", x = 1.5, y = plot_ylims[2] + 0.01*diff(plot_ylims),
               label = test_df$p_val_text, size=5.5) +
      coord_cartesian(ylim = c(plot_ylims[[1]], plot_ylims[[2]] + 0.09*diff(plot_ylims)))
  } 
}