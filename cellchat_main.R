# cellchat_main.R
# Author: Justin Reimertz

#' Main script for generation of CellChat analysis reports containing all
#' functions used in each Rmd file and libraries sourced

# Load packages
library(Seurat)
library(glue)
library(patchwork)
library(scCustomize)
library(qs)
library(raster)
library(beepr)
library(hdf5r)
library(ggbeeswarm)
library(ggrepel)
library(ggridges)
library(RColorBrewer)
library(ComplexHeatmap)
library(devtools)
library(Matrix)
library(umap)
library(magrittr)
library(data.table)
library(CellChat)
library(NMF)
library(ggalluvial)
library(VennDiagram)
library(gridExtra)
library(circlize)
library(migest)
library(CCPlotR)
library(tidyverse)

# Set future maxsize
options(future.globals.maxSize=1258291200)

##### Data Preprocessing #####

#' Function to retrieve a filtered Seurat object for the specified cell and 
#' tissue combination. If the Seurat file exists read in the file. If the file 
#' does not exist, apply the appropriate filters to the full Seurat object and
#' save the filtered object to the provided file path
#' 
#' @param meta_file (str): file path for meta data if required
#' @param sc_file (str): file path for filtered Seurat object if required
#' @param full_sc_file (str): file path for full Seurat object if required
#' @param group_label (str): column name to use when building cell groups 
#' @param tissues (vec): vector of tissues to include in CellChat analysis
#' @param cells (vec): vector of cells to include in CellChat analysis
#' @param timepoints (vec): Vector of timepoints to include in CellChat analysis
#'
#' @return returns a cellchat object that was either loaded in from the given
#' file path or that was built using the given specifications

get_filtered_data <- function(meta_file, sc_file, full_sc_file, group_label,
                              tissues=NULL, cells=NULL, timepoints=NULL) {
  # Read in the meta data file
  meta <- read_csv(meta_file, show_col_types = F)
  # Read in the Seurat object and apply filters to the full file if the filtered
  # file doesn't exist yet
  if (file.exists(sc_file)) {
    seurat_obj <- qread(sc_file)
  } else {
    # Read in the full Seurat object 
    full_seurat_obj <- qread(full_sc_file)
    
    # Filter based on tissue and/or cell type
    seurat_obj <- full_seurat_obj %>%
      label_clusters(meta, group_label) %>%
      apply_filters(tissues=tissues, cells=cells, timepoints=timepoints) %T>%
      qsave(sc_file)
  }
}

#' Function to label Seurat clusters using meta data information
#' 
#' @param seurat_obj (obj): Seurat object containing unlabeled clusters
#' @param meta (tibble): Meta data information for the seurat object including
#' cell types assigned to each cluster
#' @param group_label (str): column name to use when building cell groups 
#' 
#' @return a seurat object with clusters labeled by cell type

label_clusters <- function(seurat_obj, meta, group_label) {
  # Set meta data column names
  meta <- meta %>% 
    mutate(cluster = as.factor(Cluster)) %>%
    drop_na()
  
  # Add additional meta data and cell labels to Seurat object
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    # Save indents to a column
    rownames_to_column("indents") %>%
    # Join new meta data with Seurat meta data
    inner_join(meta, join_by(seurat_clusters == cluster)) %>%
    # rename columns to remove spaces
    rename(any_of(
      c(cell_labels = "Cluster labels", 
        lineage = "Cell class",
        cell_type = "Subtype")
      )) %>%
    # Restore indents as column names
    column_to_rownames("indents")
    #select(-c(anno_cellClass, anno_cluster))
  
  # Get the specified column to label cell groups
  group <- get(group_label, seurat_obj@meta.data)
  # Create a column in the Seurat object meta data with cell groups
  seurat_obj@meta.data$cell_groups <- group
  return(seurat_obj)
}


#' Function to apply filters to a seurat object depending on what tissue or
#' cell types should be included in further analysis
#' 
#' @param seurat_obj (obj): Seurat object containing labeled clusters
#' @param tissues (vec): vector of tissue types to include
#' @param cells (vec): vector of cell types to include
#' @param timepoints (vec): Vector of timepoints to include
#' 
#' @return Seurat object

apply_filters <- function(seurat_obj, tissues=NULL, cells=NULL, timepoints=NULL) {
  if (!is.null(tissues)) {
    Idents(seurat_obj) <- "Tissues"
    seurat_obj <- subset(seurat_obj, subset = Tissue %in% tissues)
  }
  if (!is.null(cells)) {
    Idents(seurat_obj) <- "lineage"
    seurat_obj <- subset(seurat_obj, subset = lineage %in% cells)
  }
  if (!is.null(timepoints)) {
    Idents(seurat_obj) <- "time"
    seurat_obj <- subset(seurat_obj, subset = time %in% timepoints)
  } 
  if (!is.null(phenotypes)) {
    Idents(seurat_obj) <- "phenotype"
    seurat_obj <- subset(seurat_obj, subset = phenotype %in% phenotypes)
  }
  Idents(seurat_obj) <- "seurat_clusters"
  return(seurat_obj)
}


##### Build Network #####

#' Function to create the CellChat object and set the appropriate database
#' 
#' @param seurat_obj (obj): Seurat object containing labeled clusters
#' @param group_label (str): column name to use when building cell groups 
#' @param sample_type (str): String to distinguish the sample type in order to
#' select the appropriate database
#' @param DB_subset (list): Named list specifying parameters to subset the
#' database by. Defaults to NULL which will result in using the entire database

build_cellchat <- function(seurat_obj, sample_type, DB_subset=NULL) {
  # Create the CellChat object using the Seurat object
  out <- createCellChat(object = seurat_obj, group.by = "cell_groups")
  # number of cells in each cell group
  groupSize <- as.numeric(table(out@idents))
  
  # Select DB to use based on given sample type
  if (sample_type == "mouse") {
    CellChatDB <- CellChatDB.mouse
  } else {
    CellChatDB <- CellChatDB.human
  }
  
  # Decide whether to use the full CellChatDB or to subset based on a given
  # parameter
  if (is.null(DB_subset)) {
    # use all CellChatDB for cell-cell communication analysis
    CellChatDB.use <- CellChatDB
  } else {
    # use a subset of CellChatDB for cell-cell communication analysis
    # use Secreted Signaling
    CellChatDB.use <- subsetDB(CellChatDB, 
                               search = DB_subset$search, 
                               key = DB_subset$key)
  }
  # set the used database in the object
  out@DB <- CellChatDB.use
  
  return(out)
}


#' Function to set up parallel computing for working with CellChat even if the
#' CellChat object was previously built
#' 
#' @param cellchat_obj (obj): cellchat object
#' @param threads (int): Number of cores to use during parallel computing
#' 
#' @return cellchat object with parallel workers allocated

set_parallel <- function(cellchat_obj, threads=4) {
  # subset the expression data of signaling genes for saving computation cost
  cellchat_obj <- subsetData(cellchat_obj)
  # do parallel processing using specified number of threads
  future::plan("multisession", workers = threads)
  cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
  cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
  
  return(cellchat_obj)
}

#' Function to infer the cell-cell communication network
#' 
#' @param cellchat_obj (obj): cellchat object
#' @param type (str): Type of cell-cell communication to infer for the given 
#' data. Defaults to "Tri-Mean"
#' @param min_cells (int): Minimum number of cells to filter cell-cell 
#' communication by. Defaults to 10

infer_comm_network <- function(cellchat_obj, type="triMean", min_cells=10) {
  #Compute the communication probability and infer cellular communication
  cellchat_obj <- computeCommunProb(cellchat_obj, type = type)
  # Filter cell-cell communication if there are only a few cells in cell groups
  cellchat_obj <- filterCommunication(cellchat_obj, min.cells = min_cells)
  # Infer the cell-cell communication at a signaling pathway level
  cellchat_obj <- computeCommunProbPathway(cellchat_obj)
  # Calculate the aggregated cell-cell communication network
  cellchat_obj <- aggregateNet(cellchat_obj)
  # Compute the network centrality scores
  cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
  
  return(cellchat_obj)
}



##### Statistical Analysis #####

#' Function to set up the appropriate pairwise distance comparisons
#' 
#' @param df (df): df containing x and y coordinates for each cell type and
#' time point information
#' @param index (int): current index in list of timepoints
#' 
#' @return tibble

build_pairwise_comparison <- function(df, index) {
  subset_time_vec <- unique(df$time)[index:length(unique(df$time))]
  
  out <- bind_rows(lapply(c(2:length(subset_time_vec)), function(i) {
    
    sub_df1 <- df %>% filter(time == subset_time_vec[1]) %>%
      dplyr::select(x,y)
    sub_df2 <- df %>% filter(time == subset_time_vec[i]) %>%
      dplyr::select(x,y)
    
    get_euc_dist(df1 = sub_df1,
                 df2 = sub_df2,
                 comparison = glue("{subset_time_vec[1]}_{subset_time_vec[i]}"))
  }))
  return(out)
}

#' Function to calculate pairwise euclidean distances between points when 
#' comparing multiple datasets
#' 
#' @param df1 (df): first dataframe to include 
#' @param df2 (df): second dataframe to include
#' @param comparison (str): label for the two datasets being compared
#' 
#' @return tibble

get_euc_dist <- function(df1, df2, comparison) {
  dist <- pointDistance(df1, df2, lonlat = F)
  
  out <- tibble(cell_types = names(dist), 
                distances = dist, 
                time = comparison)
  
  return(out)
}


##### Visualize Network #####

#' Function to plot aggregated cell-cell communication network
#' 
#' @param cellchat (obj): CellChat object
#' @param pathway (vec): vector of pathways to plot
#' @param type (str): string to specify type of plot to generate
#' 
#' @return circle plots of the aggregated cell-cell communication network for
#' the given CellChat object

plot_netVis <- function(cellchat, pathway=cellchat@netP$pathways, type="circle") {
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathway, layout = type)
}


#' Function to plot the individual cell contributions for the cell-cell
#' communication network as a circle plot
#' 
#' @param cellchat (obj): CellChat object
#' 
#' @return circle plots for each individual cell present in the aggregated
#' cell-cell communication network

plot_netVis_circlInd <- function(cellchat) {
  mat <- cellchat@net$weight
  par(mfrow = c(2,2), xpd=TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(
      mat2, vertex.weight = groupSize, weight.scale = T, 
      edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}

#' Function to plot cell contributions in a single pathway as a heatmap
#' 
#' @param cellchat (obj): CellChat object
#' @param pathway (vec): vector with specific pathway to visualize
#' @param color (str): heatmap color scheme to use 
#' 
#' @return heatmap of individual cell contributions for a given pathway

plot_netVis_heat <- function(cellchat, pathway, color) {
  par(mfrow=c(1,1))
  netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
}


#' Function to plot bubble plot of the significant interactions (L-R pairs) 
#' from some cell groups (defined by 'sources.use') to other cell groups 
#' (defined by 'targets.use')
#' 
#' @param cellchat (obj): CellChat object
#' @param source (vec): vector of cell group indexes for the sources of 
#' cell-cell interactions
#' @param target (vec): vector of cell group indexes for the targets of
#' cell-cell interactions
#' @param save_plot (bool): Boolean value to determine whether or not to save
#' the generate plot. Defaults to `FALSE`
#' @param plot_file (str): filename to save the generated plot to. Defaults to
#' `NULL`
#' @param height (int): Specified height to save the generated plot to. Defaults
#' to `NULL`
#' @param width (int): Specified width to save the generated plot to. Defaults
#' to `NULL`
#' @param device (str): File type to save the plot as. Defaults to `"pdf"`
#' 
#' @return bubble plot of significant L-R pairs from specified sources to the
#' specified targets

plot_netVis_bubble <- function(cellchat, source, target, save_plot=F, 
                               plot_file=NULL, height=NULL, width=NULL, 
                               device="pdf") {
  p <- netVisual_bubble(cellchat, sources.use = source, targets.use = target, 
                   remove.isolate = F)
  if (save_plot) {
    ggsave(plot = p, path = "../plots", device = device, 
           filename = plot_file, height = height, width = width) 
  }
  
  return(p)
}


#' Function adapted from CellChat::netAnalysis_signalingRole_scatter() for 2D 
#' visualization of dominant senders (sources) and receivers (targets)
#'
#' @description
#' This function builds a dataframe of the dominant senders (sources) and receivers 
#' (targets) for a single CellChat object so that they can be plotted in 2D space
#'
#' @param cellchat (obj): CellChat object
#' @param signaling (vec): char vector containing signaling pathway names. 
#' signaling = NULL: signaling role analysis on the aggregated cell-cell 
#' communication network from all signaling pathways
#' @param slot.name the slot name of object that is used to compute centrality 
#' measures of signaling networks. Defaults to "netP"
#' @param x.measure The measure used as x-axis. This measure should be one of 
#' `names(slot(cellchat, slot.name)$centr[[1]])` computed from 
#' `netAnalysis_computeCentrality`
#'
#' Default = "outdeg" is the weighted outgoing links (i.e., outgoing interaction 
#' strength). If setting as "outdeg_unweighted", it represents the total number 
#' of outgoing signaling.
#'
#' @param y.measure The measure used as y-axis. This measure should be one of 
#' `names(slot(cellchat, slot.name)$centr[[1]])` computed from 
#' `netAnalysis_computeCentrality`
#'
#' Default = "indeg" is the weighted incoming links (i.e., incoming interaction 
#' strength). If setting as "indeg_unweighted", it represents the total number 
#' of incoming signaling.
#' 
#' @return dataframe
#' @export
#'
build_netAnalysis_signalingRole_df <- function(cellchat, 
                                              signaling = NULL, 
                                              slot.name = "netP", 
                                              x.measure = "outdeg", 
                                              y.measure = "indeg") {
  
  if (length(slot(cellchat, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  if (sum(c(x.measure, y.measure) %in% names(slot(cellchat, slot.name)$centr[[1]])) !=2) {
    stop(paste0("`x.measure, y.measure` should be one of ", 
                paste(names(slot(cellchat, slot.name)$centr[[1]]),collapse=", "), 
                '\n', "`outdeg_unweighted` is only supported for version >= 1.1.2"))
  }
  centr <- slot(cellchat, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(cellchat@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(cellchat@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(cellchat@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[,i] <- centr[[i]][[x.measure]]
    incoming[,i] <- centr[[i]][[y.measure]]
  }
  if (is.null(signaling)) {
    message("Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways")
  } else {
    message("Signaling role analysis on the cell-cell communication network from user's input")
    signaling <- signaling[signaling %in% cellchat@netP$pathways]
    if (length(signaling) == 0) {
      stop('There is no significant communication for the input signaling. All the significant signaling are shown in `cellchat@netP$pathways`')
    }
    outgoing <- outgoing[ , signaling, drop = FALSE]
    incoming <- incoming[ , signaling, drop = FALSE]
  }
  outgoing.cells <- rowSums(outgoing)
  incoming.cells <- rowSums(incoming)
  
  num.link <- aggregateNet(cellchat, signaling = signaling, return.object = FALSE,
                           remove.isolate = FALSE)$count
  num.link <- rowSums(num.link) + colSums(num.link)-diag(num.link)
  df <- data.frame(x = outgoing.cells, y = incoming.cells, 
                   labels = names(incoming.cells),
                   Count = num.link)
  df$labels <- factor(df$labels, levels = names(incoming.cells))
  return(df)
}


#' Function adapted from CellChat::netAnalysis_signalingRole_scatter() for 2D 
#' visualization of dominant senders (sources) and receivers (targets)
#'
#' @description
#' This scatter plot shows the dominant senders (sources) and receivers 
#' (targets) in a 2D space. x-axis and y-axis are respectively the total outgoing 
#' or incoming communication probability associated with each cell group.
#' Dot size is proportional to the number of inferred links (both outgoing and 
#' incoming) associated with each cell group.
#' Dot colors indicate different cell groups. Dot shapes indicate different 
#' categories of cell groups if `group`` is defined.
#'
#' @param df (df): dataframe returned from `build_netAnalysis_signalingRole_df()`
#' @param color.use defining the color for each cell group
#' @param group (vec): a vector to categorize the cell groups, e.g., categorize the 
#' cell groups into two major categories: immune cells and fibroblasts
#' @param weight.MinMax the Minmum/maximum weight, which is useful to control 
#' the dot size when comparing multiple datasets
#' @param point.shape point shape when group is not NULL
#' @param label.size font size of the text
#' @param dot.alpha transparency
#' @param dot.size a range defining the size of the symbol
#' @param xlabel label of x-axis
#' @param ylabel label of y-axis
#' @param title (str): main title of the plot
#' @param font.size (int): font size of the text
#' @param font.size.title (int): font size of the title
#' @param guide_color whether or not to include labels in legend
#' @param do.label (bool): label the each point
#' @param show.legend (bool): whether show the legend
#' @param show.axes (bool): whether show the axes
#' @param axis.lims (list): named list specifying the upper and lower limits of
#' the x and y axes
#' @return ggplot object
#' @export
#'
build_netAnalysis_signalingRole_scatter <- function(df, 
                                               color.use = NULL, 
                                               group = NULL, 
                                               weight.MinMax = NULL, 
                                               dot.size = c(2, 6), 
                                               point.shape = c(21, 22, 24, 23, 25, 8, 3), 
                                               label.size = 3, 
                                               dot.alpha = 0.6,
                                               xlabel = "Outgoing interaction strength", 
                                               ylabel = "Incoming interaction strength", 
                                               title = NULL,
                                               font.size = 10, 
                                               font.size.title = 10, 
                                               guide_color = FALSE,
                                               do.label = T, 
                                               show.legend = T, 
                                               show.axes = T,
                                               axis.lims = list(x=c(0,30), y=c(0,30))) {
  
  if (!is.null(group)) {
    df$Group <- group
  }
  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(cellchat@idents))
  }
  if (!is.null(group)) {
    gg <- ggplot(data = df, aes(x, y)) +
      geom_point(aes(size = Count, colour = labels, fill = labels, shape = Group))
  } else {
    gg <- ggplot(data = df, aes(x, y)) +
      geom_point(aes(size = Count, colour = labels, fill = labels))
  }
  
  gg <- gg + geom_abline(aes(slope = 1, intercept = 0, alpha = 0.5), show.legend = F) + 
    CellChat_theme_opts() + theme(text = element_text(size = font.size), 
          legend.key.height = grid::unit(0.15, "in")) +
    # guides(colour = guide_legend(override.aes = list(size = 3)))+
    labs(title = title, x = xlabel, y = ylabel) + 
    theme(plot.title = element_text(size= font.size.title, face="plain")) +
    # theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
    theme(axis.line.x = element_line(size = 0.25), 
          axis.line.y = element_line(size = 0.25))
  gg <- gg + scale_fill_manual(
    values = ggplot2::alpha(color.use, alpha = dot.alpha), drop = FALSE) + 
    guides(fill=FALSE)
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE) + 
    guides(colour=guide_color) + lims(x=axis.lims$x, y=axis.lims$y)
  # gg <- gg + scale_colour_manual(values = ggplot2::alpha(color.use, alpha = dot.alpha), 
  # drop = FALSE) + guides(colour=FALSE)
  # gg <- gg + scale_shape_manual(values = point.shape[1:length(prob)])
  if (!is.null(group)) {
    gg <- gg + scale_shape_manual(values = point.shape[1:length(unique(df$Group))])
  }
  if (is.null(weight.MinMax)) {
    gg <- gg + scale_size_continuous(range = dot.size)
  } else {
    gg <- gg + scale_size_continuous(limits = weight.MinMax, range = dot.size)
  }
  if (do.label) {
    gg <- gg + ggrepel::geom_text_repel(
      mapping = aes(label = cell_labels, colour = cell_labels), size = label.size, 
      show.legend = F,segment.size = 0.2, segment.alpha = 0.5, box.padding = 3,
      max.overlaps = 60)
  }
  
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  
  if (!show.axes) {
    gg <- gg + theme_void()
  }
  
  gg
  
}



#' Function to build a chord diagram using the CCPlotR package
#' 
#' @param cc_net (df): Dataframe network generated by CellChat
#' @param avg_exp (tibble/file): Tibble with average expression values for each 
#' gene and cell_group as denoted in the Seurat object used for CellChat 
#' analysis or filepath where the average expression can be loaded in from
#' @param plot_type (str): Value to provide to `cc_circos` to determine the type
#' of chord diagram to plot. Possible values include: `"A"`, `"B"`, `"C"`, `"D"`,
#' or `"E"`
#' @param cell_vec1 (vec): Vector of cell types to include
#' @param cell_vec2 (vec): Vector of cell types to include
#' @param interactions_vec (vec): Vector of ligand-receptor reactions to subset
#' `cc_net` by
#' @param direction (str): String to denote whether interactions should include
#' cell_vec1 - cell_vec2 (`"forward"`), cell_vec2 - cell_vec1 (`"reverse"`), 
#' or both ways (`"both"`). Defaults to `"forward"`
#' @param top_n_rows (int): Number of rows to subset the interactions scores by as
#' required for `CCPlotR::cc_circos()`. Defaults to `1000`
#' @param seurat_obj (obj): Seuart object used for CellChat analysis. Only 
#' required if `avg_exp` needs to be calculated. Defaults to `NULL`
#' @param exp_group (str): Column in `seurat_obj` to calculate average expression by.
#' Only required if `avg_exp` needs to be calculated. Defaults to `NULL`
#' @param special_char (bool): Boolean value to determine whether or not special
#' character values in cell types should be removed. Defaults to `FALSE`
#' @param color_list (vec): Named list of cell types and corresponding colors.
#' Defaults to `NULL`
#' @param show_legend (bool): Boolean value to determine whether or not the 
#' legend should be displayed with the chord diagram output. Defaults to `FALSE`
#' @param cex (int): cex value to determine the text size for the chord diagram.
#' Defaults to `1`
#' @param subtype_buffer (int): Value to add to size of cell subtype boxes when
#' building chord diagram
#' @param big_gap (int): Space between groups in chord diagram. Correlates with
#' `big.gap` from `circlize::chordDiagram`
#' @param small_gap (int): Space between groups in chord diagram. Correlates with
#' `small.gap` from `circlize::chordDiagram`
#' @param palette (str): Which colour palette to use to show the mean expression. 
#' Should be one of the RColorBrewer sequential palettes.
#' 
#' @return chord diagram

build_chord_diagram <- function(cc_net, avg_exp, plot_type, cell_vec1, cell_vec2, 
                                interactions_vec, direction = "forward", 
                                top_n_rows = 1000, seurat_obj = NULL, 
                                exp_group = NULL, special_char = FALSE, 
                                color_list = NULL, show_legend = FALSE,
                                cex = 1, subtype_buffer = 0.1, 
                                big_gap = 8, small_gap = 1, palette = "BuPu") {
  # Get the average expression for each gene and cell type
  if(!is_tibble(avg_exp)) {
    avg_exp <- calc_avg_expression(avg_exp, seurat_obj, group)
  }
  
  # Set up the CellChat network to the correct format and subset the appropriate
  # cell types and interactions
  if (direction %in% c("forward", "both")) {
    interactions1 <- cc_net %>%
      rename(score = prob) %>%
      filter(str_detect(source, paste(cell_vec1, collapse = "|")),
             str_detect(target, paste(cell_vec2, collapse = "|")),
             interaction_name_2 %in% interactions_vec) %>%
      mutate(source = str_remove(source, "\\+"),
             target = str_remove(target, "\\+"))
  }
  
  if (direction %in% c("reverse", "both")) {
    interactions2 <- cc_net %>%
      rename(score = prob) %>%
      filter(str_detect(target, paste(cell_vec1, collapse = "|")),
             str_detect(source, paste(cell_vec2, collapse = "|")),
             interaction_name_2 %in% interactions_vec) %>%
      mutate(source = str_remove(source, "\\+"),
             target = str_remove(target, "\\+"))
  }
  
  if (direction == "both") {
    interactions_all <- bind_rows(interactions1, interactions2)
  } else {
    if (direction == "forward") {
      interactions_all <- interactions1
    } else if (direction == "reverse") {
      interactions_all <- interactions2
    }
  }
  
  # Account for if any ligand/receptor pairs include gene complexes
  interactions_all <- interactions_all %>%
    separate_wider_regex(interaction_name_2, 
                         patterns = c("^.+- \\(*",
                                      receptor_1 = ".+",
                                      "\\+",
                                      receptor_2 = ".+",
                                      "\\)$"),
                         too_few = "align_start",
                         cols_remove = F) %>%
    separate_wider_regex(interaction_name_2, 
                         patterns = c(ligand_1 = "^.+",
                                      "  -.*$"),
                         too_few = "align_start",
                         cols_remove = F) %>%
    separate_wider_regex(ligand_1, 
                         patterns = c("^\\(*",
                                      ligand_1 = ".+",
                                      "\\+",
                                      ligand_2 = ".+",
                                      "\\)$"),
                         too_few = "align_start") %>%
    mutate(receptor = ifelse(!is.na(receptor_2), 
                             paste0(receptor_1," ", receptor_2),
                             receptor),
           ligand = ifelse(!is.na(ligand_2), 
                           paste0(ligand_1," ", ligand_2),
                           ligand)) %>%
    pivot_longer(cols = c(receptor_1, receptor_2), names_to = "receptor_num",
                 values_to = "single_receptor", values_drop_na = T) %>%
    pivot_longer(cols = c(ligand_1, ligand_2), names_to = "ligand_num",
                 values_to = "single_ligand", values_drop_na = T)
  
  # Create groups for the general cell types and appropriate subtypes
  cell_groups <- interactions_all %>%
    pivot_longer(cols = c(source, target), names_to = "direction", 
                 values_to = "cell_type") %>%
    separate_wider_delim(cell_type, "_", names = c("cell_group", "cluster"),
                         too_few = "align_start", cols_remove = F) %>%
    select(cell_type, cell_group) %>%
    distinct()
  
  # Assign cell subtype identities specific colors
  cell_colors <- lapply(names(color_list), function (c) {
    if (c %in% cell_groups[["cell_group"]]) {
      # Determine the number of clusters for the given cell_type
      clusters <- cell_groups %>%
        filter(cell_group == c) %>%
        pull(cell_type) %>%
        unique() 
      # Stop if not enough colors for the number of clusters
      stopifnot(length(color_list[[c]]) >= length(clusters))
      # Subset the colors vector for the number of clusters and then set the names
      # to the cluster ids
      set_names(colors[[c]][1:length(clusters)], clusters)
    }
  }) %>% unlist()
  
  # Generate the chord diagram
  if (plot_type %in% c("A", "B", "C")) {
    out <- cc_circos(interactions_all, option = plot_type, 
                     n_top_ints = top_n_rows, exp_df = avg_exp, 
                     cell_cols = cell_colors, show_legend = show_legend)
  } else if (plot_type %in% c("D", "E")) {
    out <- build_cc_circos(interactions_all, n_top_ints = top_n_rows, 
                           colors_vec = cell_colors, avg_exp = avg_exp,
                           plot_type = plot_type, show_legend = show_legend,
                           cex = cex, subtype_buffer = subtype_buffer,
                           big_gap = big_gap, small_gap = small_gap, 
                           palette = palette)
  }
  return(out)
}


#' Function adapted from CCPlotR::cc_circos() to make a chord diagram with the
#' same general layout but split by ligand/receptor genes rather than cell type
#' 
#' @param interactions_df (tibble): Tibble of subsetted interactions to include with
#' source, target, and respective score information
#' @param n_top_ints (int): Number of rows to subset interactions_df rows by
#' @param plot_type (str): Value to provide to `circlize_plot()` to determine 
#' the type of chord diagram to plot. Possible values include: `"D"` or `"E"`
#' @param colors_vec (vec): Named vector of colors where names are cell subtypes
#' and values are color assignments
#' @param avg_exp (tibble): Tibble with average expression values for each 
#' gene and cell_group as denoted in the Seurat object used for CellChat 
#' analysis. 
#' @param include_exp (bool): Boolean to denote whether or not expression 
#' information should be added to the chord diagram. Defaults to `FALSE`
#' @param special_char (bool): Boolean value to determine whether or not special
#' character values in cell types should be removed. Defaults to `FALSE`
#' @param show_legend (bool): Boolean value to determine whether or not the 
#' legend should be displayed in the final chord diagram output. Defaults to 
#' `FALSE`
#' @param cex (int): cex value to determine the text size for the chord diagram.
#' Defaults to `1`
#' @param subtype_buffer (int): Value to add to size of cell subtype boxes when
#' building chord diagram
#' @param big_gap (int): Space between groups in chord diagram. Correlates with
#' `big.gap` from `circlize::chordDiagram`
#' @param small_gap (int): Space between groups in chord diagram. Correlates with
#' `small.gap` from `circlize::chordDiagram`
#' @param palette (str): Which colour palette to use to show the mean expression. 
#' Should be one of the RColorBrewer sequential palettes.
#' 
#' @return chord diagram

build_cc_circos <- function(interactions_df, n_top_ints, plot_type, colors_vec, 
                            avg_exp, include_exp = FALSE, special_char = FALSE,
                            show_legend = FALSE, cex = 1, subtype_buffer = 0.1,
                            big_gap = 8, small_gap = 1, palette = "BuPu") {
  # Create a dataframe with columns for each cell type and gene pairing
  input_df <- interactions_df %>%
    slice_max(order_by = score, n = n_top_ints) %>%
    mutate(source = str_remove(source, "\\+"),
           target = str_remove(target, "\\+")) %>%
    separate_wider_delim(source, "_", names = c("source_group", "source_cluster"),
                         too_few = "align_start", cols_remove = F) %>%
    separate_wider_delim(target, "_", names = c("target_group", "target_cluster"),
                         too_few = "align_start", cols_remove = F) %>%
    mutate(
      source_lig = paste0(source, "|", ligand),
      target_rec = paste0(target, "|", receptor),
      source_group_lig = paste0(source_group, "|", ligand),
      target_group_rec = paste0(target_group, "|", receptor)
    ) %>%
    arrange(source)
  # Create weights for linking arrows
  arr_wd <- (((input_df$score - min(input_df$score)
               ) / (max(input_df$score) - min(input_df$score))) * (4)) + 1
  
  # Assign colors for linking arrows
  link_cols <- c()
  for (i in input_df$source_lig) {
    link_cols <- c(link_cols, colors_vec[str_extract(i, "[^|]+")])
  }
  
  # Create segments for chord diagrams using ligands and receptors
  segments <- unique(c(paste0(input_df$source, "|", input_df$ligand), 
                       paste0(input_df$target, "|", input_df$receptor)))
  # Extract the cell names from the segments in order to make the group
  grp <- ifelse(str_extract(segments, "[^|]+") %in% Bcells, 
         paste0("B cells", "|", str_extract(segments, "[^|]+$")), 
         paste0("stroma", "|", str_extract(segments, "[^|]+$")))
  #grp <- str_extract(segments, "[^|]+$")
  names(grp) <- segments
  
  # Use expression info to create a gene df to assign inner segment colors
  gene_df <- as.data.frame(
    bind_rows(combine_mean_exp(input_df, avg_exp, "ligand"),
              combine_mean_exp(input_df, avg_exp, "receptor")) %>% 
      separate_wider_delim(cell_type, "_", names = c("cell_group", "cluster"),
                           too_few = "align_start", cols_remove = F) %>%
      mutate(cell_gene = paste0(cell_type, "|", gene),
             cell_group_gene = paste0(cell_group, "|", gene)) %>% 
      filter(cell_gene %in% segments)
  )
  rownames(gene_df) <- gene_df$cell_gene

  
  # Assign colors for each cell gene pairing segment
  inner_cols <- colors_vec %>%
    enframe(name = "cell_type", value = "colors") %>%
    merge(gene_df, by = "cell_type") %>%
    select(cell_gene, colors) %>%
    distinct() %>%
    deframe()
  
  # Create named vectors for assigning gaps
  subtype_gaps <- segments %>%
    enframe() %>%
    separate_wider_delim(cols = value, delim = "|", 
                         names = c("celltype", "gene")) %>%
    mutate(celltype_group = str_extract(celltype, "[^_]+")) %>%
    group_by(gene, celltype_group) %>%
    mutate(gap = ifelse(row_number() == n(), 2, 0)) %>%
    pull(gap)
  names(subtype_gaps) <- segments
  
  gene_gaps <- str_extract(segments, "[^|]$")
  names(gene_gaps) <- segments
  
  brks <- scales::pretty_breaks(n = 5)(c(floor(min(gene_df$mean_exp)), 
                                         ceiling(max(gene_df$mean_exp))))
  gene_col_fun <- colorRamp2(brks, RColorBrewer::brewer.pal(length(brks), palette))
  
  exp_inner_cols <- setNames(gene_col_fun(gene_df[segments, "mean_exp"]), segments)
  
  
  lgd1 <- Legend(
    labels = unique(c(input_df$source, input_df$target)),
    title = "Cell type",
    type = "points",
    title_gp = gpar(fontsize = 14 * 1),
    labels_gp = gpar(fontsize = 12 * 1),
    legend_gp = gpar(col = "transparent"),
    background = colors_vec[unique(c(input_df$source, input_df$target))],
    direction = "horizontal"
  )
  
  lgd2 <- Legend(
    title_gp = gpar(fontsize = 14 * 1),
    labels_gp = gpar(fontsize = 12 * 1),
    direction = "horizontal", at = brks,
    col_fun = gene_col_fun, title = "Mean exp."
  )
  
  if (plot_type == "D") {
    # Function to build the chord diagram
    circlize_plot <- function() {
      par(cex = cex)
      circos.par(start.degree = 90, points.overflow.warning = F, 
                 gap.after = subtype_gaps)
      chordDiagram(
        input_df %>%
          select(source_lig, target_rec, score),
        directional = 1, group = grp, link.sort = FALSE, diffHeight = 0.005, 
        scale = F, direction.type = c("arrows"), link.arr.type = "triangle", 
        annotationTrack = c(), big.gap = big_gap, small.gap = small_gap, 
        transparency = 1, 
        preAllocateTracks = list(list(track.height = 0.15), 
                                 list(track.height = 0.01), 
                                 list(track.height = 0.045)),
        link.arr.lwd = arr_wd, link.arr.col = link_cols, 
        link.arr.length = 0.4, link.arr.width = 0.35
      )
      
      for (g in unique(str_extract(segments, "[^_|]+$"))) {
        highlight.sector(segments[str_detect(segments, paste0("\\|", g, "$"))], 
                         track.index = 1, col = "white", text = g, cex = 1.3, 
                         text.col = "black", niceFacing = T)
      }
      
      circos.track(track.index = 2, panel.fun = function(x,y) {
        #for (l in unique(str_extract(segments, "[^_|]+$"))) {
          #highlight.sector(segments[str_detect(segments, paste0("\\|", l, "$"))], 
          #track.index = 2, col = "black") 
          
        #}
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        sector_name <- get.cell.meta.data("sector.index")
        sector_index <- get.cell.meta.data("sector.numeric.index")
        cell_width <- get.cell.meta.data("cell.width")
        xplot <- get.cell.meta.data("xplot")
        track_index <- get.current.track.index()
        

        sectors_list <- lapply(unique(str_extract(
          get.cell.meta.data("sector.index"), "[^|]+$")), function(g) {
            #print(g)
            get.cell.meta.data("sector.numeric.index")[
              str_extract(get.cell.meta.data("sector.numeric.index"),
                          "[^|]+$") == g
            ]
          })
        names(sectors_list) <- unique(str_extract(
          get.cell.meta.data("sector.index"), "[^|]+$"))
        #print(names(sectors_list))
        #print(sectors_list)
        #print(track_index)
        #print(sector_index)
        #print(sector_name[[1]])
        #print(get.cell.meta.data("xlim"))
        
        circos.rect(CELL_META$xlim[1]-subtype_buffer*1.5, 
                    CELL_META$ylim[1],
                    CELL_META$xlim[2]+subtype_buffer*1.5, 
                    CELL_META$ylim[2],
                    sector.index = CELL_META$sector.index, col = "black")
      })
        
        circos.track(track.index = 3, panel.fun = function(x, y) {
          
          gap_buff <- subtype_gaps[subtype_gaps == CELL_META$sector.index]
          
          circos.rect(mean(CELL_META$xcenter)-subtype_buffer, 
                      CELL_META$ylim[1], 
                      mean(CELL_META$xcenter)+subtype_buffer, 
                      1, 
                      sector.index = CELL_META$sector.index, 
                      col = inner_cols[CELL_META$sector.index]
                      )
          }, #cell.padding = c(gap_buff, 0, gap_buff, 0), 
          bg.border = NA)
      if (show_legend == TRUE) {
        draw(packLegend(lgd1, direction = "vertical"), 
             just = c("left", "bottom"), x = unit(4.75, "mm"), 
             y = unit(4.75, "mm"))
      }
      circos.clear()
    }
    out <- circlize_plot()
  } else if (plot_type == "E") {
    # Function to build the chord diagram
    circlize_plot <- function() {
      par(cex = cex)
      circos.par(start.degree = 90, points.overflow.warning = F, 
                 gap.after = subtype_gaps)
      chordDiagram(
        input_df %>%
          select(source_lig, target_rec, score),
        directional = 1, group = grp, link.sort = FALSE, diffHeight = 0.005, 
        scale = F, direction.type = c("arrows"), link.arr.type = "triangle", 
        annotationTrack = c(), big.gap = 8, transparency = 1, 
        preAllocateTracks = list(list(track.height = 0.15), 
                                 list(track.height = 0.01), 
                                 list(track.height = 0.045),
                                 list(track.height = 0.045)),
        link.arr.lwd = arr_wd, link.arr.col = link_cols, 
        link.arr.length = 0.4, link.arr.width = 0.35
      )
      
      for (g in unique(str_extract(segments, "[^_|]+$"))) {
        highlight.sector(segments[str_detect(segments, paste0("\\|", g, "$"))], 
                         track.index = 1, col = "white", text = g, cex = 1.3, 
                         text.col = "black", niceFacing = T)
      }
      
      circos.track(track.index = 2, panel.fun = function(x,y) {
        #for (l in unique(str_extract(segments, "[^_|]+$"))) {
        #highlight.sector(segments[str_detect(segments, paste0("\\|", l, "$"))], 
        #track.index = 2, col = "black") 
        
        #}
        xlim <- get.cell.meta.data("xlim")
        ylim <- get.cell.meta.data("ylim")
        sector_name <- get.cell.meta.data("sector.index")
        sector_index <- get.cell.meta.data("sector.numeric.index")
        cell_width <- get.cell.meta.data("cell.width")
        xplot <- get.cell.meta.data("xplot")
        track_index <- get.current.track.index()
        
        
        sectors_list <- lapply(unique(str_extract(
          get.cell.meta.data("sector.index"), "[^|]+$")), function(g) {
            #print(g)
            get.cell.meta.data("sector.numeric.index")[
              str_extract(get.cell.meta.data("sector.numeric.index"),
                          "[^|]+$") == g
            ]
          })
        names(sectors_list) <- unique(str_extract(
          get.cell.meta.data("sector.index"), "[^|]+$"))
        #print(names(sectors_list))
        #print(sectors_list)
        #print(track_index)
        #print(sector_index)
        #print(sector_name[[1]])
        #print(get.cell.meta.data("xlim"))
        
        circos.rect(CELL_META$xlim[1]-subtype_buffer*1.5, 
                    CELL_META$ylim[1],
                    CELL_META$xlim[2]+subtype_buffer*1.5, 
                    CELL_META$ylim[2],
                    sector.index = CELL_META$sector.index, col = "black")
      })
      
      circos.track(track.index = 3, panel.fun = function(x, y) {
        
        gap_buff <- subtype_gaps[subtype_gaps == CELL_META$sector.index]
        
        circos.rect(mean(CELL_META$xcenter)-subtype_buffer, 
                    CELL_META$ylim[1], 
                    mean(CELL_META$xcenter)+subtype_buffer, 
                    1, 
                    sector.index = CELL_META$sector.index, 
                    col = inner_cols[CELL_META$sector.index]
        )
      }, #cell.padding = c(gap_buff, 0, gap_buff, 0), 
      bg.border = NA)
      
      circos.track(track.index = 4, panel.fun = function(x, y) {
        circos.rect(CELL_META$xcenter-subtype_buffer, 
                    CELL_META$ylim[1], 
                    CELL_META$xcenter+subtype_buffer, 
                    1, 
                    sector.index = CELL_META$sector.index, 
                    col = exp_inner_cols[CELL_META$sector.index]
        )
      }, bg.border = NA)
      if (show_legend == TRUE) {
        draw(packLegend(lgd1, lgd2, direction = "vertical"), 
             just = c("left", "bottom"), x = unit(4.75, "mm"), 
             y = unit(4.75, "mm"))
      }
      circos.clear()
    }
    out <- circlize_plot()
    print(out)
  }
  
  return(out)
}

#' Function to calculate the average expression of each gene for each cell type
#' in a given Seurat object
#' 
#' @param seurat_obj (obj): Seuart object used for CellChat analysis
#' @param group (str): Column in `seurat_obj` to calculate average expression by
#' @param exp_file (file): Filepath to look for average_expression if it has
#' previously been calculated or to save the results to
#' 
#' @return Tibble with average expression results as expected by 
#' `build_chord_diagram()`

calc_avg_expression <- function(seurat_obj, group, exp_file) {
  if (!file.exists(exp_file)) {
    # Calculate the average expression for each gene and cell type
    out <- AverageExpression(seurat_obj, group.by = group) %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(cols = -gene, names_to = "cell_type", 
                   values_to = "mean_exp") %>%
      mutate(cell_type = str_remove(cell_type, "RNA\\."),
             cell_type = str_replace_all(cell_type, "\\.", "_")) %T>%
      write_delim(exp_file, delim = ",")
  } else {
    out <- read_delim(exp_file, show_col_types = F)
  }
  return(out)
}


#' Helper function to combine the mean expression values for ligands or receptors
#' that belong to a gene complex
#' 
#' @param interactions_df (tibble): Tibble of subsetted interactions to include 
#' with source, target, and respective score information
#' @param avg_exp (tibble): Tibble with average expression values for each 
#' gene and cell_group as denoted in the Seurat object used for CellChat 
#' analysis. 
#' @param signal (str): String to indicate whether the gene complex is a ligand
#' or receptor
#' 
#' @return Tibble of mean gene expression values for each gene and cell group

combine_mean_exp <- function(interactions_df, avg_exp, signal) {
  # Check that signal is valid
  stopifnot(signal %in% c("ligand", "receptor"))
  
  # Set up column variables for the given signal
  sig_num <- glue("{signal}_num")
  single_sig <- glue("single_{signal}")
  
  # Combine the mean expression for all genes in the complex
  out <- avg_exp %>%
    # Merge the relevant interactions_df columns to get the genes that are
    # part of a complex
    inner_join(select(interactions_df, 
                      c({{signal}}, {{sig_num}}, {{single_sig}})) %>%
                 distinct(), 
               join_by(gene == {{single_sig}})) %>%
    # Combine the mean expressions for each ligand/receptor and cell type
    group_by(cell_type, {{signal}}) %>%
    mutate(mean_exp = mean(mean_exp),
           gene = get(signal)) %>%
    ungroup() %>%
    select(gene, cell_type, mean_exp) %>%
    distinct()
  
  return(out)
}


