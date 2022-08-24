rm(list=ls())
gc()

ReQ_packages = c("dplyr", "tidyr", "magrittr", "data.table", "ggpubr", "prismatic", "ggplot2", 
                 "RColorBrewer", "pheatmap", "VennDiagram")

for (pack in ReQ_packages) {
  if(pack %in% rownames(installed.packages()) == FALSE) {
    BiocManager::install(pack)
    install.packages(pack)
    suppressPackageStartupMessages(library(pack, character.only = TRUE))
  } else {
    suppressPackageStartupMessages(library(pack, character.only = TRUE))
  }
}

## Auxiliary ##

aya_palette <- c("#E38A42FF", "#42B3E3FF",
                 "#E35DBFFF", "#97C75DFF",
                 "#C75D84FF", "#E6D250FF",
                 "#B84818FF", "#A477C9FF",
                 "#5D9483FF", "#DBD086FF",
                 "#86B2DBFF", "#BFEB9DFF")

## Functions ##

# Show aya color palette and perform colorblind test 
aya_colors <- function() {
  
  print("Arriba y abajo color palette:")
  print(color(aya_palette))
  print("Colorblind checks:")
  print(color(dichromat::dichromat(aya_palette, "protan")))
  print(color(dichromat::dichromat(aya_palette, "deutan")))
  print(color(dichromat::dichromat(aya_palette, "tritan")))
  
}

# Read Arriba output files
read_fusions <- function(working_dir, pattern = "_fusions.tsv$") {
  
  fusion_files <- list.files(working_dir, pattern = pattern, 
                             recursive = TRUE, full.names = TRUE)
  cat("Fusion files found:\n")
  print(fusion_files)
  
  for (i in 1:length(fusion_files)) {
    
    # Read Arriba results
    path2data <- fusion_files[i]
    name <- gsub(".*/|_fusions.tsv","", path2data, perl=TRUE)
    data <- fread(as.character(path2data), header=TRUE)
    colnames(data)[1] <- "gene1"
    
    # Create column with name
    data$Sample <- name
    data$gene1_gene2 <- paste0(data$gene_id1, "_", data$gene_id2)
    
    # Attach to concatenated table
    ifelse(i == 1, fusions.df <- data, fusions.df <- rbind(fusions.df, data))
    
  }
  
  return(fusions.df)
  
}

# Create fusion heatmap
fusion_heatmap <- function(data, metadata = NULL, fusion_mode = "fusion", top_features = 50, custom_order_par = NULL) {
  
  ## metadata: data.frame with rownames = Sample IDs, and columns the metadata parameters to be annotated
  
  # Initial parameter setting
  ifelse(top_features > 50, show_rownames = FALSE, show_rownames = TRUE)
  
  if (fusion_mode == "fusion") {
    
    # Pre-processing
    fusions.tbl <- data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(fusion_mode_par = paste0(gene1, "_", gene2))
    fusions.tbl <- unique(fusions.tbl[,c("Sample", "fusion_mode_par")])
    fusions.tbl <- fusions.tbl %>%
      group_by(Sample, fusion_mode_par) %>%
      dplyr::summarise(num = n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fusion_mode_par) %>%
      dplyr::mutate(tot_num = sum(num)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(tot_num)) %>%
      dplyr::filter(tot_num > 1)
    
  } else if (fusion_mode == "gene") {
    
    # Pre-processing
    fusions.tbl <- data[,c("Sample", "gene1")]
    colnames(fusions.tbl)[2] <- "gene2"
    fusions.tbl <- rbind(fusions.tbl, data[,c("Sample", "gene2")])
    colnames(fusions.tbl)[2] <- "fusion_mode_par"
    fusions.tbl <- unique(fusions.tbl)
    fusions.tbl <- fusions.tbl %>%
      dplyr::filter(fusion_mode_par != ".") %>%
      group_by(Sample, fusion_mode_par) %>%
      dplyr::summarise(num = n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fusion_mode_par) %>%
      dplyr::mutate(tot_num = sum(num)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(tot_num)) %>%
      dplyr::filter(tot_num > 1)
    
  } else {
    
    stop('Error: `fusion_mode` parameter not specfified correctly; choose from: "fusion", "gene"')
    
  }
  
  # Top feature to plot selection
  topfeat <- fusions.tbl[, c("fusion_mode_par", "tot_num")]
  topfeat <- unique(topfeat)
  topfeat <- topfeat$fusion_mode_par[1:top_features]
  
  # Pivot table
  fusions.tbl.piv <- fusions.tbl[fusions.tbl$fusion_mode_par %in% topfeat,] %>%
    dplyr::select(-tot_num) %>%
    tidyr::pivot_wider(names_from = "Sample", values_from = "num") %>%
    dplyr::mutate(across(.cols = everything(), ~replace_na(., 0)))
  fusions.mat <- as.matrix(fusions.tbl.piv[,2:dim(fusions.tbl.piv)[2]])
  rownames(fusions.mat) <- fusions.tbl.piv$fusion_mode_par
  
  # Pre-processing for annotation
  if (!is.null(metadata)) {
    
    # Create annotation color scheme for heatmap
    counter <- 1
    anno_colors <- list()
    for (i in 1:dim(metadata)[2]) {
      j <- length(unique(metadata[,i]))
      
      if (j == 1) {
        assign(paste0("var",i), aya_palette[counter])
      } else {
        if (counter < 12 && (counter+j-1) <= 12) {
          assign(paste0("var",i), aya_palette[counter:(counter+j-1)])
        } else if (counter <= 12 && (counter+j-1) > 12) {
          assign(paste0("var",i), aya_palette[c(counter:12, 1:((counter+j-1) %% 12))])
        }
      }
      
      anno_colors$tmp <- get(paste0("var",i))
      names(anno_colors)[i] <- colnames(metadata)[i]
      
      counter <- (counter + j) %% 12
      
    }
    
    # Custom ordering
    if (!is.null(custom_order_par)) {
      
      custom_ord <- metadata$Sample[sort(metadata[,custom_order_par], index.return = TRUE)$ix]
      custom_ord <- custom_ord[custom_ord %in% colnames(fusions.mat)]
      
      fusions.mat <- fusions.mat[,custom_ord]
      
    }
    
    # Draw heatmap w/ annotation
    return(
      pheatmap(fusions.mat, 
               cluster_rows = TRUE, 
               cluster_cols = FALSE, 
               annotation_col = metadata, 
               annotation_colors = anno_colors,
               show_rownames = show_rownames,
               breaks = c(0,0.1,1), 
               color = c("white", "#a3c771")
      )
    )
    
  } else {
    
    # Draw heatmap w/out annotation
    return(
      pheatmap(fusions.mat, 
               cluster_rows = TRUE, 
               cluster_cols = TRUE,
               show_rownames = show_rownames,
               breaks = c(0,0.1,1), 
               color = c("white", "#a3c771")
      )
    )
    
  }
}

# Perform fusion statistical tests (Fisher's)
fusion_stat_tests <- function(data, metadata, metadata_par, fusion_mode = "fusion", top_features = 100, pval = 0.1) {
  
  if (fusion_mode == "fusion") {
    
    # Pre-processing
    fusions.tbl <- data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(fusion_mode_par = paste0(gene1, "_", gene2))
    fusions.tbl <- unique(fusions.tbl[,c("Sample", "fusion_mode_par")])
    fusions.tbl <- fusions.tbl %>%
      group_by(Sample, fusion_mode_par) %>%
      dplyr::summarise(num = n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fusion_mode_par) %>%
      dplyr::mutate(tot_num = sum(num)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(tot_num)) %>%
      dplyr::filter(tot_num > 1)
    
  } else if (fusion_mode == "gene") {
    
    # Pre-processing
    fusions.tbl <- data[,c("Sample", "gene1")]
    colnames(fusions.tbl)[2] <- "gene2"
    fusions.tbl <- rbind(fusions.tbl, data[,c("Sample", "gene2")])
    colnames(fusions.tbl)[2] <- "fusion_mode_par"
    fusions.tbl <- unique(fusions.tbl)
    fusions.tbl <- fusions.tbl %>%
      dplyr::filter(fusion_mode_par != ".") %>%
      group_by(Sample, fusion_mode_par) %>%
      dplyr::summarise(num = n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fusion_mode_par) %>%
      dplyr::mutate(tot_num = sum(num)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(tot_num)) %>%
      dplyr::filter(tot_num > 1)
    
  } else {
    
    stop('Error: `fusion_mode` parameter not specfified correctly; choose from: "fusion", "gene"')
    
  }
  
  # Pivot table
  fusions.tbl.piv <- fusions.tbl %>%
    dplyr::select(Sample, fusion_mode_par, num) %>%
    tidyr::pivot_wider(names_from = "fusion_mode_par", values_from = "num") %>%
    dplyr::mutate(across(.cols = everything(), ~replace_na(., 0)))
  
  # Attach metadata
  fusions.meta <- left_join(fusions.tbl.piv, metadata, by = "Sample")
  
  for (i in 2:(top_features+1)) {
    
    if (fisher.test(table(unlist(fusions.meta[,i]), unlist(fusions.meta[,metadata_par])))$p.value < pval) {
      
      print(paste0("Found significant result for:", colnames(fusions.meta)[i]))
      print(table(unlist(fusions.meta[,i]), unlist(fusions.meta[,metadata_par])))
      print(fisher.test(table(unlist(fusions.meta[,i]), unlist(fusions.meta[,metadata_par]))))
      cat("\n")
      
    }
  }
}


