library(tidyverse)
library(pheatmap)
library(purrr) # importante para reduce()

base_input_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files"
save_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/4_Metrics_Analyses"

cell_lines <- c("HeLA", "A549", "MCF7", "HUHT")
methods <- c("scran", "sct")
rank_types <- c("SumRank", "RP")

# Función para extraer dataframes con las columnas de rankings, manejando formatos distintos
get_ranking_columns <- function(obj, rank_type) {
  
  if (rank_type == "SumRank") {
    # En SumRank, los elementos son dataframes con columnas tipo rank_1, rank_3, rank_5
    rank_names <- c("rank_1", "rank_3", "rank_5")
    # Filtrar para quedarnos con los objetos que existan y tengan esa columna
    dfs <- map(rank_names, function(name) {
      if (name %in% names(obj)) {
        # Es un tibble completo, extraemos solo columna gene y columna de ranking global
        df <- obj[[name]]
        # Buscar la columna que contenga 'global_rank' o 'rank' para ranking
        col_rank <- intersect(c("global_rank", "rank_CV2_total", "rank"), names(df))
        if (length(col_rank) == 0) {
          # si no hay ranking, intentar con la primera columna numérica después de gene
          col_rank <- names(df)[which(names(df) != "gene")[1]]
        }
        df %>%
          select(gene, !!sym(col_rank[1])) %>%
          rename(rank = !!sym(col_rank[1]))
      } else {
        NULL
      }
    }) %>% compact()
    
    if(length(dfs) < 2) return(NULL)
    
    # Merge por 'gene'
    df_merged <- purrr::reduce(dfs, full_join, by = "gene")
    colnames(df_merged)[-1] <- rank_names[1:length(dfs)]
    
  } else if (rank_type == "RP") {
    # En RP, los elementos one_metric, three_metrics, five_metrics son dataframes con 'Rank'
    rank_names <- c("one_metric", "three_metrics", "five_metrics")
    dfs <- map(rank_names, function(name) {
      if (name %in% names(obj)) {
        df <- obj[[name]]
        if(!"Rank" %in% names(df)) {
          stop(paste("No se encontró columna 'Rank' en", name))
        }
        df %>%
          select(gene, Rank) %>%
          rename(rank = Rank)
      } else {
        NULL
      }
    }) %>% compact()
    
    if(length(dfs) < 2) return(NULL)
    
    df_merged <- purrr::reduce(dfs, full_join, by = "gene")
    colnames(df_merged)[-1] <- rank_names[1:length(dfs)]
  } else {
    return(NULL)
  }
  
  return(df_merged)
}

calculate_spearman_matrix <- function(df) {
  mat <- as.matrix(df[,-1])
  cor(mat, method = "spearman", use = "pairwise.complete.obs")
}

for (cell in cell_lines) {
  for (method in methods) {
    for (rank_type in rank_types) {
      
      file_name <- if (rank_type == "RP") {
        paste0("RP_Rank_", method, "_", cell, ".rds")
      } else {
        paste0("SumRank_", method, "_", cell, ".rds")
      }
      
      input_path <- file.path(base_input_dir, cell, rank_type, file_name)
      
      if (!file.exists(input_path)) {
        message("Archivo no encontrado: ", input_path)
        next
      }
      
      rank_obj <- readRDS(input_path)
      rank_df <- get_ranking_columns(rank_obj, rank_type)
      
      if (is.null(rank_df)) {
        message("No hay suficientes columnas para calcular correlación en: ", input_path)
        next
      }
      
      spearman_matrix <- calculate_spearman_matrix(rank_df)
      
      sample_dir <- file.path(save_dir, cell)
      dir.create(file.path(sample_dir, "Files"), recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(sample_dir, "Plots"), recursive = TRUE, showWarnings = FALSE)
      
      csv_path <- file.path(sample_dir, "Files", paste0("spearman_", method, "_", rank_type, ".csv"))
      write.csv(round(spearman_matrix, 4), csv_path, row.names = TRUE)
      
      plot_path <- file.path(sample_dir, "Plots", paste0("heatmap_", method, "_", rank_type, ".png"))
      png(plot_path, width = 1000, height = 800)
      pheatmap(spearman_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
               main = paste("Spearman Correlation -", toupper(rank_type), "-", method),
               fontsize = 14,
               color = colorRampPalette(c("#deebf7","#08306B"))(50))
      dev.off()
      
      message("Procesado: ", cell, " - ", method, " - ", rank_type)
    }
  }
}