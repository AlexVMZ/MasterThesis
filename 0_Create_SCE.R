# ==============================================================================
# Creación de objetos SingleCellExperiment para múltiples muestras
# Estudio: HN00234479
# Autor: Alejandro Virues Morales
# ==============================================================================

# Limpiar el entorno
rm(list = ls())

# Cargar librerías necesarias
library(SingleCellExperiment)
library(DropletUtils)
library(dplyr)
library(BiocParallel)

# Definir rutas de trabajo
base_dir <- "/home/avirues/TFM/Data/HN00234479/Study/0_CreateSCE/"
results_dir <- file.path(base_dir, "Results")

# Lista de muestras con sus rutas y directorios de matrices
samples <- list(
  A549  = list(path = "Data/Samples/A549/H5K72DSXF/SCE_Data", matrix = "A549_raw_feature_bc_matrix"),
  HeLA  = list(path = "Data/Samples/HeLA/H5K72DSXF/SCE_Data", matrix = "HeLA_raw_feature_bc_matrix"),
  HUHT  = list(path = "Data/Samples/HUHT/H5K72DSXF/SCE_Data", matrix = "HUHT_raw_feature_bc_matrix"),
  MCF7  = list(path = "Data/Samples/MCF7/H5K72DSXF/SCE_Data", matrix = "MCF7_raw_feature_bc_matrix")
)

# Lista donde se guardarán los objetos SingleCellExperiment
sce_list <- list()

# ==============================================================================
# Función para procesar cada muestra
# ==============================================================================
procesar_muestra <- function(sample_name, sample_info) {
  
  cat("\n========================\n")
  cat("Procesando muestra:", sample_name, "\n")
  cat("========================\n")
  
  # Establecer el directorio de trabajo
  setwd(file.path(base_dir, sample_info$path))
  print(paste("Directorio de trabajo:", getwd()))
  
  # Leer datos de 10X Genomics
  sce <- read10xCounts(sample_info$matrix, col.names = FALSE, row.names = "symbol")
  rownames(colData(sce)) <- sce$Barcode
  
  # Renombrar metadatos de genes
  rowData(sce) <- rowData(sce)[, -3]
  colnames(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
  
  # Agregar metadatos de la muestra
  sce$Sample <- sample_name
  
  # ==================================================================
  # Evaluación de barcodes
  # ==================================================================
  br.out <- barcodeRanks(counts(sce))
  
  # Ordenar datos para trazado
  o <- order(br.out$rank)
  rank_sorted <- br.out$rank[o]
  total_sorted <- br.out$total[o]
  
  # Obtener umbrales
  knee_value <- metadata(br.out)$knee
  inflection_value <- metadata(br.out)$inflection
  threshold_value <- 100
  
  # Calcular número de células por encima de cada umbral
  cells_above_knee <- sum(total_sorted >= knee_value)
  cells_above_inflection <- sum(total_sorted >= inflection_value)
  cells_above_threshold <- sum(total_sorted >= threshold_value)
  
  # Guardar gráfico de barcode ranks
  jpeg(file = file.path(results_dir, paste0(sample_name, "_BCRanks.jpeg")),
       width = 8, height = 6, units = "in", res = 300)
  plot(br.out$rank, br.out$total, log = "xy", xlab = "Rank", ylab = "Total")
  lines(rank_sorted, total_sorted, col = "red")
  abline(h = c(knee_value, inflection_value, threshold_value), 
         col = c("dodgerblue", "forestgreen", "red"), lty = 2)
  
  # Identificar puntos de corte en la curva
  knee_rank <- rank_sorted[which.min(abs(total_sorted - knee_value))]
  inflection_rank <- rank_sorted[which.min(abs(total_sorted - inflection_value))]
  threshold_rank <- rank_sorted[which.min(abs(total_sorted - threshold_value))]
  
  # Etiquetas en el gráfico con número de células en lugar del umbral
  text(knee_rank, knee_value, labels = paste0("Knee: ", cells_above_knee, " cells"), pos = 4, col = "dodgerblue")
  text(inflection_rank, inflection_value, labels = paste0("Inflection: ", cells_above_inflection, " cells"), pos = 4, col = "forestgreen")
  text(threshold_rank, threshold_value, labels = paste0("Threshold: ", cells_above_threshold, " cells"), pos = 4, col = "red")
  
  # Leyenda con el número de células
  legend("bottomleft", lty = 2, col = c("dodgerblue", "forestgreen", "red"),
         legend = c(paste0("Knee: ", round(knee_value, 2)), 
                    paste0("Inflection: ", round(inflection_value, 2)), 
                    paste0("Threshold: ", threshold_value)))
  dev.off()
  
  # ==================================================================
  # Eliminación de gotas vacías
  # ==================================================================
  print("Dimensiones antes de eliminar gotas vacías:")
  print(dim(sce))

  set.seed(1234)
  e.out <- emptyDrops(counts(sce))

  # Evaluar filtrado
  print(summary(e.out$LogProb))
  print(summary(e.out$FDR))
  is.cell <- e.out$FDR <= 0.001
  table(Limited = e.out$Limited, Significant = is.cell)

  # Aplicar filtrado
  sce <- sce[, which(e.out$FDR <= 0.001)]

  print("Dimensiones tras eliminar gotas vacías:")
  print(dim(sce))

  # Guardar gráfico de distribución de gotas vacías
  jpeg(file = file.path(results_dir, paste0("Empty_Cells_", sample_name, ".jpeg")),
       width = 8, height = 6, units = "in", res = 300)
  plot(e.out$Total, -e.out$LogProb, col = ifelse(is.cell, "red", "black"),
       xlab = "Total UMI count", ylab = "-Log Probability")
  dev.off()
  
  # ==================================================================
  # Guardar objeto SingleCellExperiment individual
  # ==================================================================
  saveRDS(sce, file.path(results_dir, paste0("sce_", sample_name, "_raw_emptyDrops.rds")))
  
  print(paste("Objeto SingleCellExperiment guardado con éxito:", sample_name))
  print(paste("Dimensiones finales:", dim(sce)[1], "features *", dim(sce)[2], "cells"))
  
  return(sce)
}

# ==================================================================
# Procesar todas las muestras y almacenarlas en una lista
# ==================================================================
sce_list <- lapply(names(samples), function(sample_name) {
  procesar_muestra(sample_name, samples[[sample_name]])
})

print("Proceso finalizado para todas las muestras.")