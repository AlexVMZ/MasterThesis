# ==============================================================================
# Normalización de muestras scRNA-seq con Seurat y scran
# Estudio: HN00234479
# Autor: Alejandro Virues Morales
# ==============================================================================

# Limpiar entorno
rm(list = ls())

# ------------------------------------------------------------------------------
# Cargar librerías necesarias
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(sctransform)
  library(scran)
  library(scater)
})

# ------------------------------------------------------------------------------
# Configuración general
# ------------------------------------------------------------------------------
input_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/1_QC/Files"
output_data_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files"
output_plot_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Plots"

samples <- c("HeLA", "A549", "HUHT", "MCF7")

# ------------------------------------------------------------------------------
# Función para normalizar con Seurat
# ------------------------------------------------------------------------------
normalize_with_seurat <- function(sce, sample_name, output_dir) {
  message("[Seurat] Normalizando muestra: ", sample_name)
  seurat_obj <- Seurat::as.Seurat(sce, data = NULL)

  seurat_obj <- Seurat::SCTransform(seurat_obj, assay = "originalexp", new.assay.name = "SCT",
                            min_cells = 0, do.correct.umi = TRUE, ncells = 5000,
                            variable.features.n = 3000, variable.features.rv.th = 1.3,
                            do.scale = FALSE, do.center = TRUE, vst.flavor = "v2",
                            conserve.memory = FALSE, return.only.var.genes = FALSE,
                            seed.use = 1448145, verbose = TRUE)
  
  saveRDS(seurat_obj, file.path(output_dir, paste0(sample_name, "_seurat_object_normalized.rds")))
}

# ------------------------------------------------------------------------------
# Función para normalizar con scran
# ------------------------------------------------------------------------------
normalize_with_scran <- function(sce, sample_name, output_dir) {
  message("[scran] Normalizando muestra: ", sample_name)
  clusters <- scran::quickCluster(sce)
  sce <- scran::computeSumFactors(sce, clusters = clusters)
  sce <- scater::logNormCounts(sce)
  saveRDS(sce, file.path(output_dir, paste0(sample_name, "_scran_object_normalized.rds")))
}

# ------------------------------------------------------------------------------
# Proceso principal por muestra
# ------------------------------------------------------------------------------
for (sample in samples) {
  message("========== Procesando muestra: ", sample, " ==========")
  
  file_path <- file.path(input_dir, paste0(sample, "_emptyDrops_processed.rds"))
  if (!file.exists(file_path)) {
    warning("Archivo no encontrado: ", file_path)
    next
  }
  
  # Leer objeto
  sce <- readRDS(file_path)
  
  # Ajustar nombres si es necesario
  rownames(sce) <- gsub("_", "-", rownames(sce))
  
  # Crear carpeta de salida para gráficos (si se usa más adelante)
  dir.create(file.path(output_plot_dir, sample), recursive = TRUE, showWarnings = FALSE)
  
  # Normalización
  normalize_with_seurat(sce, sample, output_data_dir)
  normalize_with_scran(sce, sample, output_data_dir)
}

message("✅ Proceso de normalización completado.")

