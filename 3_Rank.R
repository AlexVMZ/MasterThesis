
# ==============================================================================
# Análisis de los diferentes métodos de normalización en scRNA-seq
# Estudio: HN00234479
# Autor: Alejandro Virues Morales
# ==============================================================================

# ------------------------------------------------------------------------------
# Limpiar el entorno
# ------------------------------------------------------------------------------

rm(list = ls())

# ------------------------------------------------------------------------------
# Cargar librerías necesarias
# ------------------------------------------------------------------------------

library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(paletteer)
library(patchwork)
library(scran)
library(dplyr)
library(RankProd)
library(readr)
library(purrr)
library(ineq)
library(tibble)

################################################################################
##### Cargar objeto SCE y Seurat normalizados para una muestra específica ######
################################################################################

# -------- HeLa --------
sce <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/HeLA_scran_object_normalized.rds")
seurat <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/HeLA_seurat_object_normalized.rds")
file_path <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/HeLA_scran_object_normalized.rds"

# -------- A549 --------
sce <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/A549_scran_object_normalized.rds")
seurat <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/A549_seurat_object_normalized.rds")
file_path <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/A549_scran_object_normalized.rds"

# -------- HUHT --------
sce <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/HUHT_scran_object_normalized.rds")
seurat <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/HUHT_seurat_object_normalized.rds")
file_path <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/HUHT_scran_object_normalized.rds"

# -------- MCF7 --------
sce <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/MCF7_scran_object_normalized.rds")
seurat <- readRDS("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/MCF7_seurat_object_normalized.rds")
file_path <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files/MCF7_scran_object_normalized.rds"

################################################################################
# Definir rutas base y extraer nombre de la muestra
################################################################################

save_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files"
base_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Plots"
sample_name <- sub("_.*$", "", basename(file_path))

sample_dir <- file.path(save_dir, sample_name)

subdirs <- c("RP", "SumRank", "Stat")
for (subdir in subdirs) {
  dir_path <- file.path(sample_dir, subdir)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  if (subdir == "Stat") {
    stat_subdirs <- c("Scran", "SCT")
    for (stat_subdir in stat_subdirs) {
      dir.create(file.path(dir_path, stat_subdir), recursive = TRUE, showWarnings = FALSE)
    }
  }
}

# ------------------------------------------------------------------------------
# Calcular número de células expresando cada gen
# ------------------------------------------------------------------------------

rowData(sce)$nCells <- Matrix::rowSums(counts(sce) > 0)

# ------------------------------------------------------------------------------
# Mostrar resumen estadístico del número de células por gen
# ------------------------------------------------------------------------------

cat({
  resumen <- paste0("Resumen del número de células por gen - ", sample_name, "\n")
  resumen <- paste0(resumen, "--------------------------------------------------------------------\n")
  s <- summary(rowData(sce)$nCells)
  resumen <- paste0(resumen, sprintf("Mínimo     : %6.1f\n", s["Min."]))
  resumen <- paste0(resumen, sprintf("1er cuartil: %6.1f\n", s["1st Qu."]))
  resumen <- paste0(resumen, sprintf("Mediana    : %6.1f\n", s["Median"]))
  resumen <- paste0(resumen, sprintf("Media      : %6.1f\n", s["Mean"]))
  resumen <- paste0(resumen, sprintf("3er cuartil: %6.1f\n", s["3rd Qu."]))
  resumen <- paste0(resumen, sprintf("Máximo     : %6.1f\n", s["Max."]))
  resumen
})


################################################################################
############ Modelado de variabilidad genética con Scran y SCT #################
################################################################################

# Calcular estadísticas de variabilidad con scran
stat_scran <- scran::modelGeneVar(sce)

# Convertir Seurat a SingleCellExperiment y calcular estadísticas para SCT
sce_conv <- Seurat::as.SingleCellExperiment(seurat)
stat_sct <- scran::modelGeneVar(sce_conv)

# ------------------------------------------------------------------------------
# Cálculo de todas las métricas
# ------------------------------------------------------------------------------

# CV² = varianza / media²
stat_scran$CV2_total <- stat_scran$total / (stat_scran$mean)^2
stat_sct$CV2_total <- stat_sct$total / (stat_sct$mean)^2

# Fano Factor = varianza / media
stat_scran$FF_total <- stat_scran$total / stat_scran$mean
stat_sct$FF_total <- stat_sct$total / stat_sct$mean

# Índice de Gini: cuantifica la desigualdad en la expresión génica
expr_mat_scran <- assay(sce, "logcounts")
expr_mat_sct <- assay(sce_conv, "logcounts")

gini_scran <- apply(expr_mat_scran, 1, function(x) ineq(x, type = "Gini"))
gini_sct <- apply(expr_mat_sct, 1, function(x) ineq(x, type = "Gini"))

stat_scran$Gini <- gini_scran
stat_sct$Gini <- gini_sct

# Dropout_Rate: Proporcion de celulas en las que un gen no fue detectado.
dropout_gene_scran <- rowSums(expr_mat_scran == 0) / ncol(expr_mat_scran)
dropout_gene_sct <- rowSums(expr_mat_sct == 0) / ncol(expr_mat_sct)

stat_scran$dropout_rate <- dropout_gene_scran
stat_sct$dropout_rate <- dropout_gene_sct

# IQR/MAD
iqr_scran <- apply(expr_mat_scran, 1, IQR, na.rm = TRUE)
mad_scran <- apply(expr_mat_scran, 1, mad, na.rm = TRUE)
iqr_mad_ratio_scran <- iqr_scran / (mad_scran + 1e-6)

iqr_sct <- apply(expr_mat_sct, 1, IQR, na.rm = TRUE)
mad_sct <- apply(expr_mat_sct, 1, mad, na.rm = TRUE)
iqr_mad_ratio_sct <- iqr_sct / (mad_sct + 1e-6)

stat_scran$IQR_MAD <- iqr_mad_ratio_scran
stat_sct$IQR_MAD <- iqr_mad_ratio_sct

# ------------------------------------------------------------------------------
# Limpieza de datos: reemplazo de NA por 0
# ------------------------------------------------------------------------------

stat_scran[] <- lapply(stat_scran, function(x) {
  if (is.numeric(x)) replace(x, is.na(x), 0) else x
})

stat_sct[] <- lapply(stat_sct, function(x) {
  if (is.numeric(x)) replace(x, is.na(x), 0) else x
})

# ------------------------------------------------------------------------------
# Guardado de métricas como CSV (MetaRanks)
# ------------------------------------------------------------------------------

# -------- Scran --------
df_cv2_scran <- data.frame(
  gene = rownames(stat_scran),
  CV2_total = stat_scran$CV2_total
)

df_ff_scran <- data.frame(
  gene = rownames(stat_scran),
  FF_total = stat_scran$FF_total
)

df_gini_scran <- data.frame(
  gene = rownames(stat_scran),
  Gini = stat_scran$Gini
)

df_dropout_scran <- data.frame(
  gene = rownames(stat_scran),
  dropout_rate = stat_scran$dropout_rate
)

df_iqmad_scran <- data.frame(
  gene = rownames(stat_scran),
  iqrmad = stat_scran$IQR_MAD
)

save_path <- file.path(save_dir, sample_name, paste0("Stat/Scran/Scran_", sample_name, "_CV2.csv"))
readr::write_csv(df_cv2_scran, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/Scran/Scran_", sample_name, "_FF.csv"))
readr::write_csv(df_ff_scran, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/Scran/Scran_", sample_name, "_GI.csv"))
readr::write_csv(df_gini_scran, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/Scran/Scran_", sample_name, "_DR.csv"))
readr::write_csv(df_dropout_scran, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/Scran/Scran_", sample_name, "_IM.csv"))
readr::write_csv(df_iqmad_scran, save_path)

# -------- SCT --------
df_cv2_sct <- data.frame(
  gene = rownames(stat_sct),
  CV2_total = stat_sct$CV2_total
)

df_ff_sct <- data.frame(
  gene = rownames(stat_sct),
  FF_total = stat_sct$FF_total
)

df_gini_sct <- data.frame(
  gene = rownames(stat_sct),
  Gini = stat_sct$Gini
)

df_dropout_sct <- data.frame(
  gene = rownames(stat_sct),
  dropout_rate = stat_sct$dropout_rate
)

df_iqmad_sct <- data.frame(
  gene = rownames(stat_sct),
  iqrmad = stat_sct$IQR_MAD
)

save_path <- file.path(save_dir, sample_name, paste0("Stat/SCT/SCT_", sample_name, "_CV2.csv"))
readr::write_csv(df_cv2_sct, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/SCT/SCT_", sample_name, "_FF.csv"))
readr::write_csv(df_ff_sct, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/SCT/SCT_", sample_name, "_GI.csv"))
readr::write_csv(df_gini_sct, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/SCT/SCT_", sample_name, "_DR.csv"))
readr::write_csv(df_dropout_sct, save_path)

save_path <- file.path(save_dir, sample_name, paste0("Stat/SCT/SCT_", sample_name, "_IM.csv"))
readr::write_csv(df_iqmad_sct, save_path)


stat_scran$tech <- NULL
stat_scran$bio <- NULL
stat_scran$p.value <- NULL
stat_scran$FDR <- NULL

df_out <- as.data.frame(stat_scran)
df_out$gene <- rownames(df_out)
df_out <- df_out[, c("gene", setdiff(names(df_out), "gene"))]
save_path <- file.path(save_dir, sample_name, "Stat", "Scran", paste0("Scran_TotalStat_", sample_name, ".csv"))
readr::write_csv(df_out, save_path)

stat_sct$tech <- NULL
stat_sct$bio <- NULL
stat_sct$p.value <- NULL
stat_sct$FDR <- NULL

df_out <- as.data.frame(stat_sct)
df_out$gene <- rownames(df_out)
df_out <- df_out[, c("gene", setdiff(names(df_out), "gene"))]
save_path <- file.path(save_dir, sample_name, paste0("Stat/SCT/SCT_TotalStat_", sample_name, ".csv"))
readr::write_csv(df_out, save_path)

################################################################################
######################### RANKING DE GENES #####################################
################################################################################

## -- RP -- ##
subsets <- list(
  one_metric = c(3,4),          
  three_metrics = c(3, 4, 5),   
  five_metrics = c(3, 4, 5, 6, 7)
)

run_rank_products <- function(stat_data, subset_cols, tag) {
  gene_names <- rownames(stat_data)
  stat_mat <- as.matrix(as.data.frame(stat_data[, subset_cols]))
  cl <- rep(1, ncol(stat_mat))
  
  RP_out <- RankProducts(stat_mat, cl, gene.names = gene_names, rand = 123)
  
  ranking <- RP_out$RPrank[, "class1 > class2"]
  rps <- RP_out$RPs[, "class1 > class2"]
  pfp <- RP_out$pfp[, "class1 > class2"]
  
  df <- data.frame(
    gene = names(ranking),
    Rank = ranking,
    RP = rps,
    PFP = pfp
  )
  
  df_sorted <- df[order(df$Rank), ]
  
  return(df_sorted)
}

results_scran <- list()
results_sct <- list()

for (name in names(subsets)) {
  subset <- subsets[[name]]
  
  results_scran[[name]] <- run_rank_products(stat_scran, subset, paste0("scran_", name))
  results_sct[[name]] <- run_rank_products(stat_sct, subset, paste0("sct_", name))
}

base_path <- file.path(save_dir, sample_name, "RP")

saveRDS(results_scran, file = file.path(base_path, paste0("RP_Rank_scran_",sample_name,".rds")))
saveRDS(results_sct, file = file.path(base_path, paste0("RP_Rank_sct_",sample_name,".rds")))


## -- SumRank -- ##

# Definir subconjuntos de métricas
metrics_list <- list(
  rank_1 = c("CV2_total"),
  rank_3 = c("CV2_total", "FF_total", "Gini"),
  rank_5 = c("CV2_total", "FF_total", "Gini", "dropout_rate", "IQR_MAD")
)

# Función para calcular ranking dada una tabla y un conjunto de métricas
compute_ranking <- function(df, metrics) {
  df %>%
    mutate(across(all_of(metrics), ~rank(-., ties.method = "average"), .names = "rank_{.col}")) %>%
    rowwise() %>%
    mutate(total_rank_score = sum(c_across(starts_with("rank_")), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(global_rank = rank(total_rank_score, ties.method = "min")) %>%
    arrange(global_rank)
}

# Preparar dataframes
stat_df_scran <- as.data.frame(stat_scran) %>%
  rownames_to_column(var = "gene")

stat_df_sct <- as.data.frame(stat_sct) %>%
  rownames_to_column(var = "gene")

# Aplicar a scran
results_rank_scran <- lapply(metrics_list, function(m) compute_ranking(stat_df_scran, m))

# Aplicar a seurat
results_rank_sct <- lapply(metrics_list, function(m) compute_ranking(stat_df_sct, m))

# Asignar nombres a las listas
names(results_rank_scran) <- names(metrics_list)
names(results_rank_sct) <- names(metrics_list)

base_path <- file.path(save_dir, sample_name, "SumRank")

saveRDS(results_rank_scran, file = file.path(base_path, paste0("SumRank_scran_",sample_name,".rds")))
saveRDS(results_rank_sct, file = file.path(base_path, paste0("SumRank_sct_",sample_name,".rds")))
