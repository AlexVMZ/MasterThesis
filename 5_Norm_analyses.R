# ==============================================================================
# Análisis de los diferentes métodos de normalización en scRNA-seq
# Estudio: HN00234479
# Autor: Alejandro Virues Morales
# ==============================================================================

# Limpiar el entorno
rm(list = ls())

# Cargar librerías necesarias
# Librerías esenciales
library(SingleCellExperiment) 
library(Seurat)             
library(scran)      
library(scater)     
library(ggplot2)           
library(patchwork)          
library(dplyr)               
library(tibble)             
library(clusterProfiler)    
library(org.Hs.eg.db)       
library(enrichplot)          
library(ggpubr) 
library(grid)   

sample_name <- "HeLA"  # Cambia al nombre que se quiera estudiar ("HeLA", "A549", "MCF7", "HUHT)

base_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files"
save_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/5_Norm_Analyses"
sample_dir <- file.path(save_dir, sample_name)

# Definir rutas separadas para Files y Plots
save_dirFiles <- file.path(save_dir, sample_name, "Files")
save_dirPlot <- file.path(save_dir, sample_name, "Plots")

# Crear los directorios
dir.create(save_dirFiles, recursive = TRUE, showWarnings = FALSE)
dir.create(save_dirPlot, recursive = TRUE, showWarnings = FALSE)

rp_dir <- file.path(base_dir, sample_name, "RP")

# Cargar RP_Rank_sct
file_sct <- file.path(rp_dir, paste0("RP_Rank_sct_", sample_name, ".rds"))
RP_sct <- readRDS(file_sct)
rank_sct <- RP_sct$three_metrics

# Cargar RP_Rank_scran
file_scran <- file.path(rp_dir, paste0("RP_Rank_scran_", sample_name, ".rds"))
RP_scran <- readRDS(file_scran)
rank_scran <- RP_scran$three_metrics


stat_dir <- file.path(base_dir, sample_name, "Stat")

# Ruta para SCT
file_sct_stat <- file.path(stat_dir, "SCT", paste0("SCT_TotalStat_", sample_name, ".csv"))
stat_sct <- read.csv(file_sct_stat, stringsAsFactors = FALSE)

# Ruta para Scran
file_scran_stat <- file.path(stat_dir, "Scran", paste0("Scran_TotalStat_", sample_name, ".csv"))
stat_scran <- read.csv(file_scran_stat, stringsAsFactors = FALSE)


################################################################################
################## Comparativa de métricas: Seurat SCT vs Scran ################
################################################################################

# Combina los dataframes con estadísticas de cada método
df_stats_sct <- data.frame(stat_sct, method = "SCT")
df_stats_scran <- data.frame(stat_scran, method = "Scran")
df_comparison <- rbind(df_stats_sct, df_stats_scran)

# ------------------------------------------------------------------------------
# Gráficos de comparación de métricas
# ------------------------------------------------------------------------------

# Media de expresión
p1 <- ggplot(df_comparison, aes(x = method, y = mean, fill = method)) +
  geom_violin() +
  labs(y = "Media (Mean)", title = paste0("Medias Scran/Seurat - ", sample_name)) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme_minimal()

# Varianza total
p2 <- ggplot(df_comparison, aes(x = method, y = total, fill = method)) +
  geom_violin() +
  labs(y = "Coeficiente de Varianza total", title = paste0("Varianzas Scran/Seurat - ", sample_name)) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme_minimal() +
  ylim(0, 1)

# Coeficiente de Variación al cuadrado (CV2)
p3 <- ggplot(df_comparison, aes(x = method, y = CV2_total, fill = method)) +
  geom_violin() +
  labs(y = "Coeficiente de Variacion al cuadrado (CV2)", title = paste0("CV2 Scran/Seurat - ", sample_name)) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme_minimal() +
  ylim(0, 30)

# Factor de Fano (FF)
p4 <- ggplot(df_comparison, aes(x = method, y = FF_total, fill = method)) +
  geom_violin() +
  labs(y = "Factor Fano (FF)", title = paste0("Factor Fano Scran/Seurat - ", sample_name)) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme_minimal()

# Índice de Gini
p5 <- ggplot(df_comparison, aes(x = method, y = Gini, fill = method)) +
  geom_violin() +
  labs(y = "Gini Index", title = paste0("Gini Index Scran/Seurat - ", sample_name)) +
  scale_fill_manual(values = c("skyblue", "orange")) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Combinación de los gráficos con patchwork
# ------------------------------------------------------------------------------

combined_plot_1 <- (p1 + p2)
combined_plot_2 <- (p3 + p4) / (p5)

# Guardado de los gráficos combinados
ggsave(filename = file.path(save_dirPlot, paste0("Metrics_combined_plot_1_", sample_name, ".png")),
       plot = combined_plot_1, width = 10, height = 8, dpi = 100, limitsize = FALSE)

ggsave(filename = file.path(save_dirPlot, paste0("Metrics_combined_plot_2_", sample_name, ".png")),
       plot = combined_plot_2, width = 10, height = 8, dpi = 100, limitsize = FALSE)

# ------------------------------------------------------------------------------
# Comparativa de la distribución de la matriz de expresión
# ------------------------------------------------------------------------------

base_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/2_Normalization/Files"

path_sce <- file.path(base_dir, paste0(sample_name, "_scran_object_normalized.rds"))
path_seurat <- file.path(base_dir, paste0(sample_name, "_seurat_object_normalized.rds"))

sce <- readRDS(path_sce)
seurat <- readRDS(path_seurat)

# Función para datos crudos (raw)
plot_density_log <- function(counts_mat, label) {
  expr_values <- as.numeric(counts_mat)
  expr_values <- expr_values[expr_values > 0]
  df <- data.frame(value = log2(expr_values + 1), method = label)
  return(df)
}

# Función para datos normalizados
plot_density_norm_matrix <- function(counts_mat, label) {
  expr_values <- as.numeric(counts_mat)
  df <- data.frame(value = expr_values, method = label)
  return(df)
}

# Generación de dataframes
df_raw <- plot_density_log(assay(sce, "counts"), "Raw")
df_scran <- plot_density_norm_matrix(assay(sce, "logcounts"), "Scran")
df_sct <- plot_density_norm_matrix(seurat@assays$SCT@data, "Seurat")
all_df <- rbind(df_raw, df_scran, df_sct)

# Guardado del gráfico de densidad
save_path <- file.path(save_dirPlot, paste0("Density_Plot_Matrix_Expresion_", sample_name, ".png"))
png(filename = save_path, width = 1000, height = 800)
expresion_counts_matrix <- ggplot(all_df, aes(x = value, color = method, fill = method)) +
  geom_density() +
  labs(x = "log2(Expression + 1)", y = "Densidad", title = paste0("Distribución de la matriz de expresión ", sample_name)) +
  coord_cartesian(xlim = c(0.15, 3)) +
  theme_minimal()
print(expresion_counts_matrix)
dev.off()

# ------------------------------------------------------------------------------
# Comparación Media vs Varianza (método log y lineal)
# ------------------------------------------------------------------------------

# Comparación lineal
save_path2 <- file.path(save_dirPlot, paste0("Comparision_mean_var_", sample_name, ".png"))
w <- ggplot() +
  geom_point(data = stat_scran, aes(x = mean, y = total), color = "skyblue") +
  geom_point(data = stat_sct, aes(x = mean, y = total), color = "orange") +
  labs(x = "Media de expresión", y = "Varianza",
       title = paste0("Comparación entre Seurat y Scran (Media vs Varianza) - ", sample_name)) +
  theme_minimal()
ggsave(filename = save_path2, plot = w, width = 10, height = 8, dpi = 300)

q <- ggplot(df_comparison, aes(x = mean, y = total, color = method)) +
  geom_point() +
  scale_color_manual(values = c( "SCT" = "orange","Scran" = "skyblue")) +
  labs(x = "Media de expresión", y = "Varianza",
       title = paste0("Comparación entre Seurat y Scran (Media vs Varianza) - ", sample_name),
       color = "Método") +
  theme_minimal()

################################################################################
##################### Cambios entre ranking ####################################
################################################################################

comparison_df_metrics <- dplyr::inner_join(
  dplyr::select(rank_scran, gene, Rank_Scran = Rank),
  dplyr::select(rank_sct, gene, Rank_SCT = Rank),
  by = "gene"
)

comparison_df_metrics <- comparison_df_metrics %>%
  dplyr::mutate(rank_shift = Rank_Scran - Rank_SCT) %>%
  dplyr::arrange(desc(abs(rank_shift)))

genes_shift <- comparison_df_metrics$gene[comparison_df_metrics$rank_shift < -5000]

comparison_df_metrics$color_group <- ifelse(comparison_df_metrics$gene %in% genes_shift, "highlight", "other")

### Guardado comparaciones ###
save_path_comp <- file.path(save_dirFiles, paste0("Cambio_ranking_metrics_", sample_name, ".csv"))
readr::write_csv(comparison_df_metrics, save_path_comp)

################################################################################
# Preparación de nombres y datos para visualización
################################################################################
df_scran <- rank_scran %>%
  left_join(stat_scran %>% select(gene, mean, total, CV2_total), by = "gene")

df_sct <- rank_sct %>%
  left_join(stat_sct %>% select(gene, mean, total, CV2_total), by = "gene")


df_sct_sub <- df_sct %>%
  select(gene, rank_sct = Rank, mean_sct = mean, CV2_sct = CV2_total)

df_scran_sub <- df_scran %>%
  select(gene, rank_scran = Rank, mean_scran = mean, CV2_scran = CV2_total)

df_combined <- df_sct_sub %>%
  inner_join(df_scran_sub, by = "gene")

nCells_df <- data.frame(
  gene = rownames(rowData(sce)),
  nCells = rowData(sce)$nCells
)

df_combined <- df_combined %>%
  left_join(nCells_df, by = "gene")

df_sct <- df_sct %>%
  left_join(nCells_df, by = "gene")

df_scran <- df_scran %>%
  left_join(nCells_df, by = "gene")

df_combined$color_group <- ifelse(df_combined$gene %in% genes_shift, "highlight", "other")
df_scran$color_group <- ifelse(df_scran$gene %in% genes_shift, "highlight", "other")
df_sct$color_group <- ifelse(df_sct$gene %in% genes_shift, "highlight", "other")

################################################################################
## VISUALIZACIÓN DE INTER-RANKINGS ENTRE MÉTODOS DE NORMALIZACIÓN
################################################################################

# Directorio de salida para los gráficos
save_path <- file.path(save_dirPlot, paste0("Comparision_global_rank_metrics_", sample_name, ".png"))

# Comparación directa de rankings globales entre Scran y Sctransform
png(filename = save_path, width = 1000, height = 800)
plot_comparacion <- ggplot(df_combined, aes(x = rank_scran, y = rank_sct)) +
  geom_point(aes(color = color_group), alpha = 0.5) +
  scale_color_manual(values = c("highlight" = "red", "other" = "grey")) +
  geom_abline(color = "red", linetype = "dashed") +
  labs(
    title = paste0("Comparación de Rankings Globales - ", sample_name),
    x = "Ranking global (Scran)",
    y = "Ranking global (Sctransform)"
  ) +
  theme_minimal(base_size = 25)
print(plot_comparacion)
dev.off()

# Unir por gen
merged_ranks <- merge(rank_scran[, c("gene", "Rank")],
                      rank_sct[, c("gene", "Rank")],
                      by = "gene",
                      suffixes = c("_scran", "_sct"))
correlation_spearman <- cor(merged_ranks$Rank_scran, merged_ranks$Rank_sct, method = "spearman")
print(correlation_spearman)

# Comparación entre el ranking y la media de expresión para Scran
save_path <- file.path(save_dirPlot, paste0("Comparision_mean_rank_scran_", sample_name, ".png"))
png(filename = save_path, width = 1000, height = 800)
plot_scran <- ggplot(df_scran, aes(x = Rank, y = mean)) +
  geom_point(aes(color = color_group), alpha = 0.5) +
  scale_color_manual(values = c("highlight" = "red", "other" = "grey")) +
  labs(
    title = paste0("Comparación de Rankings Media - Rank - ", sample_name),
    x = "Ranking Scran",
    y = "Media Expresión"
  ) +
  theme_minimal(base_size = 25)
print(plot_scran)
dev.off()

# Comparación entre el ranking y la media de expresión para Sctransform
save_path <- file.path(save_dirPlot, paste0("Comparision_mean_rank_sct_", sample_name, ".png"))
png(filename = save_path, width = 1000, height = 800)
plot_sct <- ggplot(df_sct, aes(x = Rank, y = mean)) +
  geom_point(aes(color = color_group), alpha = 0.5) +
  scale_color_manual(values = c("highlight" = "red", "other" = "grey")) +
  labs(
    title = paste0("Comparación de Rankings Media - Rank - ", sample_name),
    x = "Ranking SCT",
    y = "Media Expresión"
  ) +
  theme_minimal()
print(plot_sct)
dev.off()


################################################################################
## ANÁLISIS FUNCIONAL (ORA) DE RANKINGS SEGÚN MÉTODO DE NORMALIZACIÓN
################################################################################

### ORA de los genes que no rankean por SCT (genes únicos de SHIFT) ####
gene_ids_shift <- clusterProfiler::bitr(genes_shift, 
                                        fromType = "SYMBOL", 
                                        toType = "ENTREZID", 
                                        OrgDb = org.Hs.eg.db)

ego_shift <- clusterProfiler::enrichGO(
  gene = gene_ids_shift$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.25,
  readable = TRUE
)

n_categories <- length(ego_shift@result$Description)
max_categories <- 10
categories_to_show <- min(n_categories, max_categories)

save_path <- file.path(save_dirPlot, paste0("ORA_Shift_genes_", sample_name, ".png"))
png(filename = save_path, width = 1000, height = 800)
barplot(ego_shift, showCategory = categories_to_show, 
        title = paste0("ORA - GO: Biological Process - Shift Gene - ", sample_name))
dev.off()


### ORA de los top 1000 genes según ranking por método ####

#### Método SCT ####
topGeneSCT <- df_sct$gene[df_sct$Rank < 1001]

top1000_SCT_gene_ids <- clusterProfiler::bitr(topGeneSCT, 
                                              fromType = "SYMBOL", 
                                              toType = "ENTREZID", 
                                              OrgDb = org.Hs.eg.db)

ego_top1000_SCT <- clusterProfiler::enrichGO(
  gene = top1000_SCT_gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.25,
  readable = TRUE
)

n_categories <- length(ego_top1000_SCT@result$Description)
categories_to_show <- min(n_categories, max_categories)

save_path <- file.path(save_dirPlot, paste0("ORA_Top1000_genes_SCT_", sample_name, ".png"))
png(filename = save_path, width = 1000, height = 800)
barplot(ego_top1000_SCT, showCategory = categories_to_show, 
        title = paste0("ORA - GO: Biological Process - Top 1000 SCT - ", sample_name))
dev.off()


#### Método Scran ####
topGeneSran <- df_scran$gene[df_scran$Rank < 1001]

top1000_Scran_gene_ids <- clusterProfiler::bitr(topGeneSran, 
                                                fromType = "SYMBOL", 
                                                toType = "ENTREZID", 
                                                OrgDb = org.Hs.eg.db)

ego_top1000_Scran <- clusterProfiler::enrichGO(
  gene = top1000_Scran_gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.25,
  readable = TRUE
)

n_categories <- length(ego_top1000_Scran@result$Description)
categories_to_show <- min(n_categories, max_categories)

save_path <- file.path(save_dirPlot, paste0("ORA_Top1000_genes_Scran_", sample_name, ".png"))
png(filename = save_path, width = 1000, height = 800)
barplot(ego_top1000_Scran, showCategory = categories_to_show, 
        title = paste0("ORA - GO: Biological Process - Top 1000 Scran - ", sample_name))
dev.off()



################################################################################
## STATS GENES - Análisis de métricas y clasificaciones de genes
################################################################################

##### --- FUNCIONES AUXILIARES

# Función: Crear tabla resumen de tipo (clasificación por conjuntos)
crear_resumen_tipo <- function(lista_sct, lista_scran, lista_shift) {
  all_genes <- unique(c(lista_sct, lista_scran, lista_shift))
  
  tibble(gene = all_genes) %>%
    dplyr::mutate(
      in_SCT = gene %in% lista_sct,
      in_Scran = gene %in% lista_scran,
      in_Shift = gene %in% lista_shift
    ) %>%
    rowwise() %>%
    dplyr::mutate(
      type = case_when(
        in_SCT & !in_Scran & !in_Shift ~ "Only_SCT",
        !in_SCT & in_Scran & !in_Shift ~ "Only_Scran",
        !in_SCT & !in_Scran & in_Shift ~ "Only_Shift",
        in_SCT & in_Scran & !in_Shift ~ "SCT&Scran",
        in_SCT & !in_Scran & in_Shift ~ "SCT&Shift",
        !in_SCT & in_Scran & in_Shift ~ "Scran&Shift",
        in_SCT & in_Scran & in_Shift ~ "SCT&Scran&Shift"
      )
    ) %>%
    ungroup()
}

# Función: Extraer tabla resumen de métricas para una lista de genes
crear_df_genes <- function(lista_genes, tipo, df_sct, df_scran, row_data) {
  lista_genes <- as.character(lista_genes)
  
  sct_info <- df_sct %>%
    dplyr::filter(gene %in% lista_genes) %>%
    dplyr::select(gene, mean_sct = mean, rank_sct = Rank)
  
  scran_info <- df_scran %>%
    dplyr::filter(gene %in% lista_genes) %>%
    dplyr::select(gene, mean_scran = mean, rank_scran = Rank)
  
  nCells_info <- as.data.frame(row_data) %>%
    dplyr::filter(SYMBOL %in% lista_genes) %>%
    dplyr::select(gene = SYMBOL, nCells)
  
  tibble(gene = lista_genes) %>%
    dplyr::left_join(sct_info, by = "gene") %>%
    dplyr::left_join(scran_info, by = "gene") %>%
    dplyr::left_join(nCells_info, by = "gene") %>%
    dplyr::mutate(type = tipo)
}


##### --- PREPARAR TABLAS DE MÉTRICAS


# Métricas de SCT
sct_info <- df_sct %>%
  dplyr::select(
    gene,
    mean_sct = mean,
    rank_sct = Rank,
    CV2_sct = CV2_total
  )

# Métricas de Scran
scran_info <- df_scran %>%
  dplyr::select(
    gene,
    mean_scran = mean,
    rank_scran = Rank,
    CV2_scran = CV2_total
  )

# Número de células en las que se detecta cada gen
nCells_info <- as.data.frame(rowData(sce)) %>%
  dplyr::select(gene = SYMBOL, nCells)


###### --- 3. CREAR TABLA DE RESUMEN (TOP GENES)


resumen_tipo <- crear_resumen_tipo(topGeneSCT, topGeneSran, genes_shift)

df_resumen <- resumen_tipo %>%
  dplyr::left_join(df_combined, by = "gene")


###### --- CREAR RESÚMENES INDIVIDUALES PARA CADA LISTA

df_top_sct <- crear_df_genes(topGeneSCT, "Top1000SCT", df_sct, df_scran, rowData(sce))
df_top_scran <- crear_df_genes(topGeneSran, "Top1000Scran", df_sct, df_scran, rowData(sce))
df_shift <- crear_df_genes(genes_shift, "Shift", df_sct, df_scran, rowData(sce))

df_resumen_2 <- dplyr::bind_rows(df_top_sct, df_top_scran, df_shift)

###### --- EXPORTAR RESULTADOS
# Guardar archivos
readr::write_csv(df_resumen, file.path(save_dirFiles, paste0("Resumen_TOP_Shift_Stats_", sample_name, ".csv")))
readr::write_csv(df_combined, file.path(save_dirFiles, paste0("Stats_TODOS_Genes_SCT_Scran_Cells_", sample_name, ".csv")))
readr::write_csv(df_resumen_2, file.path(save_dirFiles, paste0("Resumen_TOP_Shift_Stats_Simple_", sample_name, ".csv")))
readr::write_csv(df_sct, file.path(save_dirFiles, paste0("Stat_Charateristic_SCT_", sample_name, ".csv")))
readr::write_csv(df_scran, file.path(save_dirFiles, paste0("Stat_Charateristic_Scran_", sample_name, ".csv")))

################################################################################
# Visualización STAT - Frec, nCells y rankings
################################################################################

######### --- Frecuencia de genes por tipo de pertenencia --- ####
df_frecuencia <- df_resumen %>%
  dplyr::count(type) %>%
  dplyr::arrange(desc(n))

save_path_freq <- file.path(save_dirPlot, paste0("Frequencies_By_Type_", sample_name, ".png"))
png(filename = save_path_freq, width = 1000, height = 800)
ggplot(df_frecuencia, aes(x = reorder(type, n), y = n, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n), hjust = -0.1, size = 5) +  # Tamaño del número sobre la barra
  coord_flip() +
  labs(
    title = paste0("Frecuencia de genes por tipo - ", sample_name),
    x = "Tipo de pertenencia", y = "Número de genes"
  ) +
  theme_minimal(base_size = 14) +  # Aumenta tamaño base
  theme(
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  ) +
  scale_fill_brewer(palette = "Set2")

dev.off()


######### --- Densidad del número de células por gen --- ####
save_path_density <- file.path(save_dirPlot, paste0("Density_nCells_By_Type_", sample_name, ".png"))
png(filename = save_path_density, width = 1000, height = 800)
ggplot(df_resumen_2, aes(x = nCells, fill = type)) +
  geom_density(alpha = 0.5) +
  theme_minimal(base_size = 16) +
  scale_x_log10() +
  scale_fill_paletteer_d("ggsci::default_nejm") +
  labs(
    title = paste0("Distribución del número de células por gen - ", sample_name),
    x = "nCells (escala log10)", y = "Densidad"
  )
dev.off()

######## --- Ranking Scran vs número de células (todo el rango) --- ####
save_path_rank_scran <- file.path(save_dirPlot, paste0("RankScran_vs_nCells_", sample_name, ".png"))
png(filename = save_path_rank_scran, width = 1000, height = 800)
ggplot(df_scran, aes(x = Rank, y = nCells)) +
  geom_point(color = "skyblue") +
  theme_minimal(base_size = 16) +
  labs(
    title = paste0("Ranking Scran vs número de células - ", sample_name),
    x = "Ranking Scran", y = "Número de células"
  )
dev.off()

######## --- Ranking Scran vs número de células (Top 500) --- ####
save_path_rank_scran_top500 <- file.path(save_dirPlot, paste0("Top500_RankScran_vs_nCells_", sample_name, ".png"))
png(filename = save_path_rank_scran_top500, width = 1000, height = 800)
ggplot(df_scran, aes(x = Rank, y = nCells)) +
  geom_point(color = "skyblue") +
  theme_minimal(base_size = 16) +
  xlim(0, 500) +
  ylim(0,200) +
  labs(
    title = paste0("Ranking Scran vs número de células (Top 500) - ", sample_name),
    x = "Ranking Scran", y = "Número de células"
  )
dev.off()

######## --- Ranking SCT vs número de células (todo el rango) --- ####
save_path_rank_sct <- file.path(save_dirPlot, paste0("RankSCT_vs_nCells_", sample_name, ".png"))
png(filename = save_path_rank_sct, width = 1000, height = 800)
ggplot(df_sct, aes(x = Rank, y = nCells)) +
  geom_point(color = "skyblue") +
  theme_minimal(base_size = 16) +
  labs(
    title = paste0("Ranking SCT vs número de células - ", sample_name),
    x = "Ranking SCT", y = "Número de células"
  )
dev.off()

######## --- Ranking SCT vs número de células (Top 500 --- ####
save_path_rank_sct_top500 <- file.path(save_dirPlot, paste0("Top500_RankSCT_vs_nCells_", sample_name, ".png"))
png(filename = save_path_rank_sct_top500, width = 1000, height = 800)
ggplot(df_sct, aes(x = Rank, y = nCells)) +
  geom_point(color = "skyblue") +
  theme_minimal(base_size = 16) +
  xlim(0, 500) +
  ylim(0,200) +
  labs(
    title = paste0("Ranking SCT vs número de células (Top 500) - ", sample_name),
    x = "Ranking SCT", y = "Número de células"
  )
dev.off()

################################################################################
##### Visualización STAT - CV2 y rankings
################################################################################

##### --- CV2 vs número de células (Scran) ---
save_path_cv2_ncells_scran <- file.path(save_dirPlot, paste0("CV2_Scran_vs_nCells_", sample_name, ".png"))
png(filename = save_path_cv2_ncells_scran, width = 1000, height = 800)
ggplot(df_scran, aes(x = CV2_total, y = nCells)) + 
  geom_point(alpha = 0.7, color = "steelblue") + 
  labs(
    title = paste0("CV2 total Scran vs número de células - ", sample_name),
    x = "CV2 total (Scran)", y = "Número de células"
  ) +
  theme_minimal(base_size = 16)
dev.off()

##### --- CV2 vs número de células (SCT) ---
save_path_cv2_ncells_sct <- file.path(save_dirPlot, paste0("CV2_sct_vs_nCells_", sample_name, ".png"))
png(filename = save_path_cv2_ncells_sct, width = 1000, height = 800)
ggplot(df_sct, aes(x = CV2_total, y = nCells)) + 
  geom_point(alpha = 0.7, color = "coral") + 
  labs(
    title = paste0("CV2 total SCT vs número de células - ", sample_name),
    x = "CV2 total (SCT)", y = "Número de células"
  ) +
  theme_minimal(base_size = 16)
dev.off()

##### --- CV2 vs ranking Scran (completo) ---
save_path_cv2_scran_all <- file.path(save_dirPlot, paste0("CV2_vs_RankScran_All_", sample_name, ".png"))
png(filename = save_path_cv2_scran_all, width = 1000, height = 800)
ggplot(df_scran, aes(x = Rank, y = CV2_total)) + 
  geom_point(alpha = 0.7, color = "steelblue") + 
  labs(
    title = paste0("CV2 total vs ranking Scran - ", sample_name),
    x = "Ranking Scran", y = "CV2 total"
  ) +
  theme_minimal(base_size = 16)
dev.off()

# --- CV2 vs ranking Scran (Top 100) ---
save_path_cv2_scran_top <- file.path(save_dirPlot, paste0("CV2_vs_RankScran_Top100_", sample_name, ".png"))
png(filename = save_path_cv2_scran_top, width = 1000, height = 800)
ggplot(df_scran, aes(x = Rank, y = CV2_total)) + 
  geom_point(alpha = 0.7, color = "steelblue") + 
  xlim(0, 100) +
  labs(
    title = paste0("CV2 total vs ranking Scran (Top 100) - ", sample_name),
    x = "Ranking Scran", y = "CV2 total"
  ) +
  theme_minimal(base_size = 16)
dev.off()

##### --- CV2 vs ranking SCT (completo) ---
save_path_CV2_sct_all <- file.path(save_dirPlot, paste0("CV2_vs_RankSCT_All_", sample_name, ".png"))
png(filename = save_path_CV2_sct_all, width = 1000, height = 800)
ggplot(df_sct, aes(x = Rank, y = CV2_total)) + 
  geom_point(alpha = 0.7, color = "coral") + 
  labs(
    title = paste0("CV2 total vs ranking SCT - ", sample_name),
    x = "Ranking SCT", y = "CV2 total"
  ) +
  theme_minimal(base_size = 16)
dev.off()

# --- CV2 vs ranking SCT (Top 100) ---
save_path_CV2_sct_top <- file.path(save_dirPlot, paste0("CV2_vs_RankSCT_Top100_", sample_name, ".png"))
png(filename = save_path_CV2_sct_top, width = 1000, height = 800)
ggplot(df_sct, aes(x = Rank, y = CV2_total)) + 
  geom_point(alpha = 0.7, color = "coral") + 
  xlim(0, 100) +
  labs(
    title = paste0("CV2 total vs ranking SCT (Top 100) - ", sample_name),
    x = "Ranking SCT", y = "CV2 total"
  ) +
  theme_minimal(base_size = 16)
dev.off()

################################################################################
######        Visualización perfil expresión Top/Last genes                #####
################################################################################

# Paso 0: Librerías necesarias
library(scater)
library(patchwork)
library(ggplot2)
library(viridis)
library(dplyr)

##### --- Paso 1: Selección de genes top y last según condiciones de nCells y ranking ---

top_gene_sct <- df_top_sct %>%
  filter(nCells > 20) %>%
  slice(1) %>%
  pull(gene)

top_gene_scran <- df_top_scran %>%
  filter(nCells > 6) %>%
  slice(1) %>%
  pull(gene)

last_gene_sct <- rank_sct %>%
  filter(Rank == max(Rank)) %>%
  pull(gene)

last_gene_scran <- rank_scran %>%
  filter(Rank == max(Rank)) %>%
  pull(gene)

##### --- Paso 2: Validar presencia de genes en la matriz de expresión ---

sce_conv <- Seurat::as.SingleCellExperiment(seurat)  # Conversión Seurat -> SCE
expr_mat_scran <- SummarizedExperiment::assay(sce, "logcounts")
expr_mat_sct <- SummarizedExperiment::assay(sce_conv, "logcounts")

stopifnot(top_gene_sct %in% rownames(expr_mat_sct))
stopifnot(last_gene_sct %in% rownames(expr_mat_sct))
stopifnot(top_gene_scran %in% rownames(expr_mat_scran))
stopifnot(last_gene_scran %in% rownames(expr_mat_scran))

##### --- Paso 3: Función para graficar la expresión de múltiples genes (sin agrupación) ---

plot_expression_single <- function(sce, genes, title, save_path) {
  # Filtrar genes válidos
  genes <- genes[genes %in% rownames(sce)]
  if (length(genes) == 0) {
    warning("Ninguno de los genes proporcionados se encuentra en el objeto sce.")
    return(NULL)
  }
  
  # Crear gráficos individuales
  plots <- lapply(genes, function(gene) {
    expr <- as.numeric(logcounts(sce)[gene, ])
    df <- data.frame(expresion = expr, gen = gene)
    
    ggplot(df, aes(x = gen, y = expresion)) +
      geom_violin(fill = "skyblue", alpha = 0.7, trim = FALSE) +
      geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +
      theme_minimal(base_size = 14) +
      ylab("Expresión log(counts)") +
      xlab("") +
      ggtitle(gene) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
  })
  
  # Combinar y guardar
  final_plot <- wrap_plots(plots, ncol = 2) + plot_annotation(title = title)
  ggsave(save_path, plot = final_plot, width = 10, height = 8, dpi = 300)
}

##### --- Paso 4: Definición de rutas y generación de gráficos ---

genes_sct <- c(top_gene_sct, last_gene_sct)
genes_scran <- c(top_gene_scran, last_gene_scran)

save_path_sct <- file.path(save_dirPlot, paste0("Expresion_SCT_", sample_name, ".png"))
save_path_scran <- file.path(save_dirPlot, paste0("Expresion_Scran_", sample_name, ".png"))

plot_expression_single(sce_conv, genes_sct, paste0("Expresión global genes SCT - ", sample_name), save_path_sct)
plot_expression_single(sce, genes_scran, paste0("Expresión global genes Scran - ", sample_name), save_path_scran)

################################################################################
## COMPARACIÓN DE CONTEOS MEDIOS POR GEN SEGÚN MÉTODO DE NORMALIZACIÓN
################################################################################

# Cálculo de medias de expresión por gen en los datos RAW
counts_mat <- assay(sce, "counts")
gene_means <- rowMeans(counts_mat)
df_gene_means <- data.frame(mean = gene_means)

mean_frequencies_RAW <- df_gene_means %>%
  dplyr::group_by(mean) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n))

# Cálculo de medias por gen para los datos normalizados con SCT
frecuency_counts_gene_sct <- dplyr::select(df_sct, gene, mean = mean)
mean_frequencies_SCT <- frecuency_counts_gene_sct %>%
  dplyr::group_by(mean) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n))

# Cálculo de medias por gen para los datos normalizados con scran
frecuency_counts_gene_scran <- dplyr::select(df_scran, gene, mean = mean)
mean_frequencies_Scran <- frecuency_counts_gene_scran %>%
  dplyr::group_by(mean) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n))

# Añadir columna de método a cada tabla
df_raw <- mean_frequencies_RAW %>% dplyr::mutate(method = "RAW")
df_sct <- mean_frequencies_SCT %>% dplyr::mutate(method = "SCT")
df_scran <- mean_frequencies_Scran %>% dplyr::mutate(method = "scran")

# Combinar en un solo dataframe
combined_df <- dplyr::bind_rows(df_raw, df_sct, df_scran)

# Guardar CSV con los conteos medios por gen y método
save_path <- file.path(save_dirFiles, paste0("Comparacion_conteos_medios_", sample_name, ".csv"))
readr::write_csv(combined_df, save_path)
