# ======================================================
# Meta-análisis de genes comunes entre líneas celulares
# Análisis de intersecciones y enriquecimiento funcional
# Estudio: HN00234479
# Autor: Alejandro Virues Morales
# ======================================================

# -----------------------------
# Limpiar entorno y cargar librerías
# -----------------------------

rm(list = ls())

library(readr)
library(dplyr)
library(VennDiagram)
library(UpSetR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(grid)
require(ggplot2)
library(tidyverse)
library(purrr)
library(SingleCellExperiment)
library(Seurat)
library(scater)
library(patchwork)
library(celda)

# -----------------------------
# Definir rutas
# -----------------------------

base_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study"

# -----------------------------
# Cargar genes top 1000 por línea celular
# -----------------------------

# Vector con los nombres de las muestras
samples <- c("HeLA", "A549", "MCF7", "HUHT")

# Directorios base
base_dir_stat <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files"
base_dir_rank <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/5_Norm_Analyses"

# Función para cargar y combinar los datos para cada muestra
load_and_merge <- function(sample) {
  stat_path <- file.path(base_dir_stat, sample, "Stat/Scran", paste0("Scran_TotalStat_", sample, ".csv"))
  rank_path <- file.path(base_dir_rank, sample, "Files", paste0("Stat_Charateristic_Scran_", sample, ".csv"))
  
  stat_df <- read_csv(stat_path, show_col_types = FALSE)
  rank_df <- read_csv(rank_path, show_col_types = FALSE)
  
  df_final <- stat_df %>%
    left_join(rank_df %>% select(gene, Rank, color_group, nCells), by = "gene") %>%
    select(gene, mean, CV2_total, FF_total, Gini, Rank, color_group, nCells) %>%
    arrange(Rank)
  
  return(df_final)
}

# Cargar todos los dataframes y asignarlos a variables individuales
list2env(
  set_names(
    map(samples, load_and_merge),
    paste0("df_final_", samples)
  ),
  envir = .GlobalEnv
)

# Función para extraer el top 5% de genes (por Rank) y quedarse con la columna gene
get_top_genes <- function(sample) {
  df <- get(paste0("df_final_", sample))
  top_n <- ceiling(nrow(df) * 0.05)
  df %>%
    slice_min(Rank, n = top_n) %>%
    pull(gene)
}

# Crear lista con los vectores de genes top 5%
top_genes_list <- set_names(
  map(samples, get_top_genes),
  paste0("top_genes_", samples)
)

# Opcional: asignarlos como variables individuales en tu entorno
list2env(top_genes_list, envir = .GlobalEnv)

# -----------------------------
# Diagramas Venn y UpSet
# -----------------------------

venn_list <- list(A549 = top_genes_A549, HeLA = top_genes_HeLA, HUHT = top_genes_HUHT, MCF7 = top_genes_MCF7)

png(file.path(base_dir, "6_Meta-analyses/Plot/Venn_Intersection_4cellLines.png"), width = 1200, height = 1000)
venn.plot <- VennDiagram::venn.diagram(
  x = venn_list,
  category.names = names(venn_list),
  filename = NULL,
  fill = c("red", "blue", "green", "purple"),
  alpha = 0.5,
  cat.cex = 1.5,
  cex = 1.5,
  margin = 0.1
)
grid::grid.draw(venn.plot)
dev.off()

png(file.path(base_dir, "6_Meta-analyses/Plot/UpSet_Plot_4cellLines.png"), width = 1200, height = 1000)
upset(fromList(venn_list), sets = names(venn_list), order.by = "freq", keep.order = TRUE, text.scale = c(2,2,2,1.5,2,2))
dev.off()

# -----------------------------
# Intersección de genes entre líneas celulares
# -----------------------------

all_genes <- unique(c(top_genes_HeLA, top_genes_A549, top_genes_HUHT, top_genes_MCF7))
gene_presence <- data.frame(
  gene = all_genes,
  A549 = all_genes %in% top_genes_A549,
  HeLA = all_genes %in% top_genes_HeLA,
  HUHT = all_genes %in% top_genes_HUHT,
  MCF7 = all_genes %in% top_genes_MCF7
)
gene_presence$n_lists <- rowSums(gene_presence[, -1])

genes_in_2 <- filter(gene_presence, n_lists == 2)
genes_in_3 <- filter(gene_presence, n_lists == 3)
genes_in_4 <- filter(gene_presence, n_lists == 4)

write_csv(genes_in_2, file.path(base_dir, "6_Meta-analyses/Files/Genes_en_2_listas.csv"))
write_csv(genes_in_3, file.path(base_dir, "6_Meta-analyses/Files/Genes_en_3_listas.csv"))
write_csv(genes_in_4, file.path(base_dir, "6_Meta-analyses/Files/Genes_en_4_listas.csv"))

cat(
  "Resumen: Longitud conjunto genes intersect\n",
  "===============================================\n",
  sprintf("Genes en 2 listas          : %3d genes\n", length(genes_in_2$gene)),
  sprintf("Genes en 3 listas          : %3d genes\n", length(genes_in_3$gene)),
  sprintf("Genes en 4 listas          : %3d genes\n", length(genes_in_4$gene)),
  sprintf("Genes en Intersect A549    : %3d genes\n", length(top_genes_A549)),
  sprintf("Genes en Intersect HeLA    : %3d genes\n", length(top_genes_HeLA)),
  sprintf("Genes en Intersect HUHT    : %3d genes\n", length(top_genes_HUHT)),
  sprintf("Genes en Intersect MCF7    : %3d genes\n", length(top_genes_MCF7)),
  "===============================================\n\n",
  sep = ""
)

# =============================
# Filtrar genes específicos por línea
# =============================

genes_especificos_A549 <- gene_presence %>%
  filter(A549 == TRUE, n_lists == 1) %>%
  pull(gene)

genes_especificos_HeLA <- gene_presence %>%
  filter(HeLA == TRUE, n_lists == 1) %>%
  pull(gene)

genes_especificos_HUHT <- gene_presence %>%
  filter(HUHT == TRUE, n_lists == 1) %>%
  pull(gene)

genes_especificos_MCF7 <- gene_presence %>%
  filter(MCF7 == TRUE, n_lists == 1) %>%
  pull(gene)


# =======================================
# Función para ejecutar y guardar ORA - GO:BP
# =======================================

realizar_ORA <- function(gene_list, nombre_archivo, titulo_grafico, base_dir, cutoff = 0.25, show_top = 10) {

  
  ##### --- Conversión de símbolos a ENTREZID ---
  gene_ids <- clusterProfiler::bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  ##### --- Enriquecimiento GO (BP) ---
  ego <- clusterProfiler::enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = cutoff,
    readable = TRUE
  )
  
  if (nrow(ego@result) == 0) {
    warning(paste("No se encontraron términos enriquecidos para:", nombre_archivo))
    return(NULL)
  }
  
  ##### --- Guardar CSV ---
  readr::write_csv(ego@result, file.path(base_dir, "6_Meta-analyses/Files", paste0(nombre_archivo, "_ORA_GO_BP.csv")))
  
  ##### --- Guardar PNG ---
  categorias <- min(show_top, nrow(ego@result))
  ruta_plot <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/6_Meta-analyses/Plot"
  nombre_png <- paste0(nombre_archivo, "_ORA_GO_BP.png")
  ruta_completa <- file.path(ruta_plot, nombre_png)
  png(filename = ruta_completa, width = 1000, height = 800)  
  enrichplot::dotplot(ego)
  dev.off()
}


# =======================================
# ORA para todos los conjuntos de genes
# =======================================

realizar_ORA(genes_especificos_A549, "Especificos_A549", "ORA - GO:BP - A549 Específicos", base_dir)
realizar_ORA(genes_especificos_HeLA, "Especificos_HeLA", "ORA - GO:BP - HeLA Específicos", base_dir)
realizar_ORA(genes_especificos_HUHT, "Especificos_HUHT", "ORA - GO:BP - HUHT Específicos", base_dir)
realizar_ORA(genes_especificos_MCF7, "Especificos_MCF7", "ORA - GO:BP - MCF7 Específicos", base_dir)

genes_especificos_MCF7_ids <- clusterProfiler::bitr(genes_especificos_MCF7, 
                                                fromType = "SYMBOL", 
                                                toType = "ENTREZID", 
                                                OrgDb = org.Hs.eg.db)

ego_genes_especificos_MCF7 <- clusterProfiler::enrichGO(
  gene = genes_especificos_MCF7_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  readable = TRUE
)

max_categories <- 10
n_categories <- length(ego_genes_especificos_MCF7@result$Description)
categories_to_show <- min(n_categories, max_categories)
save_dirPlot <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/6_Meta-analyses/Plot"
save_path <- file.path(save_dirPlot, "Especificos_MCF7_ORA_GO_BP.png")
# Seleccionar las top categorías y ordenar por p.adjust
go_results <- ego_genes_especificos_MCF7@result

top_go <- go_results %>%
  arrange(p.adjust) %>%           # Ordenar por p.ajustado (opcional)
  head(n = categories_to_show)    # Seleccionar las top N categorías

# Crear el gráfico con ggplot2
p <- ggplot(top_go, aes(x = reorder(Description, p.adjust), 
                        y = p.adjust, 
                        fill = p.adjust)) +
  geom_bar(stat = "identity", color = "black", fill = "steelblue") +
  coord_flip() +  # Barras horizontales para mejor lectura
  labs(
    title = "ORA - GO:BP - MCF7 Específicos",
    x = "GO Biological Process",
    y = "-log10(p.adjusted)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

# Guardar el gráfico
ggsave(
  filename = save_path,
  plot = p,
  width = 12,  # Ajustar ancho (en pulgadas)
  height = 8,  # Ajustar altura
  dpi = 300    # Resolución
)

# -----------------------------
# Rutas base
# -----------------------------
base_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study"
fig_dir <- file.path(base_dir, "6_Meta-analyses", "Plots")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

lines <- c("HeLA", "A549", "HUHT", "MCF7")

# -----------------------------
# Función para categorizar tamaños
# -----------------------------
categorize_tamaño <- function(sce) {
  sums <- colSums(logcounts(sce))
  breaks <- quantile(sums, probs = seq(0, 1, 0.2), na.rm = TRUE)
  labels <- c("Muy Pequeña", "Pequeña", "Mediana", "Grande", "Muy Grande")
  tamaño <- cut(sums, breaks = breaks, include.lowest = TRUE, labels = labels)
  colData(sce)$tamaño <- tamaño
  colData(sce)$sum <- sums
  return(sce)
}

# -----------------------------
# Función para graficar genes
# -----------------------------
plot_gene_expression <- function(sce, gene, line) {
  p <- scater::plotExpression(sce, features = gene, x = "tamaño", exprs_values = "logcounts", color_by="tamaño") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position = c(0.05,0.85)) +  xlab("Tamaño")+
    ggtitle(paste0("Density plot by - ", gene, " of cell line - ", line)) + scale_color_viridis_d("tamaño")
  
  return(p)
}

plot_umap_gene <- function(sce, list_genes, line) { 
  sce <- runPCA(sce, exprs_values = "logcounts")
  sce <- runUMAP(sce, dimred = "PCA")
  umap <- reducedDim(sce, "UMAP")
  celda::plotDimReduceFeature(logcounts(sce),
                              dim1 = umap[, 1],
                              dim2 = umap[, 2],
                              features =list_genes,
                              zscore = FALSE,
                              ncol = 3)
}

# -----------------------------
# Procesar cada línea celular
# -----------------------------
for (line in lines) {
  message("Procesando línea celular: ", line)
  
  # Cargar estadísticas
  stat_path <- file.path(base_dir, "3_Rank", "Files", line, paste0("Stats_TODOS_Genes_SCT_Scran_Cells_", line, ".csv"))
  stats <- readr::read_csv(stat_path)
  
  # Cargar SCE
  sce_path <- file.path(base_dir, "2_Normalization", "Files", paste0(line, "_scran_object_normalized.rds"))
  sce <- readRDS(sce_path)
  sce <- categorize_tamaño(sce)
  
  # Obtener top 10 genes por rank_scran y rank_sct
  top_scran <- stats |> arrange(rank_scran) |> slice_head(n = 9)
  # Graficar y guardar
  for (i in 1:9) {
    g1 <- plot_gene_expression(sce, top_scran$gene[i], line)
    final_plot <- g1
    
    file_name <- paste0("V2_TopGene_", line, "_", i, "_", top_scran$gene[i], ".png")
    ggsave(filename = file.path(fig_dir, line, file_name), plot = final_plot, width = 12, height = 5)
  }
  g3 <- plot_umap_gene(sce, top_scran$gene, line) +
    labs(title = paste0("UMAP - ", line))
  
  final_plot_v2 <- (g3)
    patchwork::plot_annotation(
      title = paste0("UMAP - ", line)
    )
  
  file_name <- paste0("UMAP_TopGenes_", line,".png")
  ggsave(filename = file.path(fig_dir, line, file_name), plot = final_plot_v2, width = 12, height = 5)
}


#### META-RANK #####
##### Cargar datos

# HeLA
CV2_HeLA <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/HeLA/Stat/Scran/Scran_HeLA_CV2.csv")
FF_HeLA <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/HeLA/Stat/Scran/Scran_HeLA_FF.csv")
GI_HeLA <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/HeLA/Stat/Scran/Scran_HeLA_GI.csv")

# A549
CV2_A549 <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/A549/Stat/Scran/Scran_A549_CV2.csv")
FF_A549 <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/A549/Stat/Scran/Scran_A549_FF.csv")
GI_A549 <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/A549/Stat/Scran/Scran_A549_GI.csv")

# HUHT
CV2_HUHT <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/HUHT/Stat/Scran/Scran_HUHT_CV2.csv")
FF_HUHT <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/HUHT/Stat/Scran/Scran_HUHT_FF.csv")
GI_HUHT <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/HUHT/Stat/Scran/Scran_HUHT_GI.csv")

# MCF7
CV2_MCF7 <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/MCF7/Stat/Scran/Scran_MCF7_CV2.csv")
FF_MCF7 <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/MCF7/Stat/Scran/Scran_MCF7_FF.csv")
GI_MCF7 <- read_csv("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files/MCF7/Stat/Scran/Scran_MCF7_GI.csv")

# Renombrar columnas
CV2_HeLA <- CV2_HeLA %>% rename(CV2_HeLA = CV2_total)
CV2_A549 <- CV2_A549 %>% rename(CV2_A549 = CV2_total)
CV2_HUHT <- CV2_HUHT %>% rename(CV2_HUHT = CV2_total)
CV2_MCF7 <- CV2_MCF7 %>% rename(CV2_MCF7 = CV2_total)

FF_HeLA <- FF_HeLA %>% rename(FF_HeLA = FF_total)
FF_A549 <- FF_A549 %>% rename(FF_A549 = FF_total)
FF_HUHT <- FF_HUHT %>% rename(FF_HUHT = FF_total)
FF_MCF7 <- FF_MCF7 %>% rename(FF_MCF7 = FF_total)

GI_HeLA <- GI_HeLA %>% rename(GI_HeLA = Gini)
GI_A549 <- GI_A549 %>% rename(GI_A549 = Gini)
GI_HUHT <- GI_HUHT %>% rename(GI_HUHT = Gini)
GI_MCF7 <- GI_MCF7 %>% rename(GI_MCF7 = Gini)

# Lista de todos los data frames
lista_dfs <- list(
  CV2_HeLA, CV2_A549, CV2_HUHT, CV2_MCF7,
  FF_HeLA, FF_A549, FF_HUHT, FF_MCF7,
  GI_HeLA, GI_A549, GI_HUHT, GI_MCF7
)

merged_df <- lista_dfs %>% 
  purrr::reduce(full_join, by = "gene")

final_matrix <- merged_df %>% 
  na.omit()

# Crear vector de clases (4 líneas celulares repetidas 3 veces)
cl <- c(1,1,1,1,1,1,1,1,1,1,1,1)

# Extraer nombres de genes y matriz de expresión
gene_names_scran <- final_matrix$gene
expr_matrix <- as.matrix(final_matrix[, -1])

# Ejecutar Rank Products
RP_scran.out <- RankProducts(expr_matrix, cl, gene.names = gene_names_scran, rand = 123)

# Extraer métricas
ranking <- RP_scran.out$RPrank[, "class1 > class2"]
rprank <- RP_scran.out$RPs[, "class1 > class2"]
pfp_vals <- RP_scran.out$pfp[, "class1 > class2"]

# Crear dataframe con resultados
df_metarank <- data.frame(
  gene = names(ranking),
  Rank = ranking,
  RP = rprank,
  PFP = pfp_vals
)

# Ordenar resultados finales por Rank_Adjust
metarankinf_final <- df_metarank[order(df_metarank$Rank), ]

save_path <- file.path(base_dir, paste0("6_Meta-analyses/Files/Ranking_Meta-analisis.csv"))
readr::write_csv(metarankinf_final, save_path)
