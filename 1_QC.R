# ==============================================================================
# Control de Calidad (QC) para cada una de las  muestras del estudio
# Estudio: HN00234479
# Autor: Alejandro Virues Morales
# ==============================================================================

# Limpiar el entorno
rm(list = ls())

# Cargar librerías necesarias
library(SingleCellExperiment)
library(ggplot2)
library(DT)
library(RColorBrewer)
library(ComplexHeatmap)
library(htmltools)
library(BiocParallel)
library(scater)
library(patchwork)
library(scDblFinder)

# ==================================================================
# Carga objeto SingleCellExperiment
# ==================================================================
cargar_sce <- function(file_path) {
  if (!file.exists(file_path)) stop("El archivo no existe.")
  sce <- readRDS(file_path)
  if (!inherits(sce, "SingleCellExperiment")) stop("El archivo no contiene un objeto SingleCellExperiment.")
  return(sce)
}

# Definir ruta del archivo y cargar SCE
file_path <- "/home/avirues/TFM/Data/HN00234479/Study/0_CreateSCE/Results/sce_HeLA_raw_emptyDrops.rds"
file_path <- "/home/avirues/TFM/Data/HN00234479/Study/0_CreateSCE/Results/sce_A549_raw_emptyDrops.rds"
file_path <- "/home/avirues/TFM/Data/HN00234479/Study/0_CreateSCE/Results/sce_HUHT_raw_emptyDrops.rds"
file_path <- "/home/avirues/TFM/Data/HN00234479/Study/0_CreateSCE/Results/sce_MCF7_raw_emptyDrops.rds"


sce <- cargar_sce(file_path)
sample_name <- sub("^sce_([^_]+)_.*$", "\\1", basename(file_path))

# Dimensiones iniciales
dim_ini <- dim(sce)
print(paste("Dimensiones iniciales:", dim_ini[1], "genes *", dim_ini[2], "células"))

# Crear directorios si no existen
dirs <- file.path("/home/avirues/TFM/Data/HN00234479/Study/1_QC", c("", "Plots", "Files"))
lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# Identificación de genes mitocondriales
mito_genes <- rowData(sce)[grep("^MT-", rowData(sce)$SYMBOL), "SYMBOL"]

# Cálculo de métricas de QC
sce <- scater::addPerCellQC(sce, subsets = list(Mito = mito_genes), BPPARAM = MulticoreParam(10))

# ==================================================================
# Visualización QC antes del filtrado
# ==================================================================

# Función modificada para incluir líneas horizontales
plot_qc_violin <- function(sce, y_var, title, y_label) {
  p <- ggplot(as.data.frame(colData(sce)), aes(x = Sample, y = .data[[y_var]])) +
    geom_violin(fill = "lightblue", alpha = 0.6) +
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red") +
    stat_summary(fun = median, geom = "text", vjust = -1.5, color = "black", aes(label = round(..y.., 2))) +
    geom_hline(yintercept = c(400, 700, 1500), linetype = "dashed", color = "blue") + 
    ggtitle(title) +
    ylab(y_label) +
    theme_minimal()
  return(p)
}

p1 <- plot_qc_violin(sce, "sum", "ARN total antes del filtrado", "Genes Detectados")
p2 <- plot_qc_violin(sce, "detected", "Genes Detectados antes del filtrado", "Genes Detectados")

# Número de células con más del 50% de genes mitocondriales
cells_above_20 <- sum(colData(sce)$subsets_Mito_percent > 20)

# Modificar gráfico de % genes mitocondriales
p3 <- ggplot(as.data.frame(colData(sce)), aes(x = Sample, y = .data[["subsets_Mito_percent"]])) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red") +
  stat_summary(fun = median, geom = "text", vjust = -1.5, color = "black", aes(label = round(..y.., 2))) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") +  # Línea en 50%
  annotate("text", x = 1, y = 55, label = paste("Células ->", cells_above_20), color = "red", size = 5) + # Número de células sobre el 50%
  ggtitle("% de genes mitocondriales antes del filtrado") +
  ylab("% Genes Mitocondriales") +
  theme_minimal()

p4 <- scater::plotColData(sce, y = "sum", x = "detected") +
    ggtitle("Tamaño librerías x número de genes detectados") +
    xlab("Número de genes detectados") +
    ylab("Tamaño de la librería (UMIs)") +
    theme_minimal()

# Guardar gráficos
jpeg(file = file.path(dirs[2], paste0(sample_name, "_QC_Before_Filtering_emptyDrops.jpeg")), width = 12, height = 10, units = "in", res = 300)
print((p1 | p2) / (p3 | p4) + plot_annotation(title = paste0("QC Metrics - ", sample_name, " (Before Filtering)")))
dev.off()


# Gráfico distribución librerias
jpeg(file = file.path(dirs[2], paste0(sample_name, "_UMIs_Distribution_Before_Filtering_emptyDrops.jpeg")), width = 12, height = 10, units = "in", res = 300)
ggplot(as.data.frame(sce$sum), aes(x = sce$sum)) + 
  geom_histogram(binwidth = 100, fill = "lightblue", color = "black") +
  ggtitle("Distribución del tamaño de la librería antes del filtrado") +
  xlab("Número de genes detectados") +
  ylab("Tamaño de la librería (UMIs)")
dev.off()

# ==================================================================
# Filtrado de células
# ==================================================================

# Identificación de dobletes
bp <- BiocParallel::MulticoreParam(3, RNGseed=1234)
sce <- scDblFinder::scDblFinder(sce, samples = sce$Sample, BPPARAM=bp)


thresholds <- list(umi = 400, feature = 1000, mito = 20)
sizelib.sup <- scater::isOutlier(sce$sum, log = TRUE, nmads = 3, type = "higher")

sizelib.ft <- sce$sum < thresholds$umi
nexprs.ft <- sce$detected < thresholds$feature
rmito.ft <- sce$subsets_Mito_percent > thresholds$mito
doublet <- sce$scDblFinder.class == "doublet"

discard.ft <- sizelib.ft | nexprs.ft | rmito.ft | doublet | sizelib.sup

# Crear data frame con el resumen del filtrado
discard_summary <- data.frame(
  Filtro = c("Libsize", "NExprs", "MitoProp", "Doublet", "Libsize_sup", "Total"),
  Cantidad = c(sum(sizelib.ft), sum(nexprs.ft), sum(rmito.ft), sum(doublet), sum(sizelib.sup), sum(discard.ft))
)

cat("Resumen del filtrado:")
print(discard_summary)

# Guardar en CSV
out_file <- file.path("/home/avirues/TFM/Data/HN00234479/Study/1_QC/Files", paste0("Cell_discard_", sample_name, "_emptyDrops.csv"))
write.csv(discard_summary, file = out_file, row.names = FALSE)

cat("Resumen del filtrado guardado en: ", out_file)

sce_gfilt <- sce
colData(sce_gfilt)$discard <- discard.ft
sce <- sce[, !discard.ft]

n_cell_discard <- sum(sce_gfilt$discard)
cat("Células eliminadas: ", n_cell_discard)
rm(sce_gfilt)
cat("Dimensiones después del filtrado: ", paste(dim(sce), collapse = " x "))

# ==================================================================
# Filtrado de genes
# ==================================================================

non_coding_rnas <- list(
  tRNA = grep("^TR[A-Z]-", rowData(sce)$SYMBOL, value = TRUE),
  miRNA = grep("^(MIR|mir-)", rowData(sce)$SYMBOL, value = TRUE),
  lncRNA = grep("(^LINC[0-9]{5}|-(AS|[DIO]T)[0-9]?$)", rowData(sce)$SYMBOL, value = TRUE),
  smRNA = grep("^RNU", rowData(sce)$SYMBOL, value = TRUE),
  snoRNA = grep("^(SNOR[DA]|SCARNA)", rowData(sce)$SYMBOL, value = TRUE),
  rRNA = grep("^RNA", rowData(sce)$SYMBOL, value = TRUE),
  Pseudogenes = grep("^RP[0-9]{2}-", rowData(sce)$SYMBOL, value = TRUE)
)
cat("Resumen de genes no codificantes:")
lapply(non_coding_rnas, length)

# Guardar resumen de genes en CSV
out_file <- paste0("/home/avirues/TFM/Data/HN00234479/Study/1_QC/Files/gene_summary_", sample_name, "_emptyDrops.csv")
write.csv(data.frame(lapply(non_coding_rnas, length)), file = out_file, row.names = FALSE)
cat("Resumen de genes guardado en: ", out_file)

# Filtrar genes de bajo soporte
# Paso 1: Calcular nCells (número de células donde se expresa cada gen)
rowData(sce)$nCells <- Matrix::rowSums(counts(sce) > 0)

# Paso 2: Filtrado inicial mínimo para quitar genes con muy bajo soporte (<= 3 células)
sce_tmp <- sce[rowData(sce)$nCells > 3, ]

# Paso 3: Calcular el primer cuartil sobre este subconjunto más limpio
cuartiles <- summary(rowData(sce_tmp)$nCells)
cell_filter <- cuartiles["1st Qu."]

# Paso 4: Filtrado final usando ese cuartil como umbral dinámico
non_coding <- unlist(non_coding_rnas)
keep_features <- rowData(sce)$nCells > cell_filter
sce <- sce[!(rownames(sce) %in% non_coding) & keep_features, ]



# ==================================================================
# Visualización QC después del filtrado
# ==================================================================
p1 <- plot_qc_violin(sce, "sum", "ARN total después del filtrado", "Genes Detectados")
p2 <- plot_qc_violin(sce, "detected", "Genes Detectados después del filtrado", "Genes Detectados")

# Modificar gráfico de % genes mitocondriales
cells_above_20 <- sum(colData(sce)$subsets_Mito_percent > 20)

p3 <- ggplot(as.data.frame(colData(sce)), aes(x = Sample, y = .data[["subsets_Mito_percent"]])) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "red") +
  stat_summary(fun = median, geom = "text", vjust = -1.5, color = "black", aes(label = round(..y.., 2))) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") +  # Línea en 50%
  annotate("text", x = 1, y = 55, label = paste("Células ->", cells_above_20), color = "red", size = 5) + # Número de células sobre el 50%
  ggtitle("% de genes mitocondriales antes del filtrado") +
  ylab("% Genes Mitocondriales") +
  theme_minimal()

p4 <- scater::plotColData(sce, y = "sum", x = "detected") +
  ggtitle("Tamaño librerías x número de genes detectados") +
  xlab("Número de genes detectados") +
  ylab("Tamaño de la librería (UMIs)") +
  theme_minimal()

# Guardar gráficos comunes
jpeg(file = file.path(dirs[2], paste0(sample_name, "_QC_After_Filtering_emptyDrops.jpeg")), width = 12, height = 10, units = "in", res = 300)
print((p1 | p2) / (p3 | p4) + plot_annotation(title = paste0("QC Metrics - ", sample_name, " (After Filtering)")))
dev.off()

# Gráfico distribución librerias
jpeg(file = file.path(dirs[2], paste0(sample_name, "_UMIs_Distribution_After_Filtering_emptyDrops.jpeg")), width = 12, height = 10, units = "in", res = 300)
ggplot(as.data.frame(sce$sum), aes(x = sce$sum)) + 
  geom_histogram(binwidth = 100, fill = "lightblue", color = "black") +
  ggtitle("Distribución del tamaño de la librería después del filtrado") +
  xlab("Número de genes detectados") +
  ylab("Tamaño de la librería (UMIs)")
dev.off()

# Generar gráfico de los genes más comunes
jpeg(file = file.path(dirs[2], paste0(sample_name, "_most_expressed_SYMBOL_emptyDrops.jpeg")), width = 12, height = 10, units = "in", res = 300)
plotHighestExprs(sce)
dev.off()

# ==================================================================
# Resumen dimensiones
# ==================================================================

dim_final <- dim(sce)

print(paste("Dimensiones iniciales:", dim_ini[1], "genes *", dim_ini[2], "células"))
print(paste("Dimensiones finales:", dim_final[1], "genes *", dim_final[2], "células"))
print(paste("Número total de células eliminadas:", n_cell_discard,"células" ))

# ==================================================================
# Guardado del objeto
# ==================================================================
saveRDS(sce, file.path(dirs[3], paste0(sample_name, "_emptyDrops_processed.rds")))

