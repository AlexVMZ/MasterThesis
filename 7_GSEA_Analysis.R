# Limpiar el entorno
rm(list = ls())

# Cargamos paquetes necesarios
pacman::p_load(topGO, stringr, dplyr, AnnotationDbi, GO.db, ReactomePA, biomaRt, KEGGREST, fgsea, mdgsa, 
               xlsx, reactome.db, KEGG.db, org.Hs.eg.db, ggplot2, openxlsx)

# Ruta de trabajo
setwd("/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study")

# Argumentos del análisis
args <- list("HeLA", "A549", "HUHT", "MCF7")
tipo_celular <- args[[3]]

ruta <- paste0("3_Rank/Files/", tipo_celular)
ruta_2 <- paste0("7_GSEA_Analysis/",tipo_celular)
dir.create(ruta_2, showWarnings = FALSE)

base_dir <- "/clinicfs/userhomes/avirues/TFM/Data/HN00234479/Study/3_Rank/Files"

rp_dir <- file.path(base_dir, tipo_celular, "RP")

# Cargar RP_Rank_scran
file_scran <- file.path(rp_dir, paste0("RP_Rank_scran_", tipo_celular, ".rds"))
RP_scran <- readRDS(file_scran)
dataset <- RP_scran$three_metrics

#### Meta-Analisis conjunto #####
dataset <- read_csv("~/TFM/Data/HN00234479/Study/3_Rank/Files/Meta-analisis/Ranking_scran_Meta-analisis.csv")

# Conectamos con Ensembl
mart <- biomaRt::useMart(
  biomart = "ensembl", 
  dataset = "hsapiens_gene_ensembl", 
  host = "https://asia.ensembl.org"
)

# Unimos con dataset original
conversion <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(dataset$gene),
  columns = "ENSEMBL",
  keytype = "SYMBOL"
) %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

dataset_converted <- dataset %>%
  left_join(conversion, by = c("gene" = "SYMBOL")) %>%
  mutate(gene = ENSEMBL) %>%        # Reemplazar la columna original
  select(-ENSEMBL)       

######################## ANOTACIÓN GO ########################
go_terms <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "go_id", "name_1006", "namespace_1003"),
                  filter = 'ensembl_gene_id',
                  values = dataset_converted$gene, 
                  mart = mart)
go_terms <- go_terms[go_terms$go_id != "", ]

# Procesamos términos por categoría
procesar_GO <- function(namespace) {
  df <- go_terms[go_terms$namespace_1003 == namespace, 1:4]
  list(
    gentoset = split(df$go_id, df$ensembl_gene_id),
    settogen = split(df$ensembl_gene_id, df$go_id)
  )
}

bp <- procesar_GO("biological_process")
mf <- procesar_GO("molecular_function")
cc <- procesar_GO("cellular_component")

######################## REACTOME ########################
reactome <- getBM(attributes = c('ensembl_gene_id', "hgnc_symbol", 'reactome'),
                  filter = 'ensembl_gene_id',
                  values = dataset_converted$gene,
                  mart = mart)
reactome <- reactome[reactome$reactome != "", ]
genetoreactome <- split(reactome$reactome, reactome$ensembl_gene_id)
reactometogen <- split(reactome$ensembl_gene_id, reactome$reactome)

################################ KEGG ##########################################
pathways.list <- keggList("pathway", "hsa")
kegg_db = data.frame(codigo = names(pathways.list),
                     summary = pathways.list)

kegg_genes <- AnnotationDbi::select(org.Hs.eg.db, keys=dataset_converted$gene,
                                    columns=c("ENSEMBL", "PATH"), 
                                    keytype="ENSEMBL")
kegg_genes = na.omit(kegg_genes)
kegg_genes$PATH = paste0("hsa",kegg_genes$PATH)

kegg = merge(kegg_db, kegg_genes, by.x = "codigo", by.y = "PATH")

gentokegg = split(kegg$codigo, kegg$ENSEMBL)
keggtogen = split(kegg$ENSEMBL, kegg$codigo)

######################## FGSEA ANALYSIS ########################

dataset_converted <- dataset_converted[!is.na(dataset_converted$gene), ]

dataset_converted$ranking <- seq_len(nrow(dataset_converted))
vector_rank <- dataset_converted$ranking
rank_invertido <- rev(vector_rank)
names(rank_invertido) <- dataset_converted$gene

rank_invertido <- rank_invertido[!is.na(names(rank_invertido))]

rank_genes <- sort(rank_invertido, decreasing = TRUE)
min_rank <- min(rank_genes)
max_rank <- max(rank_genes)

# Convertir a z-score en el rango [-1, 1] (mejor ranking tendrá valor 1)
rank_genes <- 2 * (rank_genes - min_rank) / (max_rank - min_rank) -1


lista_universo <- list(bp$settogen, mf$settogen, cc$settogen, reactometogen, keggtogen)
list_of_datasets <- list()
names_lista <- c("BP", "MF", "CC", "REACTOME","KEGG")

for (i in seq_along(lista_universo)) {
  rutas <- lista_universo[[i]]
  resultado <- fgsea(pathways = rutas, 
                     stats = rank_genes,
                     eps = 0.0,
                     minSize = 25,
                     maxSize = 500)
  resultado <- resultado[order(resultado$padj), ]
  list_of_datasets[[i]] <- resultado
}

names(list_of_datasets) <- names_lista

######################## NOMBRES DE TÉRMINOS ########################
list_of_datasets[[1]]$description <- getGOnames(list_of_datasets[[1]]$pathway)

# Reactome
res_reactome <- list_of_datasets[["REACTOME"]]
res_reactome$nombre <- res_reactome$pathway
notacion_reactome <- as.list(reactomePATHID2NAME)
notacion_reactome <- unlist(notacion_reactome)
notacion_reactome <- notacion_reactome[res_reactome$pathway]
df_reactome <- data.frame(ruta = notacion_reactome, nombre = names(notacion_reactome))
res_reactome <- merge(res_reactome, df_reactome, by = "nombre")
res_reactome <- res_reactome[order(res_reactome$padj), ]
list_of_datasets[["REACTOME"]] <- res_reactome

# Obtenemos los nombres de las rutas KEGG
res_kegg = list_of_datasets[["KEGG"]]
res_kegg$nombre = (res_kegg$pathway)

notacion_kegg <- unlist(as.list(KEGGPATHID2NAME))
names(notacion_kegg) = paste0("hsa", names(notacion_kegg))
notacion_kegg = notacion_kegg[(res_kegg$pathway)]

notacion_kegg = data.frame(ruta = notacion_kegg,
                           nombre = names(notacion_kegg))

list_of_datasets[["KEGG"]] = merge(res_kegg, notacion_kegg, by.x = "nombre", by.y = "nombre")
list_of_datasets[["KEGG"]] =   list_of_datasets[["KEGG"]][order(list_of_datasets[["KEGG"]]$padj),]

######################## GUARDAMOS RESULTADOS ########################
openxlsx::write.xlsx(list_of_datasets, file = paste0(ruta_2, "/resultados_fgsea_rank.xlsx"), rowNames = TRUE)
saveRDS(list_of_datasets, file = paste0(ruta_2, "/resultados_fgsea_rank.rds"))

######################## PLOTS ########################
ruta_3 <- paste0(ruta_2, "/plots")
dir.create(ruta_3, showWarnings = FALSE)

for (nombre_cat in names(list_of_datasets)) {
  res <- list_of_datasets[[nombre_cat]]
  enriquecidos <- res[res$padj < 0.05, ]
  ruta_4 <- paste0(ruta_3, "/",nombre_cat)
  dir.create(ruta_4, showWarnings = FALSE)
  if (nrow(enriquecidos) > 0) {
    for (i in 1:nrow(enriquecidos)) {
      id <- enriquecidos$pathway[i]
      # id <- "GO:0006351"
      plot_file <- paste0(ruta_4, "/", gsub("[^a-zA-Z0-9]", "_", id), "_", nombre_cat, "_rank.png")
      tryCatch({
        plot <- plotEnrichment(lista_universo[[which(names_lista == nombre_cat)]][[id]], rank_genes) +
          labs(title = id)
        ggsave(plot_file, plot = plot, width = 15, height = 10, bg = "white")
      }, error = function(e) {
        message("Error plot: ", id)
      })
    }
  }
}

######################## BARPLOTS NES POSITIVOS Y NEGATIVOS ########################
ruta_5 <- paste0(ruta_2, "/plots_barplot")
dir.create(ruta_5, showWarnings = FALSE)

for (nombre_cat in names(list_of_datasets)) {
  res <- list_of_datasets[[nombre_cat]]
  enriquecidos <- res[res$padj < 0.05 & !is.na(res$NES), ]
  
  if (nrow(enriquecidos) > 0) {
    # Ordenar por NES
    enriquecidos <- enriquecidos[order(enriquecidos$NES), ]
    
    # Tomamos 10 más positivos y 10 más negativos
    top_terms <- rbind(
      head(enriquecidos, 10),
      tail(enriquecidos, 10)
    )
    
    # Añadimos columna con nombres (GO o Reactome)
    if (nombre_cat %in% c("BP", "MF", "CC")) {
      top_terms$nombre <- top_terms$description
    } else if (nombre_cat == "REACTOME") {
      top_terms$nombre <- top_terms$ruta
    } else if (nombre_cat == "KEGG") {
      top_terms$nombre <- top_terms$ruta
    }
    
    # Asegurar que NES es numérico
    top_terms$NES <- as.numeric(top_terms$NES)
    
    # Crear columna lógica explícita para dirección
    top_terms$direction <- factor(top_terms$NES > 0, levels = c(FALSE, TRUE))
    
    # Orden para gráfico
    top_terms$nombre <- factor(top_terms$nombre, levels = top_terms$nombre[order(top_terms$NES)])
    
    # Barplot
    p <- ggplot(top_terms, aes(x = nombre, y = NES, fill = direction)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(
        values = c("FALSE" = "red", "TRUE" = "blue"),
        labels = c("FALSE" = "NES < 0", "TRUE" = "NES > 0"),
        drop = FALSE
      ) +
      labs(title = paste("Top términos enriquecidos (", nombre_cat, ") - ", tipo_celular),
           x = "", y = "NES", fill = "Dirección") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Guardar
    ggsave(paste0(ruta_5, "/", nombre_cat, "_barplot_NES.png"),
           plot = p, width = 12, height = 8, bg = "white")
  }
}

