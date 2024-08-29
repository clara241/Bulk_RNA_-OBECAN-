library("limma")
library("edgeR")
library("RColorBrewer")
library("gplots")
library("FactoMineR")
library("factoextra")
library(cowplot)
#library("dplyr")
library(EnhancedVolcano)
library(openxlsx)

setwd("./Desktop/Bulk_RNA_(OBECAN)")


count_matrix_no <- read.csv("./03-Counts/rnaseq_salmon_workshop_counts_no_obese.txt", header = T, sep = "\t")
count_matrix_no <- count_matrix_no[-1,]
dim(count_matrix_no) # number of genes 

count_matrix_ob <- read.csv("./03-Counts/rnaseq_salmon_workshop_counts.txt", header = T, sep = "\t")
count_matrix_ob <- count_matrix_ob[-1,]
dim(count_matrix_ob) # number of genes 

# Unir los dataframes por los nombres de las filas
df_merged <- merge(df1, df2, by = 0, all = TRUE)


nom_filas <- rownames(count_matrix)

# Convertir las columnas a numéricas usando apply
count_matrix <- apply(count_matrix, 2, as.numeric)

# Reasignar los nombres de las filas
rownames(count_matrix) <- nom_filas

d0 <- DGEList(count_matrix)

d0 <- calcNormFactors(d0)

cutoff <- 3
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

logcpm <- cpm(d, prior.count=2, log=TRUE)
write.table(logcpm,"./rnaseq_workshop_normalized_count_matrix.txt",sep="\t",quote=F)

# Definir el diseño experimental
group <- factor(c(c(rep(1, 23)),c(rep(2, 68))))
plotMDS(d, col = as.numeric(group))

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

# Assign the column names
colnames(mm) <- c("No_obese", "Obese")
fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(FlightVsG_control = Obese - No_obese, levels = mm)

tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# Generate a volcano plot to visualize differential expression
volcanoplot(tmp)

top.table <- topTable(tmp, adjust.method = "BH")
head(top.table, 30)

# Save results
write.fit(tmp, file = "./DEG_limma.txt", sep = "\t", adjust = "BH", row.names = T)



#########
## PCA ##
#########
counts_normalized <- read.csv("./rnaseq_workshop_normalized_count_matrix.txt", header = T, sep = "\t")

## Seleccionamos solo los valores de las cuantificaciones normalizadas
count_matrix_df <- as.data.frame(counts_normalized)

# Ver si hay datos vacíos
table(is.na(count_matrix_df))

count_matrix_df_t <- t(count_matrix_df)

##  Añadimos una columna especificando el grupo/condición a la que pertenece cada muestra
count_matrix_matrix_groups <- factor(c(paste(rep("No_obese",23)),paste(rep("Obese",68))))

count_matrix_df_t <- cbind.data.frame(count_matrix_df_t,"Groups"=count_matrix_matrix_groups)

## Antes de hacer el PCA tenemos que quitar la columna de grupo. No se incluye el grupo/condición en el PCA
res.pca <- PCA(X = count_matrix_df_t[,-19346],scale.unit = T, graph = F)

ggpubr::ggpar(fviz_pca_ind(X = res.pca,habillage = count_matrix_df_t$Groups, palette = c("slateblue","firebrick"),pointsize=4,pointshape=,
                           addEllipses = ,legend.title="Genotype",ellipse.type=,mean.point=F,repel = T,
                           geom.ind = c("point","text"), title = "Principal Component Analysis", subtitle = "No obese vs Obese RNA seq",
                           xlab = "PC1", ylab = "PC2", ggtheme = theme_cowplot()))

ggsave("./results/PCA.png", width = 10, height = 8, units = "in", dpi = 300)




########################
## Volcano plot Bayes ##
########################
de_genes <- read.csv("./DEG_limma.txt", sep = "\t")

# Guardar los nombres de las filas
genes <- list()
ID <- list()
for (cadena in de_genes[,1]){
  # Dividir la cadena en partes usando el delimitador '|'
  partes <- strsplit(cadena, "\\|")[[1]]
  # Extraer la primera parte
  ID <- append(ID, partes[1])
  genes <- append(genes, partes[6])
}

rownames(de_genes) <- ID
gene_name <- unlist(genes)
de_genes <- cbind(de_genes, gene_name)

# Eliminamos la columna con los nombres largos
de_genes <- de_genes[,-1]


# Nos aseguramos de que las columnas son numéricas
names(de_genes)[names(de_genes) == "Coef"] <- "logFC"
de_genes$logFC <- as.numeric(de_genes$logFC)
de_genes$P.value.adj <- as.numeric(de_genes$P.value.adj)

# Asignamos colores a los genes de acuerdo a su FC y su pval
keyvals <- ifelse(
  de_genes$logFC < -1 & de_genes$P.value.adj < 0.05, "navy",
  ifelse(de_genes$logFC  > 1 & de_genes$P.value.adj < 0.05,"#D82632",
         "grey"))
names(keyvals)[keyvals == "navy"] <- "Downregulated"
names(keyvals)[keyvals == "#D82632"] <- "Upregulated"
names(keyvals)[keyvals == "grey"] <- "NS"
#de_genes$gene <- rownames(de_genes)
#selected_labs <- c(slice_min(.data = de_genes, order_by = P.value.adj, n = 10)$gene_name)
volcano_plot <- EnhancedVolcano(toptable = de_genes,lab = de_genes$gene_name,
                                #selectLab = selected_labs,
                                x = "logFC",y = "P.value.adj",pCutoff = 0.05,FCcutoff = 1,
                                ylim = range(-log10(de_genes$P.value.adj)), xlim = range(de_genes$logFC), 
                                drawConnectors = TRUE,widthConnectors = 0.5,typeConnectors = ,
                                endsConnectors = 'first',labSize = 4,gridlines.minor = F,gridlines.major = F
                                ,pointSize = 3,colAlpha = 0.5, title = "Volcano plot",
                                subtitle = 
                                "Adj p-value cutoff (dashed line): p<0.05
                                Log2 FC cutoff (dashed line): 1",colCustom = keyvals,
                                labFace = "italic",boxedLabels = TRUE, max.overlaps = Inf)
ggsave(path = "./results/",filename = "Volcano.png",
       plot = volcano_plot,device = png,
       width = 9, height = 6,dpi = 300)



##########################################
## 9) Analysis of Pathways with EnrichR ##
##########################################

library(enrichR)

listEnrichrSites()

#setEnrichrSite("Enrichr")

dbs <- listEnrichrDbs()


#dbs <- c("GO_Biological_Process_2021",
#         "KEGG_2019_Mouse","KEGG_2021_Human","MSigDB_Hallmark_2020","WikiPathways_2019_Mouse","BioPlanet_2019",
#        "Elsevier_Pathway_Collection","OMIM_Disease")

dbs <- c("BioPlanet_2019","KEGG_2019_Mouse","KEGG_2021_Human","WikiPathway_2023_Human","WikiPathways_2019_Mouse",
         "Elsevier_Pathway_Collection","MSigDB_Hallmark_2020","BioCarta_2016","HumanCyc_2016","Panther_2016",
         "OMIM_Disease", "GO_Biological_Process_2021")

# Para los datos de Bayes
de_genes_filt <- de_genes %>% filter(P.value.adj<0.05)

enriched <- enrichr(genes = de_genes_filt$gene_name, dbs)

# Para los datos con Bayes
plotEnrich(df = enriched[[1]], showTerms = 30, numChar = 70, y = "Ratio", orderBy = "Adjusted.P.value")

for (i in 1:length(dbs)) {
  write.xlsx(x = enriched[[i]], 
             file = paste("./results/Enrichment/","_",names(enriched[i]),
                          ".xlsx"))
}

for (i in 1:length(dbs)) {
  title_plot <- paste(names(enriched[i]),"liver",sep = "_")
  enrich_plot <-plotEnrich(df = enriched[[i]],
                           showTerms = 30, numChar = 70,y = "Ratio",orderBy = "P.value",
                           title = title_plot) + theme(text = element_text(size=14)) #+ scale_fill_continuous(limits=c(0,0.01), low="red", high="blue")
  ggsave(plot = enrich_plot, filename = paste("./results/Enrichment/",title_plot,".png",sep = ""),
         dpi = 300, width = 11, height = 8)
}

