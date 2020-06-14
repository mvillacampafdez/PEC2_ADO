
library(DESeq2)
library("magrittr")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library(stats)
library("ggbeeswarm")
library("genefilter")
library(stringr)
library("AnnotationDbi")
library(org.Hs.eg.db)
require("biomaRt")
library(gplots)
library(hgu133plus2.db)
library(ReactomePA)
library(limma)

if(!file.exists("results"))  dir.create("results", recursive=TRUE)

##### PREPARACIÓN DE LOS DATOS #####

# En primer lugar, se leen las dos tablas 'targets.csv' y 'counts.csv':
counts <- read.csv('counts.csv', sep=';')
targets <- read.csv('targets.csv')
dim(counts) # 56202 293
dim(targets) # 292 9

# La primera columna del conjunto 'counts' corresponde con los genes, por lo que se modifica este conjunto de datos para tener esta columna como
# nombres de cada fila:
rownames(counts) <- counts[,1]; counts <- counts[,-1]
dim(counts) # 56202 292
head(rownames(counts)) # [1] "ENSG00000223972.4" "ENSG00000227232.4" "ENSG00000243485.2" "ENSG00000237613.2" "ENSG00000268020.2" "ENSG00000240361.1"

# Ahora ya se tienen el mismo número de columnas en 'counts' (correspondientes a las muestras) que filas en el conjunto 'targets'

#####

# Seleccionar 10 muestras por cada grupo de muestreo:
# Seleccionar las muestras que corresponden a cada grupo:
table(targets$Group) # ELI: 14  NIT: 236  SFI: 42

muestrasELI <- which(targets$Group=='ELI')
length(muestrasELI) # 14

muestrasNIT <- which(targets$Group=='NIT')
length(muestrasNIT) # 236

muestrasSFI <- which(targets$Group=='SFI')
length(muestrasSFI) # 42

set.seed(123456)
(random_ELI <- muestrasELI[order(runif(muestrasELI))][1:10])
set_ELI <- counts[,random_ELI]

set.seed(123456)
(random_NIT <- muestrasNIT[order(runif(muestrasNIT))][1:10])
set_NIT <- counts[,random_NIT]

set.seed(123456)
(random_SFI <- muestrasSFI[order(runif(muestrasSFI))][1:10])
set_SFI <- counts[,random_SFI]

# Juntar los tres conjuntos
datos <- cbind(set_ELI,set_NIT,set_SFI)
dim(datos) # 56202 30

datos <- counts[,c(random_ELI,random_NIT,random_SFI)]
samples <- targets[c(random_ELI,random_NIT,random_SFI),]

ddsMat <- DESeqDataSetFromMatrix(countData = datos,
                                 colData = samples,
                                 design = ~ Group)
ddsMat


##### PREPROCESADO DE LOS DATOS ##### 

dim(ddsMat)
samples$Group
# A partir del conjunto datos, las primeras 10 columnas corresponden con ELI, las siguientes 10 con NIT y las últimas 10 con SFI


### Pre-filtrado del conjunto de datos:
ddsMat <- ddsMat[rowSums(counts(ddsMat)) > 1,]
dim(ddsMat) # 43658 30     

### Transformación estabilizadora de la varianza:
vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

# Transformación rlog
rld <- rlog(ddsMat, blind = FALSE)
head(assay(rld), 3)

# Efecto de la transformación:
dds <- estimateSizeFactors(ddsMat)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))


colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

### Samples distances     
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

### PCA plot       
plotPCA(vsd, intgroup = c('Group'))

### MDS plot      
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Group, shape = Group)) +
  geom_point(size = 3) + coord_fixed()

# Gene clustering (los genes con expresión más variable)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno,show_rownames = F,show_colnames = F)


##### DEG ##### 
dds <- DESeq(dds, parallel =TRUE)

## ELI-NIT      
resEN <- results(dds, contrast=c("Group","ELI","NIT"))
resEN
mcols(resEN, use.names = TRUE)
summary(resEN)

## SFI-NIT
resSN <- results(dds, contrast=c("Group","SFI","NIT"))
resSN
mcols(resSN, use.names = TRUE)
summary(resSN)

## ELI-SFI
resES <- results(dds, contrast=c("Group","ELI","SFI"))
resES
mcols(resES, use.names = TRUE)
summary(resES)

# Considerando una fracción de 10% de falsos positivos: número de DEG:
sum(resEN$padj < 0.1, na.rm=TRUE) # 5872
sum(resES$padj < 0.1, na.rm=TRUE) # 7998
sum(resSN$padj < 0.1, na.rm=TRUE) # 794

# 'Subset' los resultados y ordenar en función de log2FC (up-regulated y down-regulated)
resSigEN <- subset(resEN, padj < 0.1)
head(resSigEN[ order(resSigEN$log2FoldChange), ])
head(resSigEN[ order(resSigEN$log2FoldChange, decreasing = TRUE), ])

resSigES <- subset(resES, padj < 0.1)
head(resSigES[ order(resSigES$log2FoldChange), ])
head(resSigES[ order(resSigES$log2FoldChange, decreasing = TRUE), ])

resSigSN <- subset(resSN, padj < 0.1)
head(resSigSN[ order(resSigSN$log2FoldChange), ])
head(resSigSN[ order(resSigSN$log2FoldChange, decreasing = TRUE), ])

### plotting results:
# Counts plot
topGeneSN <- rownames(resSN)[which.min(resSN$padj)]
topGeneEN <- rownames(resEN)[which.min(resEN$padj)]
topGeneES <- rownames(resES)[which.min(resES$padj)]

# SN
geneCounts <- plotCounts(dds, gene = topGeneSN, intgroup = c("Group"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + ggtitle('Between SFI and NIT')

# EN
geneCounts <- plotCounts(dds, gene = topGeneEN, intgroup = c("Group"),returnData = TRUE)
ggplot(geneCounts, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + ggtitle('Between ELI and NIT')

# ES
geneCounts <- plotCounts(dds, gene = topGeneES, intgroup = c("Group"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Group, y = count, color = Group, group = Group)) +
  scale_y_log10() + geom_point(size = 3) + geom_line() + ggtitle('Between ELI and SFI')


# GENE CLUSTERING

topSigGenes <- unique(c(head(rownames(resSigEN[order(resSigEN$padj),]),10),
                        head(rownames(resSigES[order(resSigES$padj),]),10),
                        head(rownames(resSigSN[order(resSigSN$padj),]),10)))
mat  <- assay(vsd)[topSigGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno, main = 'Genes más significativos de los tres contrastes',show_rownames = F,show_colnames = F)


mat  <- assay(vsd)[head(rownames(resSigEN[order(resSigEN$padj),]),10), ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno, main = 'Genes más significativos ELIvsNIT',show_rownames = F,show_colnames = F)


mat  <- assay(vsd)[head(rownames(resSigES[order(resSigES$padj),]),10), ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno, main = 'Genes más significativos ELIvsSFI',show_rownames = F,show_colnames = F)


mat  <- assay(vsd)[head(rownames(resSigSN[order(resSigSN$padj),]),10), ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno, main = 'Genes más significativos SFIvsNIT',show_rownames = F,show_colnames = F)


##### ANOTACIÓN #####
resEN$symbol <- mapIds(org.Hs.eg.db,
                     keys=gsub("\\.[0-9]*$", "", row.names(resEN)),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

resEN$entrez <- mapIds(org.Hs.eg.db,
                     keys=gsub("\\.[0-9]*$", "", row.names(resEN)),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resENOrdered <- resEN[order(resEN$pvalue),]
head(resENOrdered)


resES$symbol <- mapIds(org.Hs.eg.db,
                       keys=gsub("\\.[0-9]*$", "", row.names(resES)),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

resES$entrez <- mapIds(org.Hs.eg.db,
                       keys=gsub("\\.[0-9]*$", "", row.names(resES)),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

resESOrdered <- resES[order(resES$pvalue),]
head(resESOrdered)


resSN$symbol <- mapIds(org.Hs.eg.db,
                       keys=gsub("\\.[0-9]*$", "", row.names(resSN)),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

resSN$entrez <- mapIds(org.Hs.eg.db,
                       keys=gsub("\\.[0-9]*$", "", row.names(resSN)),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

resSNOrdered <- resSN[order(resSN$pvalue),]
head(resSNOrdered)


## MULTIPLE COMPARISONS

## DIAGRAMA DE VENN
res_EN.genes <- unique(rownames(resEN)[resEN$padj<0.1])
res_ES.genes <- unique(rownames(resES)[resES$padj<0.1])
res_SN.genes <- unique(rownames(resSN)[resSN$padj<0.1])

# Combinación de las tres listas
comb <- unique(c(res_SN.genes,res_EN.genes,res_ES.genes))

# Comparandolas
res_EN.genes.2 <- comb %in% res_EN.genes
res_ES.genes.2 <- comb %in% res_ES.genes  
res_SN.genes.2 <- comb %in% res_SN.genes  

# Generando conteos Venn para el diagrama
counts.genes <- cbind(res_SN.genes.2, res_EN.genes.2,res_ES.genes.2)
results.genes <- vennCounts(counts.genes)
vennDiagram(results.genes, cex = 1,names = c("SFIvsNIT","ELIvsNIT","ELIvsSFI"), circle.col = c("red", "blue","green"))


## HEATMAP
# Selección de los gees significativos para representar la expresión:
res_EN.genes <- unique(rownames(resEN)[resEN$padj<0.01])
res_ES.genes <- unique(rownames(resES)[resES$padj<0.01])
res_SN.genes <- unique(rownames(resSN)[resSN$padj<0.01])
probesInHeatmap <- unique(c(res_SN.genes,res_EN.genes,res_ES.genes))
HMdata <- assay(ddsMat)[rownames(assay(ddsMat)) %in% probesInHeatmap,]

my_palette <- colorRampPalette(c("blue", "red"))(n = 299)
heatmap.2(HMdata,
     Rowv = TRUE,
     Colv = TRUE,
     dendrogram = 'both',
     main = "Genes diferencialmente expresados",
     scale = "row",
     col = my_palette,
     sepcolor = "white",
     sepwidth = c(0.05,0.05),
     cexRow = 0.5,
     cexCol = 0.9,
     key = TRUE,
     keysize = 1.5,
     density.info = "histogram",
     ColSideColors = c(rep("red",10),rep("blue",10), rep("green",10)),
     tracecol = NULL,
     srtCol = 30)


##### ANÁLISIS DE SIGNIFICACIÓN BIOLÓGICA #####

listOfTables <- list(ELIvsNIT = resENOrdered, 
                     ELIvsSFI = resESOrdered, 
                     SFIvsNIT = resSNOrdered)
listOfSelected <- list()

for (i in 1:length(listOfTables)){
  topTab <- listOfTables[[i]]
  listOfSelected[[i]] <- na.omit(topTab$entrez[topTab$padj<0.1])
}
sapply(listOfSelected, length)

mapped_genes2GO <- mappedkeys(org.Hs.egGO)
mapped_genes2KEGG <- mappedkeys(org.Hs.egPATH)
mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)

listOfData <- listOfSelected[1:3]
comparisonsNames <- c('ELIvsNIT','ELIvsSFI','SFIvsNIT')
universe <- mapped_genes
enrichment <- list()
for (i in 1:length(listOfData)){
  genesIn <- listOfData[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichPathway(gene = genesIn,
                                 pvalueCutoff = 0.05,
                                 readable = T,
                                 pAdjustMethod = "BH",
                                 organism = "human",
                                 universe = universe)
  enrichment[[i]] <- enrich.result
  cat("##################################")
  cat("\nComparison: ", comparison,"\n")
  print(head(enrich.result[,-8]))
  
  if (length(rownames(enrich.result@result)) != 0) {
    write.csv(as.data.frame(enrich.result), 
              file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
              row.names = FALSE)
    
    pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
    print(barplot(enrich.result, showCategory = 15, font.size = 4, 
                  title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
    dev.off()
    
    pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
                   vertex.label.cex = 0.75))
    dev.off()
  }
}



