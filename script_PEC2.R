### 7.1 Extraer datos

counts = read.csv("counts.csv",header=T,sep = ";")
target = read.csv("targets.csv",header=T,sep = ",")

NIT <- subset(target,grepl("^(NIT)", target$Group))
SFI <- subset(target,grepl("^(SFI)", target$Group))
ELI <- subset(target,grepl("^(ELI)", target$Group))

set.seed(12)

NIT_10 <- sample_n(NIT, 10)
SFI_10 <- sample_n(SFI, 10)
ELI_10 <- sample_n(ELI, 10)

columnas <- rbind(NIT_10,SFI_10,ELI_10)
grupos <- as.factor(columnas$Group)
colNIT_SFI <- rbind(NIT_10,SFI_10)
g1 <- factor(colNIT_SFI$Group)
colSFI_ELI <- rbind(SFI_10,ELI_10)
g2 <- factor(colSFI_ELI$Group)
colNIT_ELI <- rbind(NIT_10,ELI_10)
g3 <- factor(colNIT_ELI$Group)

nit <- counts[,NIT_10$Sample_Name]
sfi <- counts[,SFI_10$Sample_Name]
eli <- counts[,ELI_10$Sample_Name]

datos <- cbind(nit,sfi,eli)
rownames(datos) <- counts[,1]

write.csv(datos,"datos.csv", row.names = FALSE)


### 7.2 Filtrado

table(is.na(datos))

med_gen<-apply(datos, 1, mean)
table(med_gen == 0)
boxplot(med_gen)

borrar<-datos[which(med_gen ==0),]
i<-intersect(rownames(borrar), rownames(datos))
datos<-datos[!rownames(datos)%in% i,]

datos_sinnorm <- datos
media_genes_elim<-apply(datos_sinnorm, 1, mean)
barplot(media_genes_elim,main = 'Barplot expresión cruda', xlim=NULL, xlab = 'Genes', ylab='Frecuencia')


### 7.3 Normalización

datos <- normalizeCounts(datos)
maPlot(datos[,1], datos[,2],
       pch=19, cex=.5, ylim=c(-8,8), 
       allCol="darkgrey", lowess=TRUE)
grid(col="black")
title("normalización TMM")

par(mfrow=c(1,2))

media_genes_elim<-apply(datos_sinnorm, 1, mean)
barplot(media_genes_elim,main = 'Barplot expresión cruda', xlim=NULL, xlab = 'Genes', ylab='Frecuencia')

media_genes_nor<-apply(datos, 1, mean)
barplot(media_genes_nor,main = 'Barplot expresión normalizada', xlim=NULL, xlab = 'Genes', ylab='Frecuencia')

### 7.4 Identificación de genes diferencialmente expresados (DE)

primer <- cbind(datos[,1:10],datos[,11:20])
segun <- cbind(datos[,11:20],datos[,21:30])

# Se aplica un refiltrado
d <- DGEList(counts = primer, group = g1)
d <- calcNormFactors(d)
m <- sweep(d$counts, 2, 1e6 / d$samples$lib.size, '*')
ridx_primer <- rowSums(m>1) >= 2
table(ridx_primer)
d <- d[ridx_primer,]
plotMDS(d)

d1 <- estimateCommonDisp(d)
dtag1 <- estimateTagwiseDisp(d1)
res.common1  <- exactTest(d1, pair=c("NIT", "SFI"), dispersion="common")
res.tagwise1 <- exactTest(dtag1, pair=c("NIT", "SFI"), dispersion="tagwise")

plotBCV(dtag1, cex=0.4)

dec1 <- decideTestsDGE(res.common1,p=0.001, adjust="BH")
dtag_primer <- rownames(d1)[as.logical(dec1)]
plotSmear(res.common1, de.tags = dtag_primer)
abline(h=c(-1,1),col="blue")

# así se calculan los genes up/down regulados:
d1_df <- as.data.frame(dec1)
d1_df[,2] <- rownames(res.common1)
up1 <- d1_df[which(d1_df$`SFI-NIT`==1),]
down1 <- d1_df[which(d1_df$`SFI-NIT`==-1),]
summary(dec1)
genes_1 <- c(up1$V2,down1$V2)

# Se aplica un refiltrado
d <- DGEList(counts = segun, group = g2)
d <- calcNormFactors(d)
m <- sweep(d$counts, 2, 1e6 / d$samples$lib.size, '*')
ridx_segun <- rowSums(m>1) >= 2
table(ridx_segun)
d <- d[ridx_segun,]
plotMDS(d)

d2 <- estimateCommonDisp(d)
dtag2 <- estimateTagwiseDisp(d2)
res.common2  <- exactTest(d2, pair=c("SFI", "ELI"), dispersion="common")
res.tagwise2 <- exactTest(dtag2, pair=c("SFI", "ELI"), dispersion="tagwise")

plotBCV(dtag2, cex=0.4)

dec2 <- decideTestsDGE(res.common2,p=0.001, adjust="BH")
dtag_segun <- rownames(d2)[as.logical(dec2)]
plotSmear(res.common2, de.tags = dtag_segun)
abline(h=c(-1,1),col="blue")

# así se calculan los genes up/down regulados:
d2_df <- as.data.frame(dec2)
d2_df[,2] <- rownames(res.common2)
up2 <- d2_df[which(d2_df$`ELI-SFI`==1),]
down2 <- d2_df[which(d2_df$`ELI-SFI`==-1),]
summary(dec2)
genes_2 <- c(up2$V2,down2$V2)

# Se aplica un refiltrado
d <- DGEList(counts = terc, group = g3)
d <- calcNormFactors(d)
m <- sweep(d$counts, 2, 1e6 / d$samples$lib.size, '*')
ridx_terc <- rowSums(m>1) >= 2
table(ridx_terc)
d <- d[ridx_terc,]
plotMDS(d)

d3 <- estimateCommonDisp(d)
dtag3 <- estimateTagwiseDisp(d3)
res.common3  <- exactTest(d3, pair=c("NIT", "ELI"), dispersion="common")
res.tagwise3 <- exactTest(dtag3, pair=c("NIT", "ELI"), dispersion="tagwise")

plotBCV(dtag3, cex=0.4)

dec3 <- decideTestsDGE(res.common3,p=0.001, adjust="BH")
dtag_terc <- rownames(d3)[as.logical(dec3)]
plotSmear(res.common3, de.tags = dtag_terc)
abline(h=c(-1,1),col="blue")

# así se calculan los genes up/down regulados:
d3_df <- as.data.frame(dec3)
d3_df[,2] <- rownames(res.common3)
up3 <- d3_df[which(d3_df$`ELI-NIT`==1),]
down3 <- d3_df[which(d3_df$`ELI-NIT`==-1),]
summary(dec3)
genes_3 <- c(up3$V2,down3$V2)

vd <- venn.diagram(x = list("Genes DE NIT-SFI" = genes_1, "Genes DE SFI-ELI" = genes_2,"Genes DE NIT-ELI" = genes_3), fill = brewer.pal(3, "Pastel2"), filename = NULL)
grid.draw(vd)

comunes <- intersect(intersect(genes_1,genes_2),genes_3)
comunes
length(comunes)

up<-intersect(intersect(up1$V2,up2$V2),up3$V2)
length(up)

down<-intersect(intersect(down1$V2,down2$V2),down3$V2)
length(down)


### 7.5 Anotación de los resultados

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

genes_comunes<-getBM(attributes =  c("hgnc_symbol"), filters = "ensembl_gene_id_version", values =comunes, mart = mart)
genes_comunes <- genes_comunes[,1]
genes_comunes

ups<-getBM(attributes =  c("hgnc_symbol","description"), filters = "ensembl_gene_id_version", values =up, mart = mart)
ups

downs<-getBM(attributes =  c("hgnc_symbol","description"), filters = "ensembl_gene_id_version", values =down, mart = mart)
downs


### 7.6 Análisis de significación biológica (Gene Enrichment Analysis)

#Debemos calcular también las anotaciones para el total de los genes con los que se ha trabajado
genes_totales<-getBM(attributes =  c("hgnc_symbol","go_id"), filters = "ensembl_gene_id_version", values =rownames(datos), mart = mart)
head(genes_totales,10)

# Eliminamos los que no cuenten con GO, para que no interfieran
borr <- which(genes_totales$go_id=="")
dim(genes_totales)
genes_totales <- genes_totales[-borr,]
dim(genes_totales)

# convertir en lista para poder calcular el objeto topGo
l_genes <- unique(genes_totales$hgnc_symbol)
lista <- list()
for (i in l_genes) {
  lista[[i]] = genes_totales[which(genes_totales$hgnc_symbol==i),]$go_id
}
head(lista,2)

gen <- names(lista)
gen_comparado <- factor(as.integer(gen %in% genes_comunes))
table(gen_comparado)
names(gen_comparado) <- gen

GO_data <- new("topGOdata", ontology="BP", allGenes=gen_comparado,annot = annFUN.gene2GO, gene2GO = lista)

# se le acplicará el test de Fisher al objeto topGo
resFisher = runTest(GO_data, algorithm = 'classic', statistic = 'fisher')
resFisher

Nodes = 25
allRes = GenTable(GO_data, classicFisher = resFisher, topNodes = Nodes)
head(allRes)


# Plots
plotEnrich = function(allRes, title){
  # Plotting!
  layout(t(1:2), widths = c(8,1))
  par(mar=c(4, .5, .7, .7), oma = c(3, 15, 3, 4), las = 1)
  
  rbPal = colorRampPalette(c('red', 'white', 'blue'))
  pvalue = as.numeric(gsub("<", "", allRes$classicFisher))
  max_value = as.integer(max(-log(pvalue))) + 1
  pv_range = exp(-seq(max_value, 0, -1))
  allRes$Color = rbPal(max_value) [cut(pvalue, pv_range)]
  
  o = order(allRes$Significant, decreasing = T)
  barplot(allRes$Significant[o], names.arg = allRes$Term[o], las = 2, horiz = T, col = allRes$Color[o],
          xlab = "Number of sequences", main = title, cex.names = 0.85)
  
  image(0, seq(1, max_value), t(seq_along(seq(1, max_value))), col = rev(rbPal(max_value)), axes = F, ann = F)
  pv_label = exp(-seq(log(1), -log(min(pvalue)), l = 6))
  pv_label = formatC(pv_label, format = "e", digits = 2)
  axis(4, at = seq(1, max_value, length = 6), labels = c(1, pv_label[2:6]), cex.axis = 0.85)
  title("p.value", cex.main = 0.6)
}

plotEnrich(allRes = allRes, title = 'Enrichment Analysis')