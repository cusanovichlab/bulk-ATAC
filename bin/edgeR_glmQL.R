# DESeq pipeline
# setup work directory
setwd("C:/Users/hzhan/Documents/Dr.Cusanovich/data/bulk_ATAC/20230306_glo1ko_kla")

library(limma)
library(edgeR)
library(RColorBrewer)
library(DESeq2)
library(tidyverse)
library(plyr)
#library(org.Hs.eg.db)
library(pheatmap)
# library(calibrate)
library(Signac)
library(ggbeeswarm)
library(ggrepel)
library(grid)
library(gridExtra)
library(lattice)
library(limma)
library(edgeR)
library(UpSetR)
options(scipen = 100)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plotPCA.custom <- function (object, intgroup = "group", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 2], PC2 = pca$x[, 3], group = group, 
                  intgroup.df, name = rownames(colData(object)))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}

runPCA = function(mat, group, subgroup = NA, pca_feature = 5000,
                  coldata){
  coldata$condition = coldata[, group,drop = T]
  peaks <- StringToGRanges(regions = rownames(mat), sep = c("_", "_"))
  
  if (is.na(subgroup)) {
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                  colData = coldata,
                                  design = ~ condition,
                                  rowRanges = peaks)
  } else {
    coldata$subgroup = coldata[, subgroup,drop = T]
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                  colData = coldata,
                                  design = ~ condition + subgroup,
                                  rowRanges = peaks)
  }
  # fit Negative Binomial model
  dds <- DESeq(dds)
  Rlog = rlog(dds)
  dfrlog = assay(Rlog)
  cor.mat = cor(dfrlog)
  
  sampleDists <- as.dist(1-cor.mat)
  n = nrow(coldata)
  
  pheatmap_col = coldata[, group, drop = F]
  message("Correlation heatmap...")
  corplot = pheatmap(cor.mat,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     #breaks = seq(0.8, 1, length = 101),
                     cellwidth = 12,
                     cellheight = 12,
                     clustering_method = "ward.D2",
                     fontsize = 12,
                     annotation_row = pheatmap_col,
                     annotation_col = pheatmap_col)

  ############################
  # PCA plot
  message("Run PCA...")
  pcaData <- plotPCA(Rlog, intgroup=c("condition"), ntop = pca_feature, returnData=TRUE)
  # pcaData_pc23 <- plotPCA.custom(Rlog, intgroup=c("condition"), ntop = 1000, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  # percentVar_pc23 <- round(100 * attr(pcaData_pc23, "percentVar"))
  
  pcaplot = ggplot(pcaData, aes(PC1, PC2, label = name)) +
    geom_point(aes(color=condition), size=5) +
    geom_text(hjust =0, size=5) +
    #scale_color_manual(values = col) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_classic() +
    theme(axis.ticks.length.y = unit(.15, "cm"),
          axis.ticks.length.x = unit(.15, "cm"),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=16),
          axis.ticks = element_line(colour = "black", linewidth = 1),
          axis.title = element_text(size=16),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16))
  return(list(corplot, pcaplot))
}

##############################################
# set variables
output = "glo2only"
# sample = "_"
# load raw counts
seqdata = read_tsv("glomq_rawcounts.txt")
seqdata$peak = paste(seqdata$`#'chr'`, seqdata$`'start'`, seqdata$`'end'`, sep = "_")

colnm = colnames(seqdata)
colnm = gsub(".nodups.bam", "", colnm)
colnm = gsub("'", "", colnm)
colnm = gsub("#", "", colnm)
colnames(seqdata) = colnm

seqdata = seqdata %>%
  select(-chr, -start, -end) %>%
  column_to_rownames("peak")

# build meta data
meta = data.frame(sample = colnames(seqdata))
meta$temp = gsub("kla", "kla_", 
                 gsub("_v", "_vech_", meta$sample))
meta = meta %>%
  separate(temp, c("genotype", "treat", "rep"))
meta = meta %>%
  column_to_rownames("sample")

all(rownames(meta) %in% colnames(seqdata))

meta$sample = rownames(meta)
meta = meta %>%
  filter(genotype != "glo1")
countdata = as.matrix(seqdata)
countdata = countdata[,rownames(meta)]
colnames(countdata)
all(rownames(meta) == colnames(countdata))

# if not, do following
# countdata <- countdata[, coldata$samples]
# all(coldata$samples == colnames(countdata))

# calculate cpm
myCPM <- cpm(countdata)

# # check cpm value that is equal to 10 reads
# plot(myCPM[,1],countdata[,1],ylim=c(0,20),xlim=c(0,3))
# abline(v=1)
# abline(h=10)

# remove gene with cpm < 0.5 over n-3 samples
# thresh <- myCPM > 1
# keep = rowSums(thresh) > 4

meta$condition = factor(paste(meta$genotype, meta$treat, sep = "_"), 
                        levels = c("wt_vech", "wt_kla",
                                   "glo2_vech", "glo2_kla"))

d <- DGEList(counts=countdata, group=meta$condition)
print(d)

# filter low read counts
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]
unnorm <- cpm(d, log=TRUE)
# normalization
d <- calcNormFactors(d, method = "TMM")
print(d)
lcpm <- cpm(d, log=TRUE)

##########################################################
# DE analysis with edgeR
design.mat <- model.matrix(~0+group, data=d$samples)
# colnames(design.mat) <- gsub("group", "", colnames(design.mat))
colnames(design.mat) <- gsub("group", "", colnames(design.mat))
design.mat
# d <- estimateDisp(d,design.mat)
# prior = 50/(6-2)
# d <- estimateDisp(d, prior.df = prior)
# d <- estimateDisp(d, trend = "none")

message("Estimate Dispersion...")
d <- estimateDisp(d, design.mat, robust=TRUE)
print(d$common.dispersion)
# fit model
# et <- exactTest(d, pair=c("wt","ko"))
message("Fitting QL model...")
fit <- glmQLFit(d, design.mat, robust=TRUE)

########################################
# correlaton heatmap
cor.mat = cor(lcpm)
sampleDists <- as.dist(1-cor.mat)

pheatmap_col = d$samples[, "group", drop = F]
message("Correlation heatmap...")
corplot = pheatmap(cor.mat,
                   clustering_distance_rows=sampleDists,
                   clustering_distance_cols=sampleDists,
                   #breaks = seq(0.8, 1, length = 101),
                   cellwidth = 12,
                   cellheight = 12,
                   clustering_method = "ward.D2",
                   fontsize = 12,
                   annotation_row = pheatmap_col,
                   annotation_col = pheatmap_col)
###########################
# PCA
message("Run PCA...")
pcadata = plotMDS(d, top = 5000,
              col=as.numeric(d$samples$group),
              gene.selection="common")

pcadf = data.frame(PC1 = pcadata$x, PC2 = pcadata$y,
                   sample = rownames(pcadata$distance.matrix.squared),
                   group = d$samples$group)

percentVar <- round(100 * pcadata$var.explained[1:2])
# percentVar_pc23 <- round(100 * attr(pcaData_pc23, "percentVar"))

pcaplot = ggplot(pcadf, aes(PC1, PC2, label = sample)) +
  geom_point(aes(color=group), size=5, alpha = 0.75) +
  geom_text(hjust =0, size=5) +
  #scale_color_manual(values = col) +
  scale_color_brewer(palette = "Paired") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  theme(axis.ticks.length.y = unit(.15, "cm"),
        axis.ticks.length.x = unit(.15, "cm"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        axis.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16))

# pcaplot = runPCA(mat = d$counts, group = "group", coldata = d$samples)
message("Plot PCA..")
pdf(file = paste0(output, "_pca_and_correlation.pdf"), width = 14, height = 5)
grid.arrange(corplot[[4]], pcaplot, nrow = 1)
dev.off()

# plot PCA
png(file = paste0(output, "_edger_qc.png"), width = 800, height = 800)
par(mfrow=c(2,2))
boxplot(unnorm, las=2, main="Unnormalized counts")
boxplot(lcpm, las=2, main="Normalized counts")
plotBCV(d)
plotQLDisp(fit)
dev.off()

########################################################
message("Run DE test...")
# qlf <- glmQLFTest(fit, coef=2)
my.contrasts <- makeContrasts(wt.TvsC = wt_kla-wt_vech, 
                              glo2.TvsC = glo2_kla-glo2_vech, 
                              c.glo2vswt = glo2_vech-wt_vech,
                              t.glo2vswt = glo2_kla-wt_kla, 
                              levels=design.mat)

qlf = vector(mode = "list", length = ncol(my.contrasts))
detable = vector(mode = "list", length = ncol(my.contrasts))

colnames(lcpm) = paste0("lcpm_", colnames(lcpm))

for (i in 1:ncol(my.contrasts)) {
  test = colnames(my.contrasts)[i]
  qlf[[i]] <- glmQLFTest(fit, contrast=my.contrasts[,test])
  res <- as.data.frame(topTags(qlf[[i]], n = nrow(qlf[[i]]$table)))
  res = left_join(res %>% 
                    rownames_to_column("peaks"), 
                  as.data.frame(lcpm) %>%
                    rownames_to_column("peaks"), 
                  by = "peaks")
  res = left_join(res,
                  as.data.frame(d$counts) %>%
                    rownames_to_column("peaks"), 
                  by = "peaks")
  res$contract = test
  detable[[i]] = res
  
  de <- decideTests(qlf[[i]], adjust.method="BH", p.value = 0.05)
  print(summary(de))
}

demerge = bind_rows(detable)

write_csv(demerge, file= paste0(output, "_detable.csv"))

# upset plot 
message("Plotting upset plot up peaks...")
desig_up = demerge %>%
  filter(FDR <= 0.05 & logFC > 0) %>%
  select(peaks, contract, FDR, logFC)

desig_down = demerge %>%
  filter(FDR <= 0.05 & logFC < 0) %>%
  select(peaks, contract, FDR, logFC)

desig = list(up = desig_up, down = desig_down)
mypal = c("#c03728", "#66a182")

for (i in 1:2) {
  desig_split = split(desig[[i]], desig[[i]]$contract)
  peak_list = lapply(desig_split, function(x) x$peaks)
  
  top_idx = order(sapply(peak_list, length, simplify = T), decreasing = T)
  top_name = names(peak_list)[top_idx]
  
  upset_plot = upset(fromList(peak_list), nsets = length(peak_list),
                     sets = rev(top_name), text.scale = 3,
                     order.by = "freq", set_size.show = T, 
                     keep.order = TRUE, point.size = 2.5,
                     number.angles = 45,
                     set_size.scale_max = 40800,
                     main.bar.color = mypal[i], 
                     matrix.color = mypal[i], 
                     mainbar.y.label = paste0("DA ", names(desig)[i], " peaks"))

  message(paste0("plot upset plot for ", names(desig)[i], " peaks..."))
  pdf(file = paste0(output, "_", names(desig)[i], "_peak_upsetplot.pdf"), width = 10, height = 6)
  print(upset_plot)
  dev.off()
}
