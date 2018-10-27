#K-MEANS PROTEOMICS DATA ON 34 PAM50 PROTEINS
library(ggplot2)
library(useful)
library(readr)
PAM50 <- read_csv("Downloads/ATIB M5/PAM50.csv", 
                  col_names = FALSE)
colnames(PAM50) = c('GeneSymbol', 'Gene')
colnames(brca_PROTEIN)[1]="GeneSymbol"
temp = brca_PROTEIN
aa = gsub("-.*","",brca_PROTEIN$GeneSymbol)
brca_PROTEIN$GeneSymbol = aa
mg = merge(x=PAM50, y=brca_PROTEIN, by= "GeneSymbol")
mg =mg[-c(27,28,29,31,34),]
mg = mg[-which(rowMeans(is.na(mg)) > 0.2), ]
# mg = na.omit(mg)
mg2 = mg[,-1]
rownames(mg2) = mg[,1]
mg2 = mg2[,-c(1)]
k = which(is.na(mg2), arr.ind=TRUE)
mg2[k] <- rowMeans(mg2, na.rm=TRUE)[k[,1]]
mg2_scale <- scale(t(mg2))
k1 = kmeans(mg2_scale, 4, iter.max = 100)
plot.kmeans(k1, data=mg2_scale, title = 'K-means proteomics using 34 PAM50 proteins')
kd = data.frame(k1$cluster)
kd$ids <- rownames(kd)

pam_ids <- read.table("~/downloads/ATIB M5/subtypes.txt", header=T, sep='\t')
pam_ids$ids <- gsub("\\..*", "", pam_ids$Sample)
pam_ids$ids <- paste0("TCGA-", pam_ids$ids)
pam_ids <- pam_ids[-c(10, 68, 74, 81,82,83),]
rownames(pam_ids) <- gsub("-", "\\.", pam_ids$ids)

kd = merge(x=kd, y=pam_ids, by ="ids")
kd$Sample = NULL
index <- c(1, 2, 3, 4)
values <- c( "Her2","Basal", "LumB", "LumA")
kd$pred <- values[match(kd$k1.cluster, index)]
pam_ids = kd
rownames(pam_ids) <- gsub("-", "\\.", pam_ids$ids)
# 51 out of 77 matches

pp <- prcomp(mg2_scale)
huh <- pp$x
pp_plot <- data.frame(PC1=pp$x[,1], PC2=pp$x[,2],
                      kclust = kd$pred,
                      subtype = kd$subtype)

ggplot(pp_plot, aes(PC1, PC2, shape = subtype, color = factor(kclust))) +
  geom_point(size = 3)


pp_plot <- pp_plot[order(pp_plot$subtype),]
mg2_order <- scale(mg2[, rownames(pp_plot)])


library(ComplexHeatmap)
library(circlize)
heat_df <- data.frame(SubType = pp_plot[,4])
pdf("newHmProt34.pdf", height = 6, width = 8) 
column_ha = HeatmapAnnotation(df = heat_df,
                              col = list(SubType=c(Basal="purple", Her2="dodgerblue2",
                                                   LumA="green3", LumB="red4")))

Heatmap(t(mg2_scale), top_annotation = column_ha,
        heatmap_legend_param = list(title = "Scaled Expression"),
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = T, cluster_columns = F,
        show_column_names = F)
dev.off()

#K-MEANS PROTEOMICS DATA ON All
library(ggplot2)
library(useful)
library(dplyr)
mg = brca_PROTEIN
mg = do.call(data.frame,lapply(mg, function(x) replace(x, is.infinite(x),NA)))
# mg = na.omit(mg)
# 6200 proteins left
mg = mg[-which(rowMeans(is.na(mg)) > 0.2), ]
# 7801 proteins left
row.names(mg) <- mg$GeneSymbol
mg[1] <- NULL
k = which(is.na(mg), arr.ind=TRUE)
mg[k] <- rowMeans(mg, na.rm=TRUE)[k[,1]]
mg = as.data.frame(t(mg))
mg_scale = na.omit(mg)
mg_scale <- as.data.frame(scale(mg))
mg_scale = na.omit(mg_scale)
k1 = kmeans(mg_scale, 4, iter.max = 100)
plot.kmeans(k1, data=mg_scale, title = 'K-means proteomics using All proteins')

kd = data.frame(k1$cluster)
kd$ids <- rownames(kd)

pam_ids <- read.table("~/downloads/ATIB M5/subtypes.txt", header=T, sep='\t')
pam_ids$ids <- gsub("\\..*", "", pam_ids$Sample)
pam_ids$ids <- paste0("TCGA.", pam_ids$ids)
pam_ids <- pam_ids[-c(10, 68, 74, 81,82,83),]
pam_ids$ids <- gsub("-", "\\.", pam_ids$ids)
rownames(pam_ids) <- gsub("-", "\\.", pam_ids$ids)

kd = merge(x=kd, y=pam_ids, by ="ids")
kd$Sample = NULL
index <- c(1, 2, 3, 4)
values <- c("LumA", "Basal", "Her2","LumB")
kd$pred <- values[match(kd$k1.cluster, index)]
pam_ids = kd
rownames(pam_ids) <- gsub("-", "\\.", pam_ids$ids)
# 41 of 77 matches

pp <- prcomp(mg_scale)
huh <- pp$x
pp_plot <- data.frame(PC1=pp$x[,1], PC2=pp$x[,2],
                      kclust = kd$pred,
                      subtype = kd$subtype)

ggplot(pp_plot, aes(PC1, PC2, shape = subtype, color = factor(kclust))) +
  geom_point(size = 3)


pp_plot <- pp_plot[order(pp_plot$subtype),]
mg = t(mg)
mg2_order <- scale(mg[, rownames(pp_plot)])

library(ComplexHeatmap)
library(circlize)
heat_df <- data.frame(pp_plot[,4])
column_ha = HeatmapAnnotation(df = heat_df)

Heatmap(mg2_order, top_annotation = column_ha,
        heatmap_legend_param = list(title = "Proteome"),
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = T, cluster_columns = F)

#GSEA
library(GSEABase)
library(gage)
colnames(brca_PROTEIN)[1]="GeneSymbol"
aa = gsub("-.*","",brca_PROTEIN$GeneSymbol)
brca_PROTEIN$GeneSymbol = aa
brca4gage <- brca_PROTEIN
brca4gage = do.call(data.frame,lapply(brca4gage, function(x) replace(x, is.infinite(x),NA)))
brca4gage = brca4gage[!duplicated(brca4gage["GeneSymbol"]),]
row.names(brca4gage) <- brca4gage$GeneSymbol
brca4gage$GeneSymbol = NULL

cnts <- as.data.frame(apply(brca4gage, 2, as.numeric))
rownames(cnts) <- rownames(brca4gage)

cnts = cnts[-which(rowMeans(is.na(cnts)) > 0.2), ]
k = which(is.na(cnts), arr.ind=TRUE)
cnts[k] <- rowMeans(cnts, na.rm=TRUE)[k[,1]]
# # 7801 proteins left
write.table(cnts, "~/downloads/ATIB M5/new_prot.csv", sep = ",", row.names=T)

# cnts = na.omit(cnts)
# # 6199 proteins left

filename="~/downloads/ATIB M5/c2.all.v6.2.symbols.gmt"
demo.gs=readList(filename)

pam_ids <- read.table("~/downloads/ATIB M5/subtypes.txt", header=T, sep='\t')

pam_ids <- pam_ids[-c(10, 68, 74, 81,82,83),]

final_results <- data.frame(p.geomean = NA, stat.mean = NA, p.val = NA, q.val = NA,
                            set.size = NA, exp1 = NA, subtype = NA, gene = NA)
subs <- unique(pam_ids$subtype)
for (k in 1:length(subs)){
  sub_of_int <- subs[k]
  
  samps <- which(pam_ids$subtype==sub_of_int)
  refs <- which(pam_ids$subtype!=sub_of_int)
  
  cnts.kegg.p <- gage(cnts, gsets = demo.gs, ref = refs,
                      samp = samps, compare ="as.group", rank.test = T)
  gg <- data.frame(cnts.kegg.p$greater)
  gg$subtype <- rep(paste0(sub_of_int), nrow(gg))
  gg$gene = rownames(gg)
  
  final_results = rbind(final_results, gg)
}
final_results<-final_results[-1,]

write.table(final_results, "~/downloads/ATIB M5/new_prot-final.csv", sep = ",", row.names=T)

Nfinal_results = na.omit(final_results)
overlap$x
olp = final_results[final_results$gene %in% overlap$x, ]
# overlap$x = paste("1", overlap$x, sep = "")
# basalolp = final_results[row.names(final_results) %in% overlap$x, ]
olp$p.geomean = NULL
olp$stat.mean = NULL
olp$p.val = NULL
olp$set.size = NULL
olp$exp1 = NULL
rownames(olp) = NULL
Her2 = olp[olp$subtype == "Her2", ]
rownames(Her2) = Her2$gene
Her2$gene = NULL
Her2$subtype = NULL
Basal = olp[olp$subtype == "Basal", ]
rownames(Basal) = Basal$gene
Basal$gene = NULL
Basal$subtype = NULL
LumA = olp[olp$subtype == "LumA", ]
rownames(LumA) = LumA$gene
LumA$gene = NULL
LumA$subtype = NULL
LumB = olp[olp$subtype == "LumB", ]
rownames(LumB) = LumB$gene
LumB$gene = NULL
LumB$subtype = NULL

mggg = merge(Basal, Her2, by="row.names")
rownames(mggg) = mggg$Row.names
mggg$Row.names = NULL
mggg = merge(mggg, LumA, by="row.names")
rownames(mggg) = mggg$Row.names
mggg$Row.names = NULL
mggg = merge(mggg, LumB, by="row.names")
rownames(mggg) = mggg$Row.names
mggg$Row.names = NULL
colnames(mggg) = c("Basal", "Her2", "LumA", "LumB")
mggg = log10(mggg)
mggg = -mggg
mggg = mggg[order(mggg$Basal),]

library(ComplexHeatmap)
library(circlize)

Heatmap(t(mggg), name = "Gage Proteomics enrichment pathways",row_title_gp = gpar(fontsize = 1),column_title_gp = gpar(fontsize = 1),
        heatmap_legend_param = list(title = "-log10 Proteome"),column_title_rot = 0,
        col = colorRamp2(c(0, 0.05, 0.1), c("blue", "white", "red")),
        cluster_rows = F, cluster_columns = T)




olp = phospho[phospho$pathway %in% overlap$x, ]
# overlap$x = paste("1", overlap$x, sep = "")
# basalolp = final_results[row.names(final_results) %in% overlap$x, ]
olp$p.geomean = NULL
olp$stat.mean = NULL
olp$p.val = NULL
olp$set.size = NULL
olp$exp1 = NULL
rownames(olp) = NULL
Her2 = olp[olp$subtype == "Her2", ]
rownames(Her2) = Her2$pathway
Her2$pathway = NULL
Her2$subtype = NULL
Basal = olp[olp$subtype == "Basal", ]
rownames(Basal) = Basal$pathway
Basal$pathway = NULL
Basal$subtype = NULL
LumA = olp[olp$subtype == "LumA", ]
rownames(LumA) = LumA$pathway
LumA$pathway = NULL
LumA$subtype = NULL
LumB = olp[olp$subtype == "LumB", ]
rownames(LumB) = LumB$pathway
LumB$pathway = NULL
LumB$subtype = NULL

mggg = merge(Basal, Her2, by="row.names")
rownames(mggg) = mggg$Row.names
mggg$Row.names = NULL
mggg = merge(mggg, LumA, by="row.names")
rownames(mggg) = mggg$Row.names
mggg$Row.names = NULL
mggg = merge(mggg, LumB, by="row.names")
rownames(mggg) = mggg$Row.names
mggg$Row.names = NULL
colnames(mggg) = c("Basal", "Her2", "LumA", "LumB")
mggg = log10(mggg)
mggg = -mggg
mggg = mggg[order(mggg$Basal),]

library(ComplexHeatmap)
library(circlize)
Heatmap(t(mggg), name = "Gage Proteomics enrichment pathways",row_title_gp = gpar(fontsize = 1),column_title_gp = gpar(fontsize = 1),
        heatmap_legend_param = list(title = "-log10 Phosphoproteome"),column_title_rot = 0,
        col = colorRamp2(c(0, 0.05, 0.1), c("blue", "white", "red")),
        cluster_rows = F, cluster_columns = T)