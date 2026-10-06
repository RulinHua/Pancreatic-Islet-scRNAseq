# SCENIC donor-level regulon analysis
# ===================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

##### 读取 HPAP AUC 文件
auc.hpap.raw <- fread(file.path(result_dir, "SCENIC/HPAP_V3_AUC.csv"),check.names = FALSE)
auc.hpap.ids <- as.character(auc.hpap.raw[[1]])
auc.hpap.temp <- as.matrix(auc.hpap.raw[, -1, with = FALSE])
storage.mode(auc.hpap.temp) <- "numeric"
rownames(auc.hpap.temp) <- auc.hpap.ids
auc.hpap <- t(auc.hpap.temp)

##### 读取 IIDP AUC 文件
auc.iidp.raw <- fread(file.path(result_dir, "SCENIC/IIDP_V3_AUC.csv"),check.names = FALSE)
auc.iidp.ids <- as.character(auc.iidp.raw[[1]])
auc.iidp.temp <- as.matrix(auc.iidp.raw[, -1, with = FALSE])
storage.mode(auc.iidp.temp) <- "numeric"
rownames(auc.iidp.temp) <- auc.iidp.ids
auc.iidp <- t(auc.iidp.temp)

##### AUC 与 Seurat 细胞匹配
meta.hpap <- beta.hpap.v3@meta.data[ colnames(auc.hpap),, drop = FALSE]
meta.iidp <- beta.iidp.v3@meta.data[ colnames(auc.iidp),, drop = FALSE]

##### 去掉细胞数过少的 donor
min.cells.per.donor <- 1
### HPAP
donor.cell.number.hpap <- table(meta.hpap$donor_accession)
keep.donors.hpap <- names(donor.cell.number.hpap[donor.cell.number.hpap >= min.cells.per.donor])
keep.cells.hpap <- rownames(meta.hpap)[meta.hpap$donor_accession %in% keep.donors.hpap]
auc.hpap <- auc.hpap[,keep.cells.hpap,drop = FALSE]
meta.hpap <- meta.hpap[keep.cells.hpap,,drop = FALSE]

### IIDP
donor.cell.number.iidp <- table(meta.iidp$donor_accession)
keep.donors.iidp <- names(donor.cell.number.iidp[donor.cell.number.iidp >= min.cells.per.donor])
keep.cells.iidp <- rownames(meta.iidp)[meta.iidp$donor_accession %in% keep.donors.iidp]
auc.iidp <- auc.iidp[,keep.cells.iidp,drop = FALSE]
meta.iidp <- meta.iidp[keep.cells.iidp,,drop = FALSE]

##### HPAP pseudobulk
counts.hpap <- GetAssayData(beta.hpap.v3,assay = "RNA",slot = "counts")
counts.hpap <- counts.hpap[,colnames(auc.hpap),drop = FALSE]
meta.hpap$library_size <- Matrix::colSums(counts.hpap) ## 增加测序深度信息

donor.factor.hpap <- factor(meta.hpap$donor_accession,levels = sort(unique(meta.hpap$donor_accession)))
aggregation.hpap <- Matrix::sparse.model.matrix(~ 0 + donor.factor.hpap)
colnames(aggregation.hpap) <- levels(donor.factor.hpap)

pseudobulk.counts.hpap <- counts.hpap %*% aggregation.hpap

### DESeq2 normalized counts
keep.genes.hpap <- Matrix::rowSums(pseudobulk.counts.hpap) > 0
pseudobulk.counts.hpap.keep <- pseudobulk.counts.hpap[keep.genes.hpap,,drop = FALSE]
coldata.hpap <- data.frame(donor_accession = colnames(pseudobulk.counts.hpap.keep),row.names = colnames(pseudobulk.counts.hpap.keep))
dds.norm.hpap <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk.counts.hpap.keep),colData = coldata.hpap,design = ~ 1)
dds.norm.hpap <- estimateSizeFactors(dds.norm.hpap)
normalized.counts.hpap <- counts(dds.norm.hpap,normalized = TRUE)
SCN9A.normalized.hpap <- normalized.counts.hpap["SCN9A",,drop = TRUE]
SCN9A.log2normalized.hpap <- log2(SCN9A.normalized.hpap + 1)
SCN9A.detect.hpap <- tapply(as.numeric(counts.hpap["SCN9A", ] > 0),donor.factor.hpap,mean)

### HPAP 按 donor 聚合 regulon AUC
auc.sum.hpap <- rowsum(t(auc.hpap),group = donor.factor.hpap)
auc.mean.hpap <- sweep(auc.sum.hpap,MARGIN = 1,STATS = as.numeric(table(donor.factor.hpap)[rownames(auc.sum.hpap)]),FUN = "/")
auc.donor.hpap <- t(auc.mean.hpap)

meta.dt.hpap <- as.data.table(meta.hpap,keep.rownames = "cell")
meta.dt.hpap[, hba1c_group := trimws(as.character(hba1c_group))]
donor.info.hpap <- meta.dt.hpap[,.(n_cells = .N,hba1c_group = unique(hba1c_group),median_library_size = median(library_size,na.rm = TRUE)),by = donor_accession]
donor.info.hpap[,SCN9A_log2normalized := as.numeric(SCN9A.log2normalized.hpap[donor_accession])]
donor.info.hpap[,SCN9A_detection_rate := as.numeric(SCN9A.detect.hpap[donor_accession])]
donor.info.hpap <- donor.info.hpap[match(colnames(auc.donor.hpap),donor_accession)]

##### IIDP pseudobulk
counts.iidp <- GetAssayData(beta.iidp.v3,assay = "RNA",slot = "counts")
counts.iidp <- counts.iidp[,colnames(auc.iidp),drop = FALSE]
meta.iidp$library_size <- Matrix::colSums(counts.iidp)
donor.factor.iidp <- factor(meta.iidp$donor_accession,levels = sort(unique(meta.iidp$donor_accession)))
aggregation.iidp <- Matrix::sparse.model.matrix(~ 0 + donor.factor.iidp)
colnames(aggregation.iidp) <- levels(donor.factor.iidp)
pseudobulk.counts.iidp <- counts.iidp %*% aggregation.iidp

### DESeq2 normalization of IIDP donor pseudobulk counts
keep.genes.iidp <- Matrix::rowSums(pseudobulk.counts.iidp) > 0
pseudobulk.counts.iidp.keep <- pseudobulk.counts.iidp[keep.genes.iidp,,drop = FALSE]
coldata.iidp <- data.frame(donor_accession = colnames(pseudobulk.counts.iidp.keep),row.names = colnames(pseudobulk.counts.iidp.keep))
dds.norm.iidp <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk.counts.iidp.keep),colData = coldata.iidp,design = ~ 1)
dds.norm.iidp <- estimateSizeFactors(dds.norm.iidp)
normalized.counts.iidp <- counts(dds.norm.iidp,normalized = TRUE)
SCN9A.normalized.iidp <- normalized.counts.iidp["SCN9A",,drop = TRUE]
SCN9A.log2normalized.iidp <- log2(SCN9A.normalized.iidp + 1)
SCN9A.detect.iidp <- tapply(as.numeric(counts.iidp["SCN9A", ] > 0),donor.factor.iidp,mean)

### IIDP 按 donor 聚合 regulon AUC
auc.sum.iidp <- rowsum(t(auc.iidp),group = donor.factor.iidp)
auc.mean.iidp <- sweep(auc.sum.iidp,MARGIN = 1,STATS = as.numeric(table(donor.factor.iidp)[rownames(auc.sum.iidp)]),FUN = "/")
auc.donor.iidp <- t(auc.mean.iidp)
meta.dt.iidp <- as.data.table(meta.iidp,keep.rownames = "cell")
meta.dt.iidp[, hba1c_group := trimws(as.character(hba1c_group))]
donor.info.iidp <- meta.dt.iidp[,.(n_cells = .N,hba1c_group = unique(hba1c_group),median_library_size = median(library_size,na.rm = TRUE)),by = donor_accession]
donor.info.iidp[,SCN9A_log2normalized := as.numeric(SCN9A.log2normalized.iidp[donor_accession])]
donor.info.iidp[,SCN9A_detection_rate := as.numeric(SCN9A.detect.iidp[donor_accession])]
donor.info.iidp <- donor.info.iidp[match(colnames(auc.donor.iidp),donor_accession)]

##### HPAP：SCN9A 连续变量模型
analysis.info.hpap <- donor.info.hpap[!is.na(hba1c_group)]
analysis.info.hpap[,hba1c_group := factor(hba1c_group,levels = c("Normal","Prediabetes","Diabetes"))]
analysis.info.hpap[,SCN9A_z := as.numeric(scale(SCN9A_log2normalized))]
analysis.info.hpap[,depth_z := as.numeric(scale(log10(median_library_size + 1)))]
auc.model.hpap <- auc.donor.hpap[,analysis.info.hpap$donor_accession,drop = FALSE]

### 去掉没有变异的 regulon
regulon.sd.hpap <- apply(auc.model.hpap,1,sd,na.rm = TRUE)
auc.model.hpap <- auc.model.hpap[is.finite(regulon.sd.hpap) & regulon.sd.hpap > 1e-12,,drop = FALSE]
auc.z.hpap <- t(scale(t(auc.model.hpap))) ## 将每个 regulon 在 donor 间标准化
auc.z.hpap <- auc.z.hpap[rowSums(is.na(auc.z.hpap)) == 0,,drop = FALSE]
design.hpap <- model.matrix(~ SCN9A_z + hba1c_group + depth_z,data = analysis.info.hpap)
if (qr(design.hpap)$rank < ncol(design.hpap)) {stop("HPAP设计矩阵不满秩，请检查HbA1c分组或协变量")}

### 运行 limma
fit.hpap <- lmFit(auc.z.hpap,design.hpap)
fit.hpap <- eBayes(fit.hpap)
result.hpap <- topTable(fit.hpap,coef = "SCN9A_z",number = Inf,sort.by = "P")
result.hpap$regulon <- rownames(result.hpap)
result.hpap <- as.data.table(result.hpap)
setnames(result.hpap,old = c("logFC","t","P.Value","adj.P.Val"),new = c("beta_HPAP","t_HPAP","P_HPAP","FDR_HPAP"))
result.hpap[,SE_HPAP := abs(beta_HPAP / t_HPAP)]
result.hpap[!is.finite(SE_HPAP),SE_HPAP := NA_real_]
fwrite(result.hpap,file.path(result_dir, "HPAP_regulon_SCN9A_continuous_limma.tsv"),sep = "\t")

##### IIDP：独立验证
analysis.info.iidp <- donor.info.iidp[!is.na(hba1c_group)]
analysis.info.iidp[,hba1c_group := factor(hba1c_group,levels = c("Normal","Prediabetes","Diabetes"))]
analysis.info.iidp[,hba1c_group := droplevels(hba1c_group)]
analysis.info.iidp[,SCN9A_z := as.numeric(scale(SCN9A_log2normalized))]
analysis.info.iidp[,depth_z := as.numeric(scale(log10(median_library_size + 1)))]

auc.model.iidp <- auc.donor.iidp[,analysis.info.iidp$donor_accession,drop = FALSE]
regulon.sd.iidp <- apply(auc.model.iidp,1,sd,na.rm = TRUE)
auc.model.iidp <- auc.model.iidp[is.finite(regulon.sd.iidp) & regulon.sd.iidp > 1e-12,,drop = FALSE]
auc.z.iidp <- t(scale(t(auc.model.iidp)))
auc.z.iidp <- auc.z.iidp[rowSums(is.na(auc.z.iidp)) == 0,,drop = FALSE]

design.iidp <- model.matrix(~ SCN9A_z + hba1c_group + depth_z,data = analysis.info.iidp)
if (qr(design.iidp)$rank < ncol(design.iidp)) {stop("IIDP设计矩阵不满秩，请检查HbA1c分组或协变量")}
fit.iidp <- lmFit(auc.z.iidp,design.iidp)
fit.iidp <- eBayes(fit.iidp)
result.iidp <- topTable(fit.iidp,coef = "SCN9A_z",number = Inf,sort.by = "P")
result.iidp$regulon <- rownames(result.iidp)
result.iidp <- as.data.table(result.iidp)
setnames(result.iidp,old = c("logFC","t","P.Value","adj.P.Val"),new = c("beta_IIDP","t_IIDP","P_IIDP","FDR_IIDP"))
result.iidp[,SE_IIDP := abs(beta_IIDP / t_IIDP)]
result.iidp[!is.finite(SE_IIDP),SE_IIDP := NA_real_]
fwrite(result.iidp,file.path(result_dir, "IIDP_regulon_SCN9A_validation_limma.tsv"),sep = "\t")

##### 合并 HPAP 和 IIDP 结果
result.combined <- merge(result.hpap[,.(regulon,beta_HPAP,SE_HPAP,t_HPAP,P_HPAP,FDR_HPAP)],
                         result.iidp[,.(regulon,beta_IIDP,SE_IIDP,t_IIDP,P_IIDP,FDR_IIDP)],
                         by = "regulon",all = TRUE)

### 判断方向是否一致
result.combined[,same_direction:=!is.na(beta_HPAP) & !is.na(beta_IIDP) & beta_HPAP * beta_IIDP > 0 ]
### 定义验证成功
result.combined[,replicated :=!is.na(FDR_HPAP) & FDR_HPAP < 0.05 & !is.na(P_IIDP) & P_IIDP < 0.05 & same_direction]

##### HPAP 与 IIDP 固定效应整合
result.combined[,meta_available :=is.finite(SE_HPAP) & SE_HPAP > 0 & is.finite(SE_IIDP) & SE_IIDP > 0]
result.combined[,beta_meta := NA_real_]
result.combined[,SE_meta := NA_real_]
result.combined[,P_meta := NA_real_]
result.combined[meta_available == TRUE,beta_meta := (beta_HPAP / SE_HPAP^2 + beta_IIDP / SE_IIDP^2) / (1 / SE_HPAP^2 + 1 / SE_IIDP^2)]
result.combined[meta_available == TRUE,SE_meta := sqrt(1 / (1 / SE_HPAP^2 + 1 / SE_IIDP^2))]
result.combined[meta_available == TRUE,P_meta := 2 * pnorm(-abs(beta_meta / SE_meta))]
result.combined[,FDR_meta := NA_real_]
result.combined[meta_available == TRUE,FDR_meta := p.adjust(P_meta,method = "BH")]

setorder(result.combined,FDR_meta,FDR_HPAP,P_IIDP)
fwrite(result.combined,file.path(result_dir, "SCN9A_regulon_combined_results.tsv"),sep = "\t")

output.xlsx <- file.path(result_dir, "SCN9A_regulon_limma_results.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "HPAP")
addWorksheet(wb, "IIDP")
addWorksheet(wb, "Combined")
writeData(wb,sheet = "HPAP",x = result.hpap)
writeData(wb,sheet = "IIDP",x = result.iidp)
writeData(wb,sheet = "Combined",x = result.combined)
saveWorkbook(wb,file = output.xlsx,overwrite = TRUE)

##### donor-level SCN9A High/Low分组和热图
SCN9A.cutoff.hpap <- median(analysis.info.hpap$SCN9A_log2normalized,na.rm = TRUE)
analysis.info.hpap[,SCN9A_group := ifelse(SCN9A_log2normalized > SCN9A.cutoff.hpap,"High","Low")]
analysis.info.hpap[,SCN9A_group := factor(SCN9A_group,levels = c("Low", "High"))]
SCN9A.cutoff.iidp <- median(analysis.info.iidp$SCN9A_log2normalized,na.rm = TRUE)
analysis.info.iidp[,SCN9A_group := ifelse(SCN9A_log2normalized > SCN9A.cutoff.iidp,"High","Low")]
analysis.info.iidp[,SCN9A_group := factor(SCN9A_group,levels = c("Low", "High"))]

### HPAP + IIDP combined/meta-supported regulons
meta.regulons <- result.combined[same_direction == TRUE & !is.na(FDR_meta) & FDR_meta < 0.05][order(FDR_meta)]
nrow(meta.regulons)
selected.regulons <- meta.regulons$regulon

### High/Low 四列热图
HPAP.Low <- rowMeans(auc.z.hpap[selected.regulons,analysis.info.hpap$donor_accession[analysis.info.hpap$SCN9A_group == "Low"],drop = FALSE])
HPAP.High <- rowMeans(auc.z.hpap[selected.regulons,analysis.info.hpap$donor_accession[analysis.info.hpap$SCN9A_group == "High"],drop = FALSE])
IIDP.Low <- rowMeans(auc.z.iidp[selected.regulons,analysis.info.iidp$donor_accession[analysis.info.iidp$SCN9A_group == "Low"],drop = FALSE])
IIDP.High <- rowMeans(auc.z.iidp[selected.regulons,analysis.info.iidp$donor_accession[analysis.info.iidp$SCN9A_group == "High"],drop = FALSE])
heatmap.matrix <- cbind(HPAP_Low = HPAP.Low,HPAP_High = HPAP.High,IIDP_Low = IIDP.Low,IIDP_High = IIDP.High)

max.abs <- max(abs(heatmap.matrix),na.rm = TRUE)
heatmap.breaks <- seq(-max.abs,max.abs,length.out = 101)
pheatmap::pheatmap(heatmap.matrix,scale = "none",cluster_rows = TRUE,cluster_cols = FALSE,
                   border_color = NA,color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   breaks = heatmap.breaks,main = "Regulon activity associated with SCN9A",
                   fontsize_row = 7,fontsize_col = 10,
                   filename =file.path(result_dir, "SCN9A_meta_supported_regulon_High_Low_heatmap.pdf"),
                   cellheight = 8,cellwidth = 20)

### 每个 donor 一列的热图
order.hpap <- order(analysis.info.hpap$SCN9A_log2normalized)
order.iidp <- order(analysis.info.iidp$SCN9A_log2normalized)
donor.heatmap.hpap <- auc.z.hpap[selected.regulons,analysis.info.hpap$donor_accession[order.hpap],drop = FALSE]
donor.heatmap.iidp <- auc.z.iidp[selected.regulons,analysis.info.iidp$donor_accession[order.iidp],drop = FALSE]
donor.heatmap <- cbind(donor.heatmap.hpap,donor.heatmap.iidp)

## donor annotation
annotation.hpap <- data.frame(Source = "HPAP",SCN9A_group = analysis.info.hpap$SCN9A_group[order.hpap],HbA1c_group = analysis.info.hpap$hba1c_group[order.hpap])
rownames(annotation.hpap) <- colnames(donor.heatmap.hpap)
annotation.iidp <- data.frame(Source = "IIDP",SCN9A_group = analysis.info.iidp$SCN9A_group[order.iidp],HbA1c_group = analysis.info.iidp$hba1c_group[order.iidp])
rownames(annotation.iidp) <- colnames(donor.heatmap.iidp)
annotation.donor <- rbind(annotation.hpap,annotation.iidp)

max.abs <- max(abs(donor.heatmap),na.rm = TRUE)
heatmap.breaks <- seq(-max.abs,max.abs,length.out = 101)

pheatmap::pheatmap(donor.heatmap,scale = "none",cluster_rows = TRUE,cluster_cols = FALSE,
                   annotation_col = annotation.donor,show_colnames = FALSE,
                   border_color = NA,color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),breaks = heatmap.breaks,
                   main = "Donor-level regulon activity ordered by SCN9A",fontsize_row = 7,fontsize_col = 8,
                   filename = file.path(result_dir, "SCN9A_meta_supported_regulon_donor_heatmap.pdf"),
                   cellheight = 8,cellwidth = 10)

### 选择HPAP前10个regulons画热图
## 选择HPAP中与SCN9A关联最显著的前10个regulons
result.hpap.plot <- copy(result.hpap)
setorder(result.hpap.plot,FDR_HPAP,P_HPAP)
top10.hpap <- head(result.hpap.plot,10)
top10.regulons.hpap <- top10.hpap$regulon
top10.hpap[,.(regulon,beta_HPAP,t_HPAP,P_HPAP,FDR_HPAP)]

## donor按照SCN9A表达量从低到高排列
order.hpap <- order(analysis.info.hpap$SCN9A_log2normalized)
donor.order.hpap <- analysis.info.hpap$donor_accession[order.hpap]

## 提取HPAP Top 10 regulons
donor.heatmap.hpap <- auc.z.hpap[top10.regulons.hpap,donor.order.hpap,drop = FALSE]

## 保持regulon按照统计显著性排序
donor.heatmap.hpap <- donor.heatmap.hpap[top10.regulons.hpap,,drop = FALSE]

## donor注释
annotation.donor.hpap <- data.frame(SCN9A_group = analysis.info.hpap$SCN9A_group[order.hpap],
                                    HbA1c_group = analysis.info.hpap$hba1c_group[order.hpap],
                                    SCN9A_expression = analysis.info.hpap$SCN9A_log2normalized[order.hpap])
rownames(annotation.donor.hpap) <- donor.order.hpap

## 注释颜色
annotation.colors.hpap <- list(SCN9A_group = c("Low" = "#4DBBD5","High" = "#E64B35"),
                               HbA1c_group = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")
)

## 对称色阶
max.abs.hpap <- max(abs(donor.heatmap.hpap),na.rm = TRUE)
heatmap.breaks.hpap <- seq(-max.abs.hpap,max.abs.hpap,length.out = 101)

## 绘制热图
pheatmap::pheatmap(donor.heatmap.hpap,scale = "none",
                   cluster_rows = T,cluster_cols = FALSE,
                   annotation_col = annotation.donor.hpap,annotation_colors = annotation.colors.hpap,
                   show_colnames = FALSE,show_rownames = TRUE,border_color = NA,
                   color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
                   breaks = heatmap.breaks.hpap,
                   main = "Top 10 SCN9A-associated regulons in HPAP",
                   fontsize_row = 9,fontsize_col = 8,width = 12,height = 4,
                   filename = file.path(result_dir, "SCN9A_HPAP_top10_regulon_donor_heatmap.pdf")
)

##### donor-level 箱线图
plot.regulons <- head(selected.regulons,10) ### 展示前 10 个 regulons
plot.data.hpap <- as.data.table(t(auc.donor.hpap[plot.regulons,analysis.info.hpap$donor_accession,drop = FALSE]),keep.rownames = "donor_accession")
plot.data.hpap <- melt(plot.data.hpap,id.vars = "donor_accession",variable.name = "regulon",value.name = "AUC")
plot.data.hpap <- merge(plot.data.hpap,analysis.info.hpap[,.(donor_accession,SCN9A_group,SCN9A_log2normalized)],by = "donor_accession")
plot.data.hpap[,cohort := "HPAP"]

plot.data.iidp <- as.data.table(t(auc.donor.iidp[plot.regulons,analysis.info.iidp$donor_accession,drop = FALSE]),keep.rownames = "donor_accession")
plot.data.iidp <- melt(plot.data.iidp,id.vars = "donor_accession",variable.name = "regulon",value.name = "AUC")
plot.data.iidp <- merge(plot.data.iidp,analysis.info.iidp[,.(donor_accession,SCN9A_group,SCN9A_log2normalized)],by = "donor_accession")
plot.data.iidp[,cohort := "IIDP"]

plot.data <- rbind(plot.data.hpap,plot.data.iidp)

p <- ggplot(plot.data,aes(x = SCN9A_group,y = AUC,fill = SCN9A_group)) +
  geom_boxplot(width = 0.6,outlier.shape = NA) +
  geom_jitter(width = 0.12,size = 1.8,alpha = 0.75) +
  facet_grid(regulon ~ cohort,scales = "free_y") +
  scale_fill_manual(values = c(Low = "#4DBBD5",High = "#E64B35")) +
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0),legend.position = "none") +
  labs(x = "Donor-level SCN9A group",y = "Mean regulon AUC")
ggsave(file.path(result_dir, "SCN9A_High_Low_regulon_boxplot.pdf"),p,width = 6,height =10)

##### 连续 SCN9A–regulon 散点图
p <- ggplot(plot.data,aes(x = SCN9A_log2normalized,y = AUC,color = SCN9A_group)) +
  geom_point(size = 2.2,alpha = 0.8) +
  geom_smooth(method = "lm",se = TRUE,color = "black",linewidth = 0.6) +
  facet_grid(regulon ~ cohort,scales = "free") +
  scale_color_manual(values = c(Low = "#4DBBD5",High = "#E64B35")) +
  theme_classic() +
  labs(x = "SCN9A pseudobulk expression (log2 normalized counts)",y = "Mean regulon AUC")
ggsave(file.path(result_dir, "SCN9A_regulon_continuous_scatter.pdf"),p,width = 9,height = max(7,length(plot.regulons) * 2))

### 每个regulon单独保存一张散点图
plot.regulons <- top10.regulons.hpap
plot.data.hpap <- as.data.table(t(auc.donor.hpap[plot.regulons,analysis.info.hpap$donor_accession,drop = FALSE]),keep.rownames = "donor_accession")
plot.data.hpap <- melt(plot.data.hpap,id.vars = "donor_accession",variable.name = "regulon",value.name = "AUC") ## 转成长数据
## 加入donor层面的SCN9A信息
plot.data.hpap <- merge(plot.data.hpap,analysis.info.hpap[,.(donor_accession,SCN9A_group,SCN9A_log2normalized)],by = "donor_accession",all.x = TRUE,sort = FALSE)
plot.data.hpap[,regulon := factor(as.character(regulon),levels = plot.regulons)]
plot.data.hpap[,SCN9A_group := factor(trimws(as.character(SCN9A_group)),levels = c("Low","High"))]

scatter.output.dir <- file.path(result_dir, "SCN9A_regulon_continuous_scatter_HPAP")
for(regulon.name in top10.regulons.hpap){
  current.data <- plot.data.hpap[regulon == regulon.name]
  
  ## 生成适合文件名的regulon名称
  safe.regulon.name <- gsub("[^A-Za-z0-9_-]+","_",regulon.name)
  
  ## 绘图
  p <- ggplot(current.data,aes(x = SCN9A_log2normalized,y = AUC,colour = SCN9A_group)) +
    geom_point(size = 2.8,alpha = 0.85) +
    geom_smooth(aes(group = 1),method = "lm",formula = y ~ x,se = TRUE,colour = "black",fill = "grey75",linewidth = 0.8) +
    scale_colour_manual(values = c("Low" = "#4DBBD5","High" = "#E64B35"),drop = FALSE) +
    labs(title = paste0(regulon.name," activity versus SCN9A expression"),
         subtitle = "HPAP V3 Beta cells; each point represents one donor",
         x = "SCN9A expression\n(log2 DESeq2-normalized pseudobulk counts + 1)",
         y = paste0(regulon.name," mean regulon AUC"),
         colour = "SCN9A group") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
          plot.subtitle = element_text(hjust = 0.5,size = 11),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 11,colour = "black"),
          legend.position = "right")
  
  ggsave(file.path(scatter.output.dir,paste0("SCN9A_HPAP_",safe.regulon.name,"_continuous_scatter.pdf")),
         plot = p,width = 6.5,height = 5.5)
}

##### HPAP Top 10 regulons：SCN9A UMAP + AUC UMAP + Violin
### 计算HPAP V3 β细胞内部UMAP
auc.hpap.plot <- auc.hpap[top10.regulons.hpap,colnames(beta.hpap.v3),drop = FALSE]
beta.hpap.v3 <- FindVariableFeatures(beta.hpap.v3,selection.method = "vst",nfeatures = 3000,verbose = FALSE)
beta.hpap.v3 <- ScaleData(beta.hpap.v3)
beta.hpap.v3 <- RunPCA(beta.hpap.v3,npcs = 50)
ElbowPlot(beta.hpap.v3,ndims = 50)
beta.hpap.v3 <- RunUMAP(beta.hpap.v3,reduction = "pca",dims = 1:50,reduction.name = "umap.hpap.beta",reduction.key = "HPAPBetaUMAP_",seed.use = 123,verbose = FALSE)

### 准备SCN9A表达量
## 提取Seurat LogNormalize后的SCN9A表达量
SCN9A.expression <- GetAssayData(beta.hpap.v3,assay = "RNA",slot = "data")["SCN9A",,drop = TRUE]
SCN9A.expression <- as.numeric(SCN9A.expression)
names(SCN9A.expression) <- colnames(beta.hpap.v3)

## 计算SCN9A阳性细胞的99%分位数
SCN9A.positive <- SCN9A.expression[SCN9A.expression > 0]
SCN9A.upper <- as.numeric(quantile(SCN9A.positive,probs = 0.99,na.rm = TRUE))

## 仅限制颜色显示范围，不修改原表达矩阵
SCN9A.display <- pmin(SCN9A.expression,SCN9A.upper)
names(SCN9A.display) <- colnames(beta.hpap.v3)
beta.hpap.v3$SCN9A_display <- SCN9A.display

### SCN9A表达UMAP
panel.output.dir <- file.path(result_dir, "SCN9A_HPAP_top10_regulon_UMAP_violin")
p <- FeaturePlot(beta.hpap.v3,features = "SCN9A_display",reduction = "umap.hpap.beta",
                 cols = c("#D9D9D9","#E64B35"),order = TRUE,pt.size = 0.12,raster = F) +
  labs(title = "SCN9A expression in HPAP V3 Beta cells",x = "UMAP 1",y = "UMAP 2",colour = "SCN9A\nexpression") +
  theme_classic() +
  theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10,colour = "black"))
ggsave(filename = file.path(panel.output.dir,"SCN9A_expression_UMAP_HPAP_V3_Beta.pdf"),p[[1]],width = 6,height = 5)

for(regulon.name in top10.regulons.hpap){
  ## 提取TF名称
  tf.name <- sub("\\(.*$","",regulon.name)
  tf.name <- sub("_extended$","", tf.name)
  
  ## 生成安全文件名
  safe.regulon.name <- gsub("[^A-Za-z0-9_-]+","_",regulon.name)
  
  ## 提取当前regulon的细胞水平AUC
  regulon.auc <- as.numeric(auc.hpap.plot[regulon.name,colnames(beta.hpap.v3)])
  names(regulon.auc) <- colnames(beta.hpap.v3)
  
  ## 限制UMAP颜色显示范围
  auc.lower <- as.numeric(quantile(regulon.auc,probs = 0.01,na.rm = TRUE))
  auc.upper <- as.numeric(quantile(regulon.auc,probs = 0.99,na.rm = TRUE))
  regulon.auc.display <- pmin(pmax(regulon.auc,auc.lower),auc.upper)
  names(regulon.auc.display) <- colnames(beta.hpap.v3)
  
  ## 将AUC加入Seurat metadata
  auc.column <- paste0("RegulonAUC_",safe.regulon.name)
  beta.hpap.v3[[auc.column]] <- regulon.auc.display
  
  ## Regulon AUC UMAP
  p.regulon <- FeaturePlot(beta.hpap.v3,features = auc.column,reduction = "umap.hpap.beta",
                           cols = c("#D9D9D9","#3B4CC0"),order = TRUE,pt.size = 0.12,raster = F) +
    labs(title = paste0(tf.name," regulon activity"),subtitle = "HPAP V3 Beta cells",
         x = "UMAP 1",y = "UMAP 2",colour = paste0(tf.name,"\nregulon AUC")) +
    theme_classic() +
    theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
          plot.subtitle = element_text(hjust = 0.5,size = 11),
          axis.title = element_text(size = 12),axis.text = element_text(size = 10,colour = "black"))
  
  ## HbA1c分组小提琴图
  violin.data <- data.table(cell = colnames(beta.hpap.v3),
                            donor_accession = beta.hpap.v3$donor_accession,
                            hba1c_group = beta.hpap.v3$hba1c_group,AUC = regulon.auc)
  violin.data <- violin.data[!is.na(hba1c_group) & is.finite(AUC)]
  
  ## 每个donor的平均regulon AUC
  donor.violin.data <- violin.data[,.(Mean_AUC = mean(AUC,na.rm = TRUE)),by = .(donor_accession,hba1c_group)]
  donor.violin.data[,hba1c_group := factor(as.character(hba1c_group),levels = c("Normal","Prediabetes","Diabetes"))]
  
  ## donor层面总体检验
  kw.result <- kruskal.test(Mean_AUC ~ hba1c_group,data = donor.violin.data)
  kw.label <- paste0("Donor-level overall Kruskal-Wallis p ",format.pval(kw.result$p.value,digits = 3,eps = 0.001))
  
  p.violin <- ggplot(violin.data,aes(x = hba1c_group,y = AUC,fill = hba1c_group)) +
    geom_violin(width = 0.82,trim = TRUE,scale = "width",alpha = 0.65,colour = "black",linewidth = 0.4) +
    geom_boxplot(width = 0.12,outlier.shape = NA,fill = "white",colour = "black",linewidth = 0.45) +
    geom_point(data = donor.violin.data,aes(x = hba1c_group,y = Mean_AUC),inherit.aes = FALSE,
               position = position_jitter(width = 0.10,height = 0),
               shape = 21,size = 2.4,stroke = 0.5,fill = "white",colour = "black",alpha = 0.9) +
    scale_fill_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35"),drop = FALSE) +
    labs(title = paste0(tf.name," regulon activity across HbA1c groups"),
         subtitle = paste0("Cell-level distribution; white points represent donor means"),
         caption = kw.label,x = NULL,y = paste0(tf.name," regulon AUC")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
          plot.subtitle = element_text(hjust = 0.5,size = 10),
          plot.caption = element_text(hjust = 0.5,size = 9),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 10,colour = "black"),
          legend.position = "none")
  
  ## 合并为两联图
  p.panel <- cowplot::plot_grid(p.regulon[[1]],p.violin,labels = c("A","B"),nrow = 1,rel_widths = c(1,1.05),align = "h",axis = "tb")
  ggsave(file.path(panel.output.dir,paste0(safe.regulon.name,"_regulon_AUC_UMAP_violin.pdf")),p.panel,width = 12,height = 5.5)
}

##### 查找直接调控 SCN9A 的候选 TF
