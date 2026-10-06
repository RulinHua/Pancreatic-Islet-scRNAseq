# CytoTRACE2 beta-cell potency analysis and replication
# =====================================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

########## CytoTRACE2
######## HPAP V3主分析
min.cells.per.donor <- 1 # donor至少需要多少个β细胞
beta.hpap.v3 <- sr[,sr$Cell_Type=="Beta" & sr$source == "HPAP" & sr$chemistry == "V3"]
beta.hpap.v3.ct2 <- CytoTRACE2::cytotrace2(beta.hpap.v3,is_seurat = TRUE,slot_type = "counts",
                                           species = "human",ncores = 20,parallelize_models = TRUE,
                                           parallelize_smoothing = TRUE,batch_size = 10000,smooth_batch_size = 1000)
ct2.meta <- beta.hpap.v3.ct2@meta.data[,c("CytoTRACE2_Score","CytoTRACE2_Potency","CytoTRACE2_Relative","preKNN_CytoTRACE2_Score","preKNN_CytoTRACE2_Potency"),drop = FALSE]
saveRDS(ct2.meta,file.path(result_dir, "CytoTRACE2/beta_hpap_v3_CytoTRACE2_metadata.rds"))
beta.hpap.v3 <- AddMetaData(beta.hpap.v3,metadata = ct2.meta[colnames(beta.hpap.v3),])
p <- DimPlot(beta.hpap.v3,reduction = "umap",group.by = "CytoTRACE2_Potency",pt.size = 0.1) +
  labs(title = "HPAP V3: predicted potency") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"))
ggsave(file.path(result_dir, "CytoTRACE2/HPAP_V3_Beta_CytoTRACE2_Potency_UMAP.pdf"),p[[1]],width = 8,height = 6)
p <- FeaturePlot(beta.hpap.v3,features = "CytoTRACE2_Score",reduction = "umap",pt.size = 0.1) +
  labs(title = "HPAP V3 Beta cells: CytoTRACE2 potency score") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"))
ggsave(file.path(result_dir, "CytoTRACE2/HPAP_V3_Beta_CytoTRACE2_Score_UMAP.pdf"),p[[1]],width = 8,height = 6)

###### 提取metadata
meta.hpap <- as.data.table(beta.hpap.v3@meta.data,keep.rownames = "cell")
meta.hpap[,donor_accession := trimws(as.character(donor_accession))]
meta.hpap[,hba1c_group := as.character(hba1c_group)]

##### 计算每个donor的β细胞数
donor.cell.number.hpap <- meta.hpap[!is.na(hba1c_group),.(n_cells = .N),by = .(donor_accession,hba1c_group)]
setorder(donor.cell.number.hpap,hba1c_group,donor_accession)
donor.cell.number.hpap

##### 筛选至少min.cells.per.donor个β细胞的donor
eligible.donors.hpap <- donor.cell.number.hpap[n_cells >= min.cells.per.donor,donor_accession]
retention.hpap <- donor.cell.number.hpap[,.(total_donors = .N,retained_donors = sum(donor_accession %in% eligible.donors.hpap)),by = hba1c_group]
retention.hpap

##### 保留合格donor对应的细胞
meta.hpap.keep <- meta.hpap[donor_accession %in% eligible.donors.hpap]
keep.cells.hpap <- colnames(beta.hpap.v3)[colnames(beta.hpap.v3) %in% meta.hpap.keep$cell]
meta.hpap.keep <- meta.hpap.keep[match(keep.cells.hpap,cell)]

##### 提取raw counts
raw.counts.hpap.keep <- beta.hpap.v3@assays$RNA@counts[,keep.cells.hpap,drop = FALSE]

##### 按donor聚合pseudobulk counts
donor.factor.hpap <- factor(meta.hpap.keep$donor_accession)
aggregation.matrix.hpap <- Matrix::sparse.model.matrix(~ 0 + donor.factor.hpap)
colnames(aggregation.matrix.hpap) <- levels(donor.factor.hpap)

##### 计算pseudobulk counts
pseudobulk.counts.hpap <- raw.counts.hpap.keep %*% aggregation.matrix.hpap

##### 使用DESeq2标准化pseudobulk counts
## 建立metadata
pseudobulk.metadata.hpap <- data.frame(donor_accession =colnames(pseudobulk.counts.hpap),row.names =colnames(pseudobulk.counts.hpap))
## 建立DESeq2对象
dds.hpap.normalization <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk.counts.hpap),colData = pseudobulk.metadata.hpap,design = ~ 1)
## 估计size factor
dds.hpap.normalization <- DESeq2::estimateSizeFactors(dds.hpap.normalization)
## 提取标准化counts
normalized.counts.hpap <- DESeq2::counts(dds.hpap.normalization,normalized = TRUE)

##### 提取每个donor的SCN9A表达
scn9a.hpap <- data.table(donor_accession = colnames(normalized.counts.hpap),
                         SCN9A_raw_count = as.numeric(pseudobulk.counts.hpap["SCN9A",]),
                         SCN9A_normalized =as.numeric(normalized.counts.hpap["SCN9A",]))
scn9a.hpap[,SCN9A_log2normalized := log2(SCN9A_normalized + 1)]

##### 计算每个donor的CytoTRACE 2分数
donor.score.hpap <- meta.hpap.keep[,.(CytoTRACE2_median = median(CytoTRACE2_Score,na.rm = TRUE),
                                      CytoTRACE2_mean = mean(CytoTRACE2_Score,na.rm = TRUE),
                                      preKNN_median = median(preKNN_CytoTRACE2_Score,na.rm = TRUE),
                                      preKNN_mean = mean(preKNN_CytoTRACE2_Score,na.rm = TRUE),
                                      n_cells = .N),
                                   by = .(donor_accession,hba1c_group)]

##### 合并SCN9A与CytoTRACE 2
donor.analysis.hpap <- merge(donor.score.hpap,scn9a.hpap,by = "donor_accession",all = FALSE,sort = FALSE)
donor.analysis.hpap[,hba1c_group := factor(hba1c_group,levels = c("Normal","Prediabetes","Diabetes"))]

##### 按donor-level SCN9A中位数分组
scn9a.median.hpap <- median(donor.analysis.hpap$SCN9A_log2normalized,na.rm = TRUE)
donor.analysis.hpap[,SCN9A_group := ifelse(SCN9A_log2normalized > scn9a.median.hpap,"SCN9A_High","SCN9A_Low")]
donor.analysis.hpap[,SCN9A_group := factor(SCN9A_group,levels = c("SCN9A_Low","SCN9A_High"))]

##### 连续相关性分析
correlation.hpap <- cor.test(donor.analysis.hpap$SCN9A_log2normalized,donor.analysis.hpap$CytoTRACE2_median,method = "spearman",exact = FALSE)
correlation.hpap

##### preKNN敏感性分析
correlation.hpap.preKNN <- cor.test(donor.analysis.hpap$SCN9A_log2normalized,donor.analysis.hpap$preKNN_median,method = "spearman",exact = FALSE)
correlation.hpap.preKNN

##### SCN9A High/Low检验
wilcox.hpap <- wilcox.test(CytoTRACE2_median ~ SCN9A_group,data = donor.analysis.hpap,exact = FALSE,alternative = "two.sided")
wilcox.hpap

##### 校正HbA1c分组
model.hpap <- lm(scale(CytoTRACE2_median) ~ scale(SCN9A_log2normalized) + hba1c_group,data = donor.analysis.hpap)
summary(model.hpap) ## 如果SCN9A的回归系数为负，表示校正HbA1c分组后，SCN9A越低仍与较高CytoTRACE 2相关。

##### donor-level相关性散点图
rho.hpap <- unname(correlation.hpap$estimate)
pvalue.hpap <- ifelse(correlation.hpap$p.value < 0.001,format(correlation.hpap$p.value,scientific = TRUE,digits = 3),sprintf("%.3f", correlation.hpap$p.value))
p <- ggplot(donor.analysis.hpap,aes(x = SCN9A_log2normalized,y = CytoTRACE2_median,colour = hba1c_group)) +
  geom_point(size = 3,alpha = 0.85) +
  geom_smooth(aes(group = 1),method = "lm",formula = y ~ x,se = TRUE,colour = "black",linewidth = 0.7) +
  scale_colour_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")) +
  labs(title = "HPAP V3: SCN9A and β-cell potency",
       subtitle = paste0("Donor-level Spearman ρ = ",round(rho.hpap, 3),"; P = ",pvalue.hpap),
       x = paste0("β-cell SCN9A expression\n","log2(DESeq2 normalized count + 1)"),
       y = "Median CytoTRACE 2 score",colour = "HbA1c group") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle =element_text(hjust = 0.5),axis.text = element_text(colour = "black"))
ggsave(file.path(result_dir, "CytoTRACE2/HPAP_V3_SCN9A_CytoTRACE2_correlation.pdf"),p,width = 6,height = 5,device = cairo_pdf)

##### donor-level High/Low箱线图
pvalue.hpap <- if (wilcox.hpap$p.value < 0.001) {
  format(wilcox.hpap$p.value,scientific = TRUE,digits = 3)} else {
    sprintf("%.3f", wilcox.hpap$p.value)
  }
p <- ggplot(donor.analysis.hpap,aes(x = SCN9A_group,y = CytoTRACE2_median,fill = SCN9A_group)) +
  geom_boxplot(width = 0.55,outlier.shape = NA,alpha = 0.75) +
  geom_jitter(aes(colour = hba1c_group),width = 0.1,size = 2.5,alpha = 0.85) +
  scale_fill_manual(values = c("SCN9A_Low" = "#4DBBD5","SCN9A_High" = "#E64B35")) +
  scale_colour_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")) +
  labs(title = "HPAP V3: CytoTRACE 2 score",subtitle = paste0("Donor-level Wilcoxon P = ",pvalue.hpap),x = NULL,y = "Median CytoTRACE 2 score",colour = "HbA1c group") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle =element_text(hjust = 0.5),axis.text.x = element_text(colour = "black",face = "bold"))
ggsave(file.path(result_dir, "CytoTRACE2/HPAP_V3_SCN9A_High_Low_donor_boxplot.pdf"),p,width = 5.5,height = 5)

######## IIDP V3验证分析
beta.iidp.v3 <- sr[,sr$Cell_Type=="Beta" & sr$source == "IIDP" & sr$chemistry == "V3"]
beta.iidp.v3.ct2 <- CytoTRACE2::cytotrace2(beta.iidp.v3,is_seurat = TRUE,slot_type = "counts",
                                           species = "human",ncores = 20,parallelize_models = TRUE,
                                           parallelize_smoothing = TRUE,batch_size = 10000,smooth_batch_size = 1000)
ct2.meta <- beta.iidp.v3.ct2@meta.data[,c("CytoTRACE2_Score","CytoTRACE2_Potency","CytoTRACE2_Relative","preKNN_CytoTRACE2_Score","preKNN_CytoTRACE2_Potency"),drop = FALSE]
saveRDS(ct2.meta,file.path(result_dir, "CytoTRACE2/beta_iidp_v3_CytoTRACE2_metadata.rds"))
beta.iidp.v3 <- AddMetaData(beta.iidp.v3,metadata = ct2.meta[colnames(beta.iidp.v3),])
p <- DimPlot(beta.iidp.v3,reduction = "umap",group.by = "CytoTRACE2_Potency",pt.size = 0.1) +
  labs(title = "IIDP V3: predicted potency") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"))
ggsave(file.path(result_dir, "CytoTRACE2/IIDP_V3_Beta_CytoTRACE2_Potency_UMAP.pdf"),p[[1]],width = 8,height = 6)
p <- FeaturePlot(beta.iidp.v3,features = "CytoTRACE2_Score",reduction = "umap",pt.size = 0.1) +
  labs(title = "IIDP V3 Beta cells: CytoTRACE2 potency score") +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"))
ggsave(file.path(result_dir, "CytoTRACE2/IIDP_V3_Beta_CytoTRACE2_Score_UMAP.pdf"),p[[1]],width = 8,height = 6)

###### 提取metadata
meta.iidp <- as.data.table(beta.iidp.v3@meta.data,keep.rownames = "cell")
meta.iidp[,donor_accession := trimws(as.character(donor_accession))]
meta.iidp[,hba1c_group := as.character(hba1c_group)]

##### 计算每个donor的β细胞数
donor.cell.number.iidp <- meta.iidp[!is.na(hba1c_group),.(n_cells = .N),by = .(donor_accession,hba1c_group)]
setorder(donor.cell.number.iidp,hba1c_group,donor_accession)
donor.cell.number.iidp

##### 筛选至少min.cells.per.donor个β细胞的donor
eligible.donors.iidp <- donor.cell.number.iidp[n_cells >= min.cells.per.donor,donor_accession]
retention.iidp <- donor.cell.number.iidp[,.(total_donors = .N,retained_donors = sum(donor_accession %in% eligible.donors.iidp)),by = hba1c_group]
retention.iidp

##### 保留合格donor对应的细胞
meta.iidp.keep <- meta.iidp[donor_accession %in% eligible.donors.iidp]
keep.cells.iidp <- colnames(beta.iidp.v3)[colnames(beta.iidp.v3) %in% meta.iidp.keep$cell]
meta.iidp.keep <- meta.iidp.keep[match(keep.cells.iidp,cell)]

##### 提取raw counts
raw.counts.iidp.keep <- beta.iidp.v3@assays$RNA@counts[,keep.cells.iidp,drop = FALSE]

##### 按donor聚合pseudobulk counts
donor.factor.iidp <- factor(meta.iidp.keep$donor_accession)
aggregation.matrix.iidp <- Matrix::sparse.model.matrix(~ 0 + donor.factor.iidp)
colnames(aggregation.matrix.iidp) <- levels(donor.factor.iidp)

##### 计算pseudobulk counts
pseudobulk.counts.iidp <- raw.counts.iidp.keep %*% aggregation.matrix.iidp

##### 使用DESeq2标准化pseudobulk counts
## 建立metadata
pseudobulk.metadata.iidp <- data.frame(donor_accession =colnames(pseudobulk.counts.iidp),row.names =colnames(pseudobulk.counts.iidp))
## 建立DESeq2对象
dds.iidp.normalization <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk.counts.iidp),colData = pseudobulk.metadata.iidp,design = ~ 1)
## 估计size factor
dds.iidp.normalization <- DESeq2::estimateSizeFactors(dds.iidp.normalization)
## 提取标准化counts
normalized.counts.iidp <- DESeq2::counts(dds.iidp.normalization,normalized = TRUE)

##### 提取每个donor的SCN9A表达
scn9a.iidp <- data.table(donor_accession = colnames(normalized.counts.iidp),
                         SCN9A_raw_count = as.numeric(pseudobulk.counts.iidp["SCN9A",]),
                         SCN9A_normalized =as.numeric(normalized.counts.iidp["SCN9A",]))
scn9a.iidp[,SCN9A_log2normalized := log2(SCN9A_normalized + 1)]

##### 计算每个donor的CytoTRACE 2分数
donor.score.iidp <- meta.iidp.keep[,.(CytoTRACE2_median = median(CytoTRACE2_Score,na.rm = TRUE),
                                      CytoTRACE2_mean = mean(CytoTRACE2_Score,na.rm = TRUE),
                                      preKNN_median = median(preKNN_CytoTRACE2_Score,na.rm = TRUE),
                                      preKNN_mean = mean(preKNN_CytoTRACE2_Score,na.rm = TRUE),
                                      n_cells = .N),
                                   by = .(donor_accession,hba1c_group)]

##### 合并SCN9A与CytoTRACE 2
donor.analysis.iidp <- merge(donor.score.iidp,scn9a.iidp,by = "donor_accession",all = FALSE,sort = FALSE)
donor.analysis.iidp[,hba1c_group := factor(hba1c_group,levels = c("Normal","Prediabetes","Diabetes"))]

##### 按donor-level SCN9A中位数分组
scn9a.median.iidp <- median(donor.analysis.iidp$SCN9A_log2normalized,na.rm = TRUE)
donor.analysis.iidp[,SCN9A_group := ifelse(SCN9A_log2normalized > scn9a.median.iidp,"SCN9A_High","SCN9A_Low")]
donor.analysis.iidp[,SCN9A_group := factor(SCN9A_group,levels = c("SCN9A_Low","SCN9A_High"))]

##### 连续相关性分析
correlation.iidp <- cor.test(donor.analysis.iidp$SCN9A_log2normalized,donor.analysis.iidp$CytoTRACE2_median,method = "spearman",exact = FALSE)
correlation.iidp

##### preKNN敏感性分析
correlation.iidp.preKNN <- cor.test(donor.analysis.iidp$SCN9A_log2normalized,donor.analysis.iidp$preKNN_median,method = "spearman",exact = FALSE)
correlation.iidp.preKNN

##### SCN9A High/Low检验
wilcox.iidp <- wilcox.test(CytoTRACE2_median ~ SCN9A_group,data = donor.analysis.iidp,exact = FALSE,alternative = "two.sided")
wilcox.iidp

##### 校正HbA1c分组
model.iidp <- lm(scale(CytoTRACE2_median) ~ scale(SCN9A_log2normalized) + hba1c_group,data = donor.analysis.iidp)
summary(model.iidp) ## 如果SCN9A的回归系数为负，表示校正HbA1c分组后，SCN9A越低仍与较高CytoTRACE 2相关。

##### donor-level相关性散点图
rho.iidp <- unname(correlation.iidp$estimate)
pvalue.iidp <- ifelse(correlation.iidp$p.value < 0.001,format(correlation.iidp$p.value,scientific = TRUE,digits = 3),sprintf("%.3f", correlation.iidp$p.value))
p <- ggplot(donor.analysis.iidp,aes(x = SCN9A_log2normalized,y = CytoTRACE2_median,colour = hba1c_group)) +
  geom_point(size = 3,alpha = 0.85) +
  geom_smooth(aes(group = 1),method = "lm",formula = y ~ x,se = TRUE,colour = "black",linewidth = 0.7) +
  scale_colour_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")) +
  labs(title = "IIDP V3: SCN9A and β-cell potency",
       subtitle = paste0("Donor-level Spearman ρ = ",round(rho.iidp, 3),"; P = ",pvalue.iidp),
       x = paste0("β-cell SCN9A expression\n","log2(DESeq2 normalized count + 1)"),
       y = "Median CytoTRACE 2 score",colour = "HbA1c group") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle =element_text(hjust = 0.5),axis.text = element_text(colour = "black"))
ggsave(file.path(result_dir, "CytoTRACE2/IIDP_V3_SCN9A_CytoTRACE2_correlation.pdf"),p,width = 6,height = 5,device = cairo_pdf)

##### donor-level High/Low箱线图
pvalue.iidp <- if (wilcox.iidp$p.value < 0.001) {
  format(wilcox.iidp$p.value,scientific = TRUE,digits = 3)} else {
    sprintf("%.3f", wilcox.iidp$p.value)
  }
p <- ggplot(donor.analysis.iidp,aes(x = SCN9A_group,y = CytoTRACE2_median,fill = SCN9A_group)) +
  geom_boxplot(width = 0.55,outlier.shape = NA,alpha = 0.75) +
  geom_jitter(aes(colour = hba1c_group),width = 0.1,size = 2.5,alpha = 0.85) +
  scale_fill_manual(values = c("SCN9A_Low" = "#4DBBD5","SCN9A_High" = "#E64B35")) +
  scale_colour_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")) +
  labs(title = "IIDP V3: CytoTRACE 2 score",subtitle = paste0("Donor-level Wilcoxon P = ",pvalue.iidp),x = NULL,y = "Median CytoTRACE 2 score",colour = "HbA1c group") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle =element_text(hjust = 0.5),axis.text.x = element_text(colour = "black",face = "bold"))
ggsave(file.path(result_dir, "CytoTRACE2/IIDP_V3_SCN9A_High_Low_donor_boxplot.pdf"),p,width = 5.5,height = 5)

