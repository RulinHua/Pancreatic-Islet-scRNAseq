# SCN9A-centered donor-level beta-cell association analysis
# =========================================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

##########
meta.beta <- as.data.table(sr@meta.data,keep.rownames = "Cell")
meta.beta <- meta.beta[Cell_Type == "Beta" & diabetes_group %in% c("T1D","T2D")]
cell.number <- meta.beta[,.(n_beta_cells = .N),by = donor_accession]
minimum.cells <- 20 ## 每个donor至少保留n个β细胞
cell.number <- cell.number[n_beta_cells >= minimum.cells]
meta.beta <- meta.beta[donor_accession %in% cell.number$donor_accession]
table(meta.beta$diabetes_group)

## 建立donor-level metadata和PB_ID
donor.info <- unique(meta.beta[,.(source,chemistry,sex,age=`age_(years)`,diabetes_group,hba1c_group,donor_accession)])
donor.info <- merge(donor.info,cell.number,by = "donor_accession",all.x = TRUE,sort = FALSE)
setorder(donor.info,diabetes_group,donor_accession)
donor.info[,PB_ID := sprintf("PB%04d",seq_len(.N))]
meta.beta <- merge(meta.beta,donor.info[,.(donor_accession,PB_ID)],by = "donor_accession",all = FALSE,sort = FALSE)

## 保证细胞顺序与表达矩阵完全一致
cells.use <- colnames(sr)[colnames(sr) %in% meta.beta$Cell]
meta.beta <- meta.beta[match(cells.use,Cell)]
raw.counts <- sr@assays$RNA@counts[,cells.use,drop = FALSE]

## 计算每个donor的SCN9A阳性β细胞比例
meta.beta$SCN9A_detected <- as.numeric(raw.counts["SCN9A",,drop = TRUE]) > 0
scn9a.detection <- meta.beta[,.(SCN9A_positive_cells =sum(SCN9A_detected),SCN9A_percent_detected =mean(SCN9A_detected) * 100),by = PB_ID]
summary(scn9a.detection$SCN9A_percent_detected)

## 构建pseudobulk原始count矩阵
pb.factor <- factor(meta.beta$PB_ID,levels = donor.info$PB_ID)
aggregation.matrix <- Matrix::sparse.model.matrix(~ 0 + pb.factor)
colnames(aggregation.matrix) <- levels(pb.factor)
pseudobulk.counts <- raw.counts %*% aggregation.matrix

## 合并pseudobulk metadata
pseudobulk.metadata <- merge(donor.info,scn9a.detection,by = "PB_ID",all.x = TRUE,sort = FALSE)
pseudobulk.metadata <- pseudobulk.metadata[match(colnames(pseudobulk.counts),PB_ID)]
pseudobulk.metadata <- as.data.frame(pseudobulk.metadata)
rownames(pseudobulk.metadata) <- pseudobulk.metadata$PB_ID
stopifnot(identical(colnames(pseudobulk.counts),rownames(pseudobulk.metadata)))

## 检查是否有模型协变量缺失的donor
pseudobulk.metadata$source <- factor(pseudobulk.metadata$source)
pseudobulk.metadata$chemistry <- factor(pseudobulk.metadata$chemistry)
pseudobulk.metadata$sex <- factor(pseudobulk.metadata$sex)
pseudobulk.metadata$diabetes_group <- factor(pseudobulk.metadata$diabetes_group)
pseudobulk.metadata$age <- as.numeric(pseudobulk.metadata$age)
model.columns <- c("source","chemistry","sex","age","diabetes_group","SCN9A_percent_detected")
complete.use <- complete.cases(pseudobulk.metadata[,model.columns,drop = FALSE])
all(complete.use)
pseudobulk.counts <- as.matrix(pseudobulk.counts)
storage.mode(pseudobulk.counts) <- "integer"

## 计算SCN9A pseudobulk标准化表达
dds.sizefactor <- DESeqDataSetFromMatrix(countData = pseudobulk.counts,colData = pseudobulk.metadata,design = ~ 1)
dds.sizefactor <- estimateSizeFactors(dds.sizefactor,type = "ratio")
normalized.counts <- counts(dds.sizefactor,normalized = TRUE)
pseudobulk.metadata$SCN9A_PB_raw_count <- pseudobulk.counts["SCN9A",rownames(pseudobulk.metadata)]
pseudobulk.metadata$SCN9A_PB_normalized_count <- normalized.counts["SCN9A",rownames(pseudobulk.metadata)]
pseudobulk.metadata$SCN9A_PB_log2 <- log2(pseudobulk.metadata$SCN9A_PB_normalized_count + 1)

summary(pseudobulk.metadata$SCN9A_PB_log2)
cor.test(pseudobulk.metadata$SCN9A_PB_log2,pseudobulk.metadata$SCN9A_percent_detected,method = "spearman")

pseudobulk.metadata$SCN9A_PB_scaled <- as.numeric(scale(pseudobulk.metadata$SCN9A_PB_log2))
pseudobulk.metadata$SCN9A_detected_scaled <- as.numeric(scale(pseudobulk.metadata$SCN9A_percent_detected))
pseudobulk.metadata$age_scaled <- as.numeric(scale(pseudobulk.metadata$age))

## 建立主要模型
design.formula.primary <- reformulate(c("source","chemistry","sex","age_scaled","diabetes_group","SCN9A_PB_scaled"))
design.formula.primary

## 检查设计矩阵是否满秩
design.matrix.primary <- model.matrix(design.formula.primary,data = pseudobulk.metadata)
design.rank <- qr(design.matrix.primary)$rank
design.ncol <- ncol(design.matrix.primary)
if (design.rank < design.ncol) {
  dependent.columns <- colnames(design.matrix.primary)[qr(design.matrix.primary)$pivot[(design.rank + 1):design.ncol]]
  print(dependent.columns)
  stop(paste0("The design matrix is rank-deficient. ","Dependent columns: ",paste(dependent.columns,collapse = ", ")))
}

## 运行主要DESeq2连续变量分析
## 常规低表达预筛选
## 预设目标基因强制保留
target.genes <- c("INS", "PDX1", "MAFA", "NKX6-1","DDIT3", "ATF4", "TXNIP")
keep.gene <- rowSums(pseudobulk.counts) >= 10 | rownames(pseudobulk.counts) %in% c(target.genes,"SCN9A")
dds.primary <- DESeqDataSetFromMatrix(countData = pseudobulk.counts[keep.gene,,drop = FALSE],colData = pseudobulk.metadata,design = design.formula.primary)
sizeFactors(dds.primary) <- sizeFactors(dds.sizefactor)[colnames(dds.primary)] ## 使用前面已经估计的size factor
dds.primary <- DESeq(dds.primary,quiet = TRUE)
resultsNames(dds.primary)
res.primary <- results(dds.primary,name = "SCN9A_PB_scaled",alpha = 0.05,pAdjustMethod = "BH",independentFiltering = TRUE)
res.primary <- as.data.table(as.data.frame(res.primary),keep.rownames = "Gene")
res.primary[,`:=`(Predictor = "SCN9A pseudobulk expression",Effect ="log2FC per 1-SD increase in SCN9A expression",Significant_FDR_0.05 = !is.na(padj) & padj < 0.05)]
res.primary <- res.primary[order(is.na(padj),padj,pvalue)]

## 提取成熟β细胞标志物和应激基因
target.definition <- data.table(Gene = target.genes,
                                Category = c(rep("Mature beta-cell marker",4),rep("Stress marker",3)),
                                Expected_direction = c(rep("Positive",4),rep("Negative",3)))
target.primary <- merge(target.definition,res.primary,by = "Gene",all.x = TRUE,sort = FALSE)
target.primary[,Direction_matches_hypothesis :=ifelse(is.na(log2FoldChange),NA,ifelse(Expected_direction =="Positive",log2FoldChange > 0,log2FoldChange < 0))]
target.primary

## 输出Excel结果
stopifnot(identical(colnames(pseudobulk.counts),colnames(normalized.counts)))
donor.names <- pseudobulk.metadata[colnames(pseudobulk.counts),"donor_accession"] ## 将PB_ID替换成donor_accession

# 建立输出副本，不修改后续分析使用的原矩阵
raw.counts.output <- pseudobulk.counts
normalized.counts.output <- normalized.counts
colnames(raw.counts.output) <- donor.names
colnames(normalized.counts.output) <- donor.names
raw.counts.table <- as.data.table(raw.counts.output,keep.rownames = "Gene")
normalized.counts.table <- as.data.table(normalized.counts.output,keep.rownames = "Gene")

openxlsx::write.xlsx(list(Primary_all_genes = res.primary,Primary_target_genes = target.primary,Raw_pseudobulk_counts = raw.counts.table,DESeq2_normalized_counts = normalized.counts.table),file = file.path(result_dir, "SCN9A_continuous_donor_pseudobulk_association_results.xlsx"),overwrite = TRUE)

## 火山图
volcano.dt <- copy(res.primary)
volcano.dt <- volcano.dt[!is.na(log2FoldChange) & !is.na(padj) & Gene!="SCN9A"]
volcano.dt[,neg_log10_FDR := -log10(pmax(padj,1e-300))]
volcano.dt[,Significance := ifelse(padj < 0.05 & log2FoldChange >= 0.5,"Positive association", ifelse(padj < 0.05 & log2FoldChange <= -0.5,"Negative association","Not significant"))]
top.positive <- head(volcano.dt[log2FoldChange > 0.5][order(padj,pvalue),Gene],5)
top.negative <- head(volcano.dt[log2FoldChange < -0.5][order(padj,pvalue),Gene],5)
genes.to.label <- unique(c(top.positive,top.negative))
p.volcano <- ggplot(volcano.dt,aes(x = log2FoldChange,y = neg_log10_FDR,colour = Significance)) +
  geom_point(size = 1.3,alpha = 0.7) +
  geom_vline(xintercept = c(-0.5,0.5),linetype = 2,colour = "grey40") +
  geom_hline(yintercept = -log10(0.05),linetype = 2,colour = "grey40") +
  ggrepel::geom_text_repel(data = volcano.dt[Gene %in% genes.to.label],aes(label = Gene),size = 3.2,max.overlaps = Inf,box.padding = 0.4) +
  scale_colour_manual(values = c("Positive association" = "#D73027","Negative association" = "#4575B4","Not significant" = "grey75")) +
  labs(title = "Genes associated with SCN9A expression in Beta cells",subtitle = "Donor-level pseudobulk analysis",x = "log2FC per 1-SD increase in SCN9A expression",y ="-log10(FDR)",colour = NULL) +
  theme_classic() +
  theme(plot.title =element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),axis.text =element_text(colour = "black"))
ggsave(file.path(result_dir, "SCN9A_continuous_association_volcano.pdf"),p.volcano,width = 7,height = 6)

## Top 50关联基因热图
heatmap.candidates <- res.primary[Gene != "SCN9A" & !is.na(padj) & !is.na(stat) & padj < 0.05]
top.positive <- heatmap.candidates[log2FoldChange > 0][order(-stat)][1:min(25, .N)]
top.negative <- heatmap.candidates[log2FoldChange < 0][order(stat)][1:min(25, .N)]
top50.genes <- c(top.positive$Gene,top.negative$Gene)

## VST转换
vsd <- varianceStabilizingTransformation(dds.primary,blind = FALSE)
vst.matrix <- assay(vsd)

##去除模型中其他协变量影响，同时保留 SCN9A 连续变量效应。这一步仅用于热图展示，不参与DESeq2显著性计算。
preserve.design <- model.matrix(~ SCN9A_PB_scaled,data = as.data.frame(colData(dds.primary)))
design.formula.primary <- reformulate(c("source","chemistry","sex","age_scaled","diabetes_group","SCN9A_PB_scaled"))

covariate.formula <- reformulate(c("source","chemistry","sex","age_scaled","diabetes_group"))
covariate.matrix <- model.matrix(covariate.formula,data = as.data.frame(colData(dds.primary)))
covariate.matrix <- covariate.matrix[,colnames(covariate.matrix) != "(Intercept)",drop = FALSE]
vst.adjusted <- limma::removeBatchEffect(vst.matrix,covariates =covariate.matrix,design = preserve.design)

## 准备热图矩阵
heatmap.matrix <- vst.adjusted[top50.genes,,drop = FALSE]

## 每个基因做Z-score
heatmap.z <- t(scale(t(heatmap.matrix)))
heatmap.z[!is.finite(heatmap.z)] <- 0

## 按SCN9A表达从低到高排列donor
donor.order <- rownames(pseudobulk.metadata)[order(pseudobulk.metadata$SCN9A_PB_log2)]
heatmap.z <- heatmap.z[,donor.order,drop = FALSE]

## 添加列注释：
annotation.col <- data.frame(SCN9A_expression = pseudobulk.metadata[donor.order,"SCN9A_PB_log2"],
                             SCN9A_positive_percent = pseudobulk.metadata[donor.order,"SCN9A_percent_detected"],
                             # Diabetes_group = pseudobulk.metadata[donor.order,"diabetes_group"],
                             # Source = pseudobulk.metadata[donor.order,"source"],
                             # Age = pseudobulk.metadata[donor.order,"age"],
                             # Sex =pseudobulk.metadata[donor.order,"sex"],
                             HbA1c_group = pseudobulk.metadata[donor.order,"hba1c_group"]
)
rownames(annotation.col) <- donor.order

## 输出热图
pheatmap::pheatmap(heatmap.z,cluster_rows = TRUE,cluster_cols = FALSE,
                   show_colnames = FALSE,show_rownames = TRUE,
                   annotation_col = annotation.col,color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),border_color = NA,
                   fontsize_row = 8,filename = file.path(result_dir, "SCN9A_continuous_association_top50_heatmap.pdf"),
                   width = 11,height = 8)

##### Top 100 SCN9A关联基因功能模块热图
heatmap.candidates <- data.table::copy(res.primary[Gene != "SCN9A" & !is.na(padj) & !is.na(stat) & padj < 0.05])
heatmap.candidates[,abs_stat := abs(stat)]
data.table::setorder(heatmap.candidates,padj,-abs_stat)

### VST转换及协变量校正
vsd <- DESeq2::varianceStabilizingTransformation(dds.primary,blind = FALSE)
vst.matrix <- SummarizedExperiment::assay(vsd)

### 保留SCN9A连续变量效应
preserve.design <- model.matrix(~ SCN9A_PB_scaled,data = as.data.frame(SummarizedExperiment::colData(dds.primary)))

### 去除其他协变量影响
covariate.formula <- reformulate(c("source","chemistry","sex","age_scaled","diabetes_group"))
covariate.matrix <- model.matrix(covariate.formula,data = as.data.frame(SummarizedExperiment::colData(dds.primary)))
covariate.matrix <- covariate.matrix[,colnames(covariate.matrix) != "(Intercept)",drop = FALSE]

vst.adjusted <- limma::removeBatchEffect(vst.matrix,covariates = covariate.matrix,design = preserve.design)

### 选择可用于热图的Top100基因
### 仅保留存在于VST矩阵中的候选基因
candidate.genes <- intersect(heatmap.candidates$Gene,rownames(vst.adjusted))

### 检查校正后表达是否存在变异
candidate.sd <- apply(vst.adjusted[candidate.genes,,drop = FALSE],1,stats::sd,na.rm = TRUE)

eligible.genes <- names(candidate.sd)[is.finite(candidate.sd) & candidate.sd > 0]
heatmap.candidates <- heatmap.candidates[Gene %in% eligible.genes]

data.table::setorder(heatmap.candidates,padj,-abs_stat)

### 选取前100个
top100.table <- head(heatmap.candidates,100)
top100.genes <- top100.table$Gene

### 准备热图矩阵
heatmap.matrix <- vst.adjusted[top100.genes,,drop = FALSE]

### 每个基因内部做Z-score
heatmap.z <- t(scale(t(heatmap.matrix)))
if(any(!is.finite(heatmap.z))){stop("Non-finite values were detected after gene-wise scaling.")}

### donor按照SCN9A表达量从低到高排列
stopifnot(setequal(rownames(pseudobulk.metadata),colnames(heatmap.z)))
donor.order <- rownames(pseudobulk.metadata)[order(pseudobulk.metadata$SCN9A_PB_log2,na.last = TRUE)]

heatmap.z <- heatmap.z[,donor.order,drop = FALSE]
stopifnot(identical(colnames(heatmap.z),donor.order))

### 仅限制绘图颜色范围
heatmap.z.display <- pmax(pmin(heatmap.z,2),-2)

### 使用Pearson相关距离进行基因聚类
gene.correlation <- stats::cor(t(heatmap.z),method = "pearson",use = "pairwise.complete.obs")
gene.correlation[!is.finite(gene.correlation)] <- 0
diag(gene.correlation) <- 1

### Pearson相关距离
gene.distance <- stats::as.dist(1 - gene.correlation)

### 平均连接法层次聚类
gene.hclust <- stats::hclust(gene.distance,method = "average")

### 使用silhouette选择模块数量
candidate.k <- 2:6

silhouette.score <- sapply(candidate.k,function(k){
  current.cluster <- stats::cutree(gene.hclust,k = k)
  mean(cluster::silhouette(current.cluster,gene.distance)[,"sil_width"])})

silhouette.result <- data.table::data.table(Number_of_modules = candidate.k,Mean_silhouette_width = silhouette.score)

best.k <- candidate.k[which.max(silhouette.score)]

cat("Selected number of gene modules:",best.k,"\n")

### 获得基因模块
raw.cluster <- stats::cutree(gene.hclust,k = best.k)

### 计算每个模块的平均表达曲线
module.profile <- do.call(rbind,lapply(sort(unique(raw.cluster)),function(current.cluster){
  colMeans(heatmap.z[raw.cluster == current.cluster,,drop = FALSE])}))

rownames(module.profile) <- sort(unique(raw.cluster))

### 按表达峰位置对模块重新命名
module.peak.position <- apply(module.profile,1,which.max)

ordered.cluster.ids <- names(sort(module.peak.position))

module.rename <- stats::setNames(paste0("M",seq_along(ordered.cluster.ids)),ordered.cluster.ids)

gene.module <- factor(unname(module.rename[as.character(raw.cluster)]),levels = paste0("M",seq_len(best.k)))
names(gene.module) <- rownames(heatmap.z)

### 对每个模块进行GO-BP富集
module.go.list <- list()
module.levels <- levels(gene.module)

module.label.map <- stats::setNames(module.levels,module.levels)

for(current.module in module.levels){
  current.genes <- names(gene.module)[gene.module == current.module]
  current.go <- clusterProfiler::enrichGO(gene = current.genes,universe = rownames(dds.primary),
                                          OrgDb = org.Hs.eg.db,keyType = "SYMBOL",ont = "BP",
                                          pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,
                                          minGSSize = 5,maxGSSize = 500,readable = TRUE)
  current.go.table <- data.table::as.data.table(as.data.frame(current.go))
  if(nrow(current.go.table) > 0){
    current.go.table[,Module := current.module]
    module.go.list[[current.module]] <- current.go.table
  }
}

module.go.result <- data.table::rbindlist(module.go.list,use.names = TRUE,fill = TRUE)

if(nrow(module.go.result) > 0){
  module.go.top5 <- module.go.result[!is.na(p.adjust) & p.adjust < 0.05][order(p.adjust)][,head(.SD,5),by = Module]
  print(module.go.top5[,.(Module,ID,Description,GeneRatio,p.adjust)])
  
  data.table::fwrite(module.go.result,file.path(result_dir, "SCN9A_top100_heatmap_module_GO_BP_results.tsv"),sep = "\t")
}

### 手动确定模块功能名称
module.label.map <- c(
  "M1" = paste0("M1: ER protein quality control / ERAD","\nNegative association with SCN9A"),
  "M2" = paste0("M2: Ca2+-regulated secretion","\nand vesicle exocytosis","\nPositive association with SCN9A")
)

### 需要在热图中标出的代表基因
genes.to.mark <- c(
  ### M1：ERAD、UPR及蛋白质质量控制
  "SGTA","RHBDD2","EDEM2","DNAJB2","ERP44","DERL2","HERPUD1",
  
  ### M2：钙转运及囊泡释放
  "HECW2","CACNA2D2","CACNA1A","CACNA1D","TRPM3","NRXN1","PCLO","RPH3A","NLGN1"
)

genes.to.mark.available <- intersect( genes.to.mark,rownames(heatmap.z.display))

missing.mark.genes <- setdiff(genes.to.mark,rownames(heatmap.z.display))

if(length(missing.mark.genes) > 0){
  message("Marked genes not present in the heatmap: ",paste(missing.mark.genes,collapse = ", "))
}

### 确定两个模块在完整聚类树中的显示顺序
global.cluster <- stats::cutree(gene.hclust,k = best.k)

### 全局聚类树中的基因顺序
dendrogram.gene.order <- gene.hclust$labels[gene.hclust$order]

### 聚类树中模块的先后顺序
slice.cluster.ids <- unique(global.cluster[dendrogram.gene.order])

### 转换为M1和M2
slice.module.ids <- unname(module.rename[as.character(slice.cluster.ids)])

### 设置模块标题
slice.title.labels <- unname(module.label.map[slice.module.ids])
stopifnot(length(slice.title.labels) == best.k,!anyNA(slice.title.labels))

### 准备donor顶部注释
annotation.metadata <- as.data.frame(pseudobulk.metadata)
annotation.metadata <- annotation.metadata[match(donor.order,rownames(annotation.metadata)),,drop = FALSE]
stopifnot(identical(rownames(annotation.metadata),colnames(heatmap.z.display)))

annotation.col <- data.frame(SCN9A_expression = as.numeric(annotation.metadata$SCN9A_PB_log2),
                             SCN9A_positive_percent = as.numeric(annotation.metadata$SCN9A_percent_detected),
                             HbA1c_group = factor(as.character(annotation.metadata$hba1c_group),
                                                  levels = c("Normal","Prediabetes","Diabetes")),
                             row.names = donor.order,check.names = FALSE)

### 设置顶部注释颜色
SCN9A.expression.breaks <- c(
  min(annotation.col$SCN9A_expression,na.rm = TRUE),
  stats::median(annotation.col$SCN9A_expression,na.rm = TRUE),
  max(annotation.col$SCN9A_expression,na.rm = TRUE))

SCN9A.percent.breaks <- c(min(annotation.col$SCN9A_positive_percent,na.rm = TRUE),
                          stats::median(annotation.col$SCN9A_positive_percent,na.rm = TRUE),
                          max(annotation.col$SCN9A_positive_percent,na.rm = TRUE))

SCN9A.expression.colour <- circlize::colorRamp2(SCN9A.expression.breaks,c("#F2F2F2","#E6A57E","#B2182B"))

SCN9A.percent.colour <- circlize::colorRamp2(SCN9A.percent.breaks,c("#F2F2F2","#80B1D3","#2166AC"))

column.annotation <- ComplexHeatmap::HeatmapAnnotation(df = annotation.col,
                                                       col = list(SCN9A_expression =SCN9A.expression.colour,
                                                                  SCN9A_positive_percent =SCN9A.percent.colour,
                                                                  HbA1c_group = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")),
                                                       na_col = "#D9D9D9",
                                                       annotation_name_side = "left",
                                                       annotation_name_gp = grid::gpar(fontsize = 9),
                                                       simple_anno_size = grid::unit(4,"mm"),
                                                       gap = grid::unit(1,"mm"),
                                                       annotation_legend_param = list(SCN9A_expression = list(title = "SCN9A expression"),
                                                                                      SCN9A_positive_percent = list(title = "SCN9A-positive cells"),
                                                                                      HbA1c_group = list(title = "HbA1c group")))

### 设置热图颜色
heatmap.colour <- circlize::colorRamp2(c(-2,0,2),c("#3B0F70","#1B1B1B","#FDE725"))

### 创建右侧空白标注区域
gene.annotation <- ComplexHeatmap::rowAnnotation(Marked_genes = ComplexHeatmap::anno_empty(border = FALSE,
                                                                                           width = grid::unit(38,"mm")),
                                                 show_annotation_name = FALSE)

### 创建热图对象
ht <- ComplexHeatmap::Heatmap(
  heatmap.z.display,
  name = "Gene-wise\nZ-score",
  
  col = heatmap.colour,
  
  top_annotation = column.annotation,
  
  ### donor不聚类
  cluster_columns = FALSE,
  
  show_column_names = FALSE,
  
  ### 使用完整的基因聚类树
  cluster_rows = gene.hclust,
  
  row_dend_reorder = FALSE,
  
  show_row_dend = TRUE,
  
  ### 在完整聚类树上切成M1和M2
  row_split = best.k,
  
  row_gap = grid::unit(3,"mm"),
  
  ### 模块标题
  row_title = slice.title.labels,
  
  row_title_rot = 0,
  
  row_title_gp = grid::gpar(fontsize = 9,fontface = "bold"),
  
  show_row_names = FALSE,
  
  ### 右侧基因标注区域
  right_annotation = gene.annotation,
  
  row_dend_width = grid::unit(22,"mm"),
  
  ### 不绘制单元格边框
  rect_gp = grid::gpar(col = NA,lwd = 0),
  
  border = FALSE,
  
  ### 使用高分辨率栅格消除PDF中的白色缝隙
  use_raster = FALSE,
  
  heatmap_legend_param = list(title = "Gene-wise\nZ-score",
                              at = c(-2,0,2),
                              labels = c("-2","0","2"),
                              legend_height = grid::unit(35,"mm"))
)

### 打开PDF并绘制热图
heatmap.pdf.file <- file.path(result_dir, "SCN9A_continuous_association_top100_module_heatmap.pdf")
grDevices::cairo_pdf(heatmap.pdf.file,width = 15,height = 8.5)

ht.drawn <- ComplexHeatmap::draw(ht,heatmap_legend_side = "right",annotation_legend_side = "right")

### 提取最终热图中的真实行顺序
final.row.order.list <- ComplexHeatmap::row_order(ht.drawn)

if(!is.list(final.row.order.list)){final.row.order.list <- list(final.row.order.list)}

final.position.list <- vector("list",length(final.row.order.list))

marked.position.list <- vector("list",length(final.row.order.list))

### 按最终位置绘制所有代表基因牵引线
for(current.slice in seq_along(final.row.order.list)){
  current.matrix.rows <- as.integer(final.row.order.list[[current.slice]])
  current.genes <- rownames(heatmap.z.display)[current.matrix.rows]
  slice.offset <- if(current.slice == 1){0} else {sum(lengths(final.row.order.list)[seq_len(current.slice - 1)])}
  
  ## 保存当前模块全部基因的最终位置
  final.position.list[[current.slice]] <- data.table::data.table(
    Heatmap_row_from_top = slice.offset + seq_along(current.genes),
    Slice_from_top =  current.slice,
    Position_in_slice = seq_along(current.genes),
    Input_matrix_row = current.matrix.rows,
    Gene = current.genes,
    Module = as.character(gene.module[current.genes])
  )
  
  ## 当前模块中需要标注的基因
  current.mark.positions <- which(current.genes %in% genes.to.mark.available)
  
  if(length(current.mark.positions) == 0){next}
  
  current.labels <- current.genes[current.mark.positions]
  
  number.rows.in.slice <- length(current.genes)
  
  ## 基因真实纵坐标
  true.y <- 1 - (current.mark.positions - 0.5) / number.rows.in.slice
  
  ## 标签初始位置保持在真实基因附近
  label.y <- true.y
  
  ## 最小标签间距
  min.label.gap <- 0.045
  if(length(label.y) > 1){
    ## 从上到下消除重叠
    for(i in 2:length(label.y)){
      maximum.allowed.position <- label.y[i - 1] - min.label.gap
      
      if(label.y[i] > maximum.allowed.position){label.y[i] <- maximum.allowed.position}
    }
    
    ## 超出下边界时整体上移
    if(min(label.y) < 0.04){label.y <- label.y + (0.04 - min(label.y))}
    
    ## 从下到上再次检查
    for(i in (length(label.y) - 1):1){
      minimum.allowed.position <- label.y[i + 1] + min.label.gap
      if(label.y[i] < minimum.allowed.position){label.y[i] <- minimum.allowed.position}
    }
    
    ## 超出上边界时整体下移
    if(max(label.y) > 0.96){label.y <- label.y - (max(label.y) - 0.96)}
  }
  
  ## 手动画短牵引线和基因名称
  ComplexHeatmap::decorate_annotation("Marked_genes",slice = current.slice,{
    ### 第一段：从真实基因位置水平引出
    grid::grid.segments(x0 = grid::unit(0,"npc"),
                        y0 = grid::unit(true.y,"npc"),
                        x1 = grid::unit(0.08,"npc"),
                        y1 = grid::unit(true.y,"npc"),
                        gp = grid::gpar(col = "grey45",lwd = 0.7)
    )
    
    ## 第二段：连接至标签
    grid::grid.segments(x0 = grid::unit(0.08,"npc"),
                        y0 = grid::unit(true.y,"npc"),
                        x1 = grid::unit(0.15,"npc"),
                        y1 = grid::unit(label.y,"npc"),
                        gp = grid::gpar(col = "grey45",lwd = 0.7)
    )
    
    ## 基因名称
    grid::grid.text(label = current.labels,
                    x = grid::unit(0.18,"npc"),y = grid::unit(label.y,"npc"),
                    just = "left",gp = grid::gpar(fontsize = 7.5,fontface = "italic",col = "black"))
  }
  )
  
  ## 保存标注基因的位置
  marked.position.list[[current.slice]] <- data.table::data.table(
    Gene =current.labels,
    Module = as.character(gene.module[current.labels]),
    Slice_from_top = current.slice,
    Position_in_slice = current.mark.positions,
    Heatmap_row_from_top = slice.offset + current.mark.positions,
    Input_matrix_row = current.matrix.rows[current.mark.positions],
    True_y_npc = true.y,
    Label_y_npc = label.y,
    Label_displacement = label.y - true.y)
}

### 关闭PDF
grDevices::dev.off()

### 整理最终基因顺序
final.gene.position <- data.table::rbindlist(final.position.list,use.names = TRUE,fill = TRUE)
data.table::setorder(final.gene.position,Heatmap_row_from_top)
marked.gene.position <- data.table::rbindlist(marked.position.list,use.names = TRUE,fill = TRUE)
data.table::setorder(marked.gene.position,Heatmap_row_from_top)

### 整理Top100基因统计结果
top100.statistics.output <- data.table::copy(top100.table)
top100.statistics.output[,Selection_rank := seq_len(.N)]
top100.statistics.output[,Module := as.character(gene.module[Gene])]
top100.statistics.output[,Module_description := unname(module.label.map[Module])]
top100.statistics.output[,Association_direction := data.table::fifelse(log2FoldChange > 0,"Positive",data.table::fifelse(log2FoldChange < 0,"Negative","No change"))]
top100.statistics.output[,Marked_in_heatmap := Gene %in% genes.to.mark.available]

statistics.columns <- intersect(c("Selection_rank","Gene","Module","Module_description",
                                  "Association_direction","Marked_in_heatmap","baseMean",
                                  "log2FoldChange","lfcSE","stat","pvalue","padj","abs_stat"),
                                colnames(top100.statistics.output))

top100.statistics.output <- top100.statistics.output[,..statistics.columns]

### 输出与热图完全一致的绘图矩阵
final.heatmap.genes <- final.gene.position$Gene

stopifnot(setequal(final.heatmap.genes,rownames(heatmap.z.display)))

heatmap.plot.output <- data.frame(Heatmap_row_from_top = final.gene.position$Heatmap_row_from_top,
                                  Slice_from_top = final.gene.position$Slice_from_top,
                                  Position_in_slice = final.gene.position$Position_in_slice,
                                  Gene = final.heatmap.genes,
                                  Module = as.character(gene.module[final.heatmap.genes]),
                                  heatmap.z.display[final.heatmap.genes,donor.order,drop = FALSE],
                                  check.names = FALSE)

### 验证Excel中的数值与热图输入矩阵完全一致
heatmap.export.matrix <- as.matrix(heatmap.plot.output[,donor.order,drop = FALSE])
storage.mode(heatmap.export.matrix) <- "numeric"
stopifnot(isTRUE(all.equal(unname(heatmap.export.matrix),unname(heatmap.z.display[final.heatmap.genes,donor.order,drop = FALSE]),tolerance = 0,check.attributes = FALSE)))

### 输出donor注释
donor.annotation.output <- data.frame(Donor = rownames(annotation.col),annotation.col,row.names = NULL,check.names = FALSE)

### 输出到同一个Excel
openxlsx::write.xlsx(list(Top100_statistics = as.data.frame(top100.statistics.output),
                          Heatmap_plot_data = heatmap.plot.output,
                          Donor_annotation =  donor.annotation.output,
                          Marked_gene_positions = as.data.frame(marked.gene.position),
                          Silhouette_results = as.data.frame(silhouette.result),
                          Module_GO_BP = as.data.frame(module.go.result)),
                     file = file.path(result_dir, "SCN9A_top100_module_heatmap_data.xlsx"),
                     overwrite = TRUE,rowNames = FALSE)

