# HbA1c-informed Monocle2 beta-cell pseudotime analysis
# =====================================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

########## monocle2
beta.hpap.v3 <- sr[,sr$Cell_Type=="Beta" & sr$source == "HPAP" & sr$chemistry == "V3"]
meta.hpap.v3 <- as.data.table(beta.hpap.v3@meta.data,keep.rownames = "cell")
meta.hpap.v3[,donor_accession := trimws(as.character(donor_accession))]
meta.hpap.v3[,hba1c_value := as.numeric(as.character(get("hba1c_(percentage)")))]

meta.hpap.v3[,sex := trimws(as.character(sex))]
meta.hpap.v3[,age_value := as.numeric(as.character(get("age_(years)")))]
meta.hpap.v3[,hba1c_group := as.character(hba1c_group)]

meta.hpap.v3 <- meta.hpap.v3[!is.na(hba1c_value)]

##### 统计每个donor的β细胞数
donor.cell.number <- meta.hpap.v3[,.(n_cells = .N),by = .(donor_accession,hba1c_group)]
setorder(donor.cell.number,hba1c_group,-n_cells)

##### 排除 β 细胞数少于min.cells.per.donor的 donor
min.cells.per.donor <- 30
keep.donors <- donor.cell.number[n_cells >= min.cells.per.donor,donor_accession]
meta.hpap.v3 <- meta.hpap.v3[donor_accession %in% keep.donors]

##### 按donor聚合pseudobulk counts
raw.counts.hpap.v3 <- beta.hpap.v3@assays$RNA@counts[,meta.hpap.v3$cell,drop = FALSE]
donor.factor <- factor(meta.hpap.v3$donor_accession)
aggregation.matrix <- Matrix::sparse.model.matrix(~ 0 + donor.factor)
colnames(aggregation.matrix) <- levels(donor.factor)
pseudobulk.counts.hpap.v3 <- raw.counts.hpap.v3 %*% aggregation.matrix

##### 建立donor metadata
donor.metadata.hpap.v3 <- meta.hpap.v3[,.(hba1c_value = hba1c_value[1],
                                          hba1c_group = hba1c_group[1],
                                          sex = sex[1],
                                          age_value = age_value[1],
                                          n_cells = .N),
                                       by = donor_accession]

##### 按照 pseudobulk count 矩阵排列 metadata
donor.metadata.hpap.v3 <- donor.metadata.hpap.v3[match(colnames(pseudobulk.counts.hpap.v3),donor_accession)]

##### 设置因子和年龄标准化
donor.metadata.hpap.v3[,sex := droplevels(factor(sex))]
donor.metadata.hpap.v3[,age_scaled := as.numeric(scale(age_value))]
donor.metadata.hpap.v3[,hba1c_group := factor(hba1c_group,levels = c("Normal","Prediabetes","Diabetes"))]

donor.metadata.hpap.v3 <- as.data.frame(donor.metadata.hpap.v3)
rownames(donor.metadata.hpap.v3) <- donor.metadata.hpap.v3$donor_accession

##### 建立DESeq2对象
pseudobulk.counts.matrix <- as.matrix(pseudobulk.counts.hpap.v3)

##### 建立模型
dds.hba1c <- DESeqDataSetFromMatrix(countData = pseudobulk.counts.matrix,colData = donor.metadata.hpap.v3,design = ~ sex + age_scaled + hba1c_value)

##### 检查设计矩阵是否满秩
design.matrix <- model.matrix(~ sex + age_scaled + hba1c_value,data = donor.metadata.hpap.v3)
c(rank = qr(design.matrix)$rank,columns = ncol(design.matrix)) ## rank 应当等于 columns

##### 过滤低表达基因
keep.genes <- rowSums(counts(dds.hba1c) >= 10) >= 3 ## 保留至少在3个 donor 中 count ≥10 的基因：
table(keep.genes)
dds.hba1c <- dds.hba1c[keep.genes,]

##### 运行连续HbA1c模型
dds.hba1c <- DESeq(dds.hba1c)

##### 提取结果
hba1c.result <- results(dds.hba1c,name = "hba1c_value",alpha = 0.05,independentFiltering = TRUE)
hba1c.result <- as.data.frame(hba1c.result)
hba1c.result$gene <- rownames(hba1c.result)
hba1c.result <- hba1c.result[!is.na(hba1c.result$padj),]
hba1c.result <- hba1c.result[order(-abs(hba1c.result$stat),hba1c.result$padj),]
head(hba1c.result, 20)

##### 生成Monocle2候选ordering genes
# hba1c.ordering.table <- head(hba1c.result,500)
hba1c.ordering.table <- hba1c.result[hba1c.result$padj < 0.05 & abs(hba1c.result$log2FoldChange) >= 0.2,]
hba1c.ordering.table <- hba1c.ordering.table[!grepl("^MT-|^RPL|^RPS|^HBA|^HBB|^HBD|^HBG",hba1c.ordering.table$gene),]
hba1c.ordering.table <- hba1c.ordering.table[hba1c.ordering.table$gene!="MALAT1",]
nrow(hba1c.ordering.table)

##### 每个 donor 最多抽取100个细胞进行轨迹分析
max.cells.per.donor <- 150
set.seed(123)
trajectory.cell.data <- meta.hpap.v3[,.SD[sample(.N,min(.N, max.cells.per.donor))],by = donor_accession]
trajectory.cells <- trajectory.cell.data$cell
length(trajectory.cells)

##### 准备Monocle2输入矩阵
expr.matrix <- beta.hpap.v3@assays$RNA@counts[,trajectory.cells,drop = FALSE]  ## 必须使用原始 counts

##### 准备细胞信息
pdata <- beta.hpap.v3@meta.data[trajectory.cells,,drop = FALSE]
pdata$donor_accession <- trimws(as.character(pdata$donor_accession))
pdata$hba1c_value <- as.numeric(as.character(pdata[["hba1c_(percentage)"]]))
pdata$age_value <- as.numeric(as.character(pdata[["age_(years)"]]))
pdata$sex <- factor(pdata$sex)
pdata$hba1c_group <- factor(pdata$hba1c_group,levels = c("Normal","Prediabetes","Diabetes"))

##### 准备基因信息
fdata <- data.frame(gene_short_name = rownames(expr.matrix),row.names = rownames(expr.matrix),stringsAsFactors = FALSE)

##### 建立 AnnotatedDataFrame：
pd <- new("AnnotatedDataFrame",data = pdata)
fd <- new("AnnotatedDataFrame",data = fdata)

##### 创建Monocle2对象
monocle_cds.hba1c <- newCellDataSet(expr.matrix,phenoData = pd,featureData = fd,expressionFamily = VGAM::negbinomial.size(),lowerDetectionLimit = 1)
monocle_cds.hba1c <- estimateSizeFactors(monocle_cds.hba1c)

##### 统计基因表达细胞数
monocle_cds.hba1c <- monocle::detectGenes(monocle_cds.hba1c,min_expr = 0.1)
summary(Biobase::fData(monocle_cds.hba1c)$num_cells_expressed)

##### 排除在线性轨迹抽样数据中表达过少的ordering_genes,这里要求至少在5%的细胞中检测到。
ordering_genes.hba1c <- hba1c.ordering.table$gene
min.expressed.cells <- ceiling(0.01 * ncol(monocle_cds.hba1c))
expressed.genes <- rownames(Biobase::fData(monocle_cds.hba1c))[Biobase::fData(monocle_cds.hba1c)$num_cells_expressed >= min.expressed.cells]
ordering_genes.hba1c <- ordering_genes.hba1c[ordering_genes.hba1c %in% expressed.genes]
hba1c.result.table <- as.data.table(hba1c.result)
hba1c.result.table[,use_for_ordering:= gene %in% ordering_genes.hba1c]

monocle_cds.hba1c <- setOrderingFilter(monocle_cds.hba1c,ordering_genes.hba1c)

##### 运行DDRTree
set.seed(1234)
system.time({monocle_cds.hba1c <- monocle::reduceDimension(monocle_cds.hba1c,
                                                           max_components = 2,
                                                           reduction_method = "DDRTree",
                                                           norm_method = "log",
                                                           pseudo_expr = 1,
                                                           relative_expr = TRUE,
                                                           scaling = TRUE,
                                                           auto_param_selection = FALSE,
                                                           ncenter = 100,
                                                           maxIter = 20,
                                                           verbose = TRUE)
})

##### 运行orderCells
system.time({monocle_cds.hba1c <- monocle::orderCells(monocle_cds.hba1c)})

##### 初步查看轨迹
p <- monocle::plot_cell_trajectory(monocle_cds.hba1c,color_by = "hba1c_value",cell_size = 0.5,show_branch_points = TRUE)
p <- monocle::plot_cell_trajectory(monocle_cds.hba1c,color_by = "State",cell_size = 0.5,show_branch_points = TRUE)

##### 选择起点
monocle_cds.hba1c <- monocle::orderCells(monocle_cds.hba1c,root_state = 6)

##### 画图
p <- monocle::plot_cell_trajectory(monocle_cds.hba1c,color_by = "Pseudotime",cell_size = 0.5,show_branch_points = TRUE)
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_Pseudotime.pdf"),p,width = 7,height = 6,device = grDevices::cairo_pdf)

p <- monocle::plot_cell_trajectory(monocle_cds.hba1c,color_by = "hba1c_value",cell_size = 0.5,show_branch_points = TRUE)
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_HbA1c.pdf"),p,width = 7,height = 6)

p <- monocle::plot_cell_trajectory(monocle_cds.hba1c,color_by = "hba1c_group",cell_size = 0.5,show_branch_points = TRUE)
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_Group.pdf"),p,width = 7,height = 6)

p <- plot_cell_trajectory(monocle_cds.hba1c, color_by = "hba1c_group",cell_size=0.001,show_branch_points=F,show_tree=T,cell_link_size=0.2)+
  facet_wrap(~hba1c_group)+
  scale_colour_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")) +
  theme_classic()+
  theme(aspect.ratio = 1,legend.position = "none",axis.title = element_text(size = 5),axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "cm"),axis.text = element_blank(),
        plot.margin = margin(0,0,0,0,unit = "cm"),axis.line = element_line(arrow = arrow(length = unit(0.1, 'cm')),size = 0.3))
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_Trajectory_by_HbA1cGroup_Facet.pdf"),p,width = 12,height = 4)

gene.normalized <- beta.hpap.v3@assays$RNA@data[c("SCN9A", "INS", "PDX1"),trajectory.cells,drop = FALSE]
pData(monocle_cds.hba1c)$SCN9A_expression <- as.numeric(gene.normalized["SCN9A",rownames(pData(monocle_cds.hba1c))])
p <- plot_cell_trajectory(monocle_cds.hba1c,color_by = "SCN9A_expression",
                          cell_size = 0.35,show_tree = TRUE,show_branch_points = TRUE,
                          cell_link_size = 0.4) +
  scale_color_viridis_c(option = "magma",name = "SCN9A\nexpression") +
  labs(title = "SCN9A expression along the β-cell trajectory",subtitle = "HPAP V3, HbA1c-informed Monocle2 trajectory") +
  theme_classic()
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_SCN9A_Trajectory.pdf"),p,width = 7,height = 6,device = grDevices::cairo_pdf)

SCN9A.expression <- pData(monocle_cds.hba1c)$SCN9A_expression
SCN9A.positive <- SCN9A.expression[SCN9A.expression > 0]
SCN9A.upper <- as.numeric(quantile(SCN9A.positive,probs = 0.99,na.rm = TRUE))
pData(monocle_cds.hba1c)$SCN9A_display <- pmin(SCN9A.expression,SCN9A.upper)
p.SCN9A.enhanced <- plot_cell_trajectory(monocle_cds.hba1c,color_by = "SCN9A_display",cell_size = 0.55,show_tree = FALSE,show_branch_points = FALSE)
p.SCN9A.enhanced$data <- p.SCN9A.enhanced$data[order(p.SCN9A.enhanced$data$SCN9A_display),,drop = FALSE] ## 将低表达细胞先画、高表达细胞后画，避免阳性细胞被遮住
p.SCN9A.enhanced <- p.SCN9A.enhanced +
  scale_color_gradientn(colours = c(
    "#D9D9D9",  # 0表达：浅灰色
    "#3B4CC0",  # 低表达：蓝色
    "#7E03A8",  # 中低表达：紫色
    "#E64B35",  # 中高表达：红色
    "#FDE725"   # 高表达：黄色
  ),
  values = c(0,0.05,0.25,0.65,1),
  limits = c(0,SCN9A.upper),
  oob = scales::squish,
  trans = "sqrt",name = "SCN9A\nexpression") +
  labs(title = "SCN9A expression along the Beta-cell trajectory",subtitle = "HPAP V3, HbA1c-informed Monocle2 trajectory",x = "Component 1",y = "Component 2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),plot.subtitle = element_text(hjust = 0.5,size = 11))
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_SCN9A_Expression_Enhanced.pdf"),p.SCN9A.enhanced,width = 7,height = 6,device = grDevices::cairo_pdf)

##### 
gene.normalized <- beta.hpap.v3@assays$RNA@data[c("SCN9A", "INS", "PDX1"),trajectory.cells,drop = FALSE]
gene.normalized <- as.matrix(gene.normalized)
trend.data <- data.table(cell = trajectory.cells,
                         Pseudotime = as.numeric(pData(monocle_cds.hba1c)[trajectory.cells,"Pseudotime"]),
                         hba1c_group = as.character(pData(monocle_cds.hba1c)[trajectory.cells,"hba1c_group"]),
                         SCN9A = as.numeric(gene.normalized["SCN9A",trajectory.cells]),
                         INS = as.numeric(gene.normalized["INS",trajectory.cells]),
                         PDX1 = as.numeric(gene.normalized["PDX1",trajectory.cells]))
trend.long <- melt(trend.data,id.vars = c("cell","Pseudotime","hba1c_group"),measure.vars = c("SCN9A","INS","PDX1"),variable.name = "Gene",value.name = "Expression")
trend.long[,Expression_scaled := as.numeric(scale(Expression)),by = Gene]
trend.long <- trend.long[is.finite(Pseudotime) & is.finite(Expression_scaled)]
trend.long[,Gene := factor(Gene,levels = c("SCN9A","PDX1","INS"))]
p <- ggplot(trend.long,aes(x = Pseudotime,y = Expression_scaled,colour = Gene)) +
  geom_hline(yintercept = 0,linetype = "dashed",linewidth = 0.35,colour = "grey70") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "cs",k = 5),method.args = list(method = "REML"),se = FALSE,linewidth = 1.4) +
  scale_colour_manual(values = c("SCN9A" = "#E64B35","PDX1" = "#4DBBD5","INS" = "#00A087")) +
  labs(title = "Gene expression along Beta-cell pseudotime",subtitle = "HPAP V3, HbA1c-informed Monocle2 trajectory",x = "Pseudotime",y = "Gene-wise scaled expression",colour = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11,colour = "black"),
        legend.position = "top",
        legend.text = element_text(size = 11,face = "italic"))
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_SCN9A_PDX1_INS_Pseudotime_Scaled.pdf"),p,width = 7,height = 5.5,device = "pdf",useDingbats = FALSE)

### SCN9A / INS / MAFA / UCN3 / ALDH1A3 expression along pseudotime
genes.pseudotime <- c("SCN9A","INS","MAFA","UCN3","ALDH1A3")

## 提取trajectory cells中的Seurat normalized expression
gene.normalized.genes <- beta.hpap.v3@assays$RNA@data[genes.pseudotime,trajectory.cells,drop = FALSE]
gene.normalized.genes <- as.matrix(gene.normalized.genes)

## 建立pseudotime趋势数据
trend.genes.data <- data.table(
  cell = trajectory.cells,
  Pseudotime = as.numeric(pData(monocle_cds.hba1c)[trajectory.cells,"Pseudotime"]),
  hba1c_group = as.character(pData(monocle_cds.hba1c)[trajectory.cells,"hba1c_group"]),
  SCN9A = as.numeric(gene.normalized.genes["SCN9A",trajectory.cells]),
  INS = as.numeric(gene.normalized.genes["INS",trajectory.cells]),
  MAFA = as.numeric(gene.normalized.genes["MAFA",trajectory.cells]),
  UCN3 = as.numeric(gene.normalized.genes["UCN3",trajectory.cells]),
  ALDH1A3 = as.numeric(gene.normalized.genes["ALDH1A3",trajectory.cells]))

trend.genes.long <- melt(trend.genes.data,id.vars = c("cell","Pseudotime","hba1c_group"),
                         measure.vars = genes.pseudotime,variable.name = "Gene",
                         value.name = "Expression") ## 转换成长表

## 每个基因分别进行Z-score,这样不同表达量级的基因可以放在同一个坐标系中比较趋势
trend.genes.long[,
  Expression_scaled := {
    current.mean <- mean(Expression,na.rm = TRUE)
    current.sd <- sd(Expression,na.rm = TRUE)
    if(!is.finite(current.sd) || current.sd == 0){
      rep(NA_real_,.N)
    } else {
      (Expression - current.mean) / current.sd}
  },by = Gene
]

trend.genes.long <- trend.genes.long[is.finite(Pseudotime) & is.finite(Expression_scaled)] ## 去除非有限值

## 固定图例和曲线显示顺序
trend.genes.long[,Gene := factor(Gene,levels = c("SCN9A","INS","MAFA","UCN3","ALDH1A3"))]

## 设置颜色
gene.colors.genes <- c("SCN9A"="#E64B35","INS"="#00A087","MAFA"="#4DBBD5","UCN3"="#3C5488","ALDH1A3"="#EFC000")

## 绘制GAM pseudotime曲线
p <- ggplot(trend.genes.long,aes(x = Pseudotime,y = Expression_scaled,colour = Gene)) +
  geom_hline(yintercept = 0,linetype = "dashed",linewidth = 0.35,colour = "grey70") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "cs",k = 5),
              method.args = list(method = "REML"),se = FALSE,linewidth = 1.35) +
  scale_colour_manual(values = gene.colors.genes,breaks = c("SCN9A","INS","MAFA","UCN3","ALDH1A3")) +
  labs(title = "Gene expression along Beta-cell pseudotime",
       subtitle = "HPAP V3, HbA1c-informed Monocle2 trajectory",
       x = "Pseudotime",y = "Gene-wise scaled expression",colour = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 11),
        axis.title = element_text(size = 12),axis.text = element_text(size = 11,colour = "black"),
        legend.position = "top",legend.text = element_text(size = 11,face = "italic"))
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_SCN9A_INS_MAFA_UCN3_ALDH1A3_Pseudotime_Scaled.pdf"),
       p,width = 6,height = 5,device = "pdf",useDingbats = FALSE)



#####
### 获取donor层面的VST表达矩阵
vst.hba1c <- varianceStabilizingTransformation(dds.hba1c,blind = FALSE)
vst.expression <- assay(vst.hba1c)

### 提取SCN9A表达量
SCN9A.vst <- as.numeric(vst.expression["SCN9A", ])
names(SCN9A.vst) <- colnames(vst.expression)

### 
hba1c.result.table[,SCN9A_spearman_rho := NA_real_]
hba1c.result.table[,SCN9A_correlation_pvalue := NA_real_]
### 逐个基因计算与SCN9A的相关性
for (i in seq_len(nrow(hba1c.result.table))) {
  current.gene <- hba1c.result.table$gene[i]
  current.expression <- as.numeric(vst.expression[current.gene, ])
  complete.samples <- is.finite(current.expression) & is.finite(SCN9A.vst)
  if (sum(complete.samples) >= 3 &&
      sd(current.expression[complete.samples]) > 0 &&
      sd(SCN9A.vst[complete.samples]) > 0
  ) {
    correlation.test <- suppressWarnings(cor.test(current.expression[complete.samples],SCN9A.vst[complete.samples],method = "spearman",exact = FALSE))
    hba1c.result.table$SCN9A_spearman_rho[i] <- as.numeric(correlation.test$estimate)
    hba1c.result.table$SCN9A_correlation_pvalue[i] <- correlation.test$p.value
  }
}
hba1c.result.table[,SCN9A_correlation_padj := p.adjust(SCN9A_correlation_pvalue,method = "BH")]
fwrite(hba1c.result.table,file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Continuous_DESeq2_All_Genes_Annotated.csv"))

#####
stress.genes <- c("DDIT3","DNAJC3","XBP1","ATF6","HSPA5","EIF2A","PPP1R15A","SOD1",
                  "PRDX1","PRDX4","GPX3","GPX7","CAT","GSR","OXR1","TXNRD1","NQO1",
                  "HMOX1","FTH1","SRXN1","VMP1","TIMP1","FAS")
stress.correlation.table <- hba1c.result.table[gene %in% stress.genes][order(-abs(SCN9A_spearman_rho),SCN9A_correlation_pvalue)]
top.stress.genes <- stress.correlation.table$gene[1:5]
gene.normalized.stress <- beta.hpap.v3@assays$RNA@data[c("SCN9A",top.stress.genes),trajectory.cells,drop = FALSE]
gene.normalized.stress <- as.matrix(gene.normalized.stress)
trend.stress.data <- as.data.table(t(gene.normalized.stress),keep.rownames = "cell")
trend.stress.data[,Pseudotime := as.numeric(pData(monocle_cds.hba1c)[cell,"Pseudotime"])]
trend.stress.long <- melt(trend.stress.data,id.vars = c("cell","Pseudotime"),measure.vars = c("SCN9A",top.stress.genes),variable.name = "Gene",value.name = "Expression")
trend.stress.long[,Expression_scaled := as.numeric(scale(Expression)),by = Gene]

gene.colors <- c("#E64B35","#4DBBD5","#00A087","#7E03A8","#EFC000","#3C5488")
names(gene.colors) <- c("SCN9A",top.stress.genes)
p <- ggplot(trend.stress.long,aes(x = Pseudotime,y = Expression_scaled,colour = Gene)) +
  geom_hline(yintercept = 0,linetype = "dashed",linewidth = 0.35,colour = "grey70") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "cs",k = 5),
              method.args = list(method = "REML"),se = FALSE,linewidth = 1.2) +
  scale_colour_manual(values = gene.colors) +
  labs(title = "Stress-related gene expression along Beta-cell pseudotime",
       subtitle = "Top stress genes associated with SCN9A",
       x = "Pseudotime",y = "Gene-wise scaled expression",colour = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11,colour = "black"),
        legend.position = "right",legend.text = element_text(size = 10,face = "italic"))
ggsave(file.path(result_dir, "monocle2/HPAP_V3_Beta_HbA1c_Monocle2_SCN9A_Top5_Stress_Genes_Pseudotime_Scaled.pdf"),p,width = 8,height = 6,device = "pdf",useDingbats = FALSE)

##### HPAP Top10 regulon activity及SCN9A表达随拟时序变化
### 对齐Monocle2、Seurat和SCENIC中的细胞
monocle.cells <- rownames(Biobase::pData(monocle_cds.hba1c))
all(monocle.cells %in% colnames(beta.hpap.v3))
all(monocle.cells %in% colnames(auc.hpap))
pseudotime.data <- data.table(cell = monocle.cells,Pseudotime = as.numeric(Biobase::pData(monocle_cds.hba1c)[monocle.cells,"Pseudotime"]))

### 提取Top10 regulon的细胞水平AUC
top10.auc.matrix <- auc.hpap[top10.regulons.hpap,monocle.cells,drop = FALSE]
regulon.trend.long <- data.table(cell = rep(monocle.cells,each = length(top10.regulons.hpap)),
                                 Feature = rep(top10.regulons.hpap,times = length(monocle.cells)),
                                 Value = as.numeric(top10.auc.matrix),
                                 Measurement = "Regulon activity") ## 转成长表

### 提取SCN9A的Seurat标准化表达量
SCN9A.expression <- as.numeric(beta.hpap.v3@assays$RNA@data["SCN9A",monocle.cells,drop = TRUE])
SCN9A.trend.long <- data.table(cell = monocle.cells,Feature = "SCN9A expression",Value = SCN9A.expression,Measurement = "Gene expression")

### 合并regulon AUC和SCN9A表达量
trend.long <- rbindlist(list(regulon.trend.long,SCN9A.trend.long),use.names = TRUE)
trend.long <- merge(trend.long,pseudotime.data,by = "cell",all.x = TRUE,sort = FALSE)

### 每个regulon及SCN9A分别进行Z-score
trend.long[,Value_scaled := {current.mean <- mean(Value,na.rm = TRUE)
current.sd <- sd(Value,na.rm = TRUE)
if(!is.finite(current.sd) || current.sd == 0){
  rep(NA_real_,.N)} else {(Value - current.mean) / current.sd}
},
by = Feature]
trend.long <- trend.long[is.finite(Pseudotime) & is.finite(Value_scaled)]

### 固定显示顺序
feature.order <- c(top10.regulons.hpap,"SCN9A expression")
trend.long[,Feature := factor(Feature,levels = feature.order)]

### 设置颜色
regulon.colors <- colorRampPalette(brewer.pal(8,"Dark2"))(length(top10.regulons.hpap))
feature.colors <- c(setNames(regulon.colors,top10.regulons.hpap),"SCN9A expression" = "#E64B35")

### 绘制所有曲线
p <- ggplot(trend.long,aes(x = Pseudotime,y = Value_scaled,colour = Feature)) +
  geom_hline(yintercept = 0,linetype = "dashed",linewidth = 0.35,colour = "grey75") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "cs",k = 5),
              method.args = list(method = "REML"),se = FALSE,linewidth = 1.15) +
  scale_colour_manual(values = feature.colors,breaks = feature.order) +
  labs(title = "Top 10 SCN9A-associated regulon activities along Beta-cell pseudotime",
       subtitle = "HPAP V3, HbA1c-informed Monocle2 trajectory",
       x = "Pseudotime",y = "Feature-wise scaled value",colour = NULL) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 11,colour = "black"),
        legend.position = "right",legend.text = element_text(size = 8))
ggsave(file.path(result_dir, "monocle2/SCN9A_HPAP_Top10_Regulon_Activities_Pseudotime.pdf"),p,width = 10,height = 6.5,device = grDevices::cairo_pdf)

### 每个regulon与SCN9A分别绘图
## 10个regulon分别与SCN9A比较
plot.list <- list()
for(i in seq_along(top10.regulons.hpap)){
  regulon.name <- top10.regulons.hpap[i]
  current.data <- copy(trend.long[as.character(Feature) %in% c(regulon.name,"SCN9A expression")])
  current.data[,Curve := ifelse(as.character(Feature) == "SCN9A expression","SCN9A expression","Regulon activity")]
  p <- ggplot(current.data,aes(x = Pseudotime,y = Value_scaled,colour = Curve)) +
    geom_hline(yintercept = 0,linetype = "dashed",linewidth = 0.3,colour = "grey75") +
    geom_smooth(method = "gam",formula = y ~ s(x,bs = "cs",k = 5),
                method.args = list(method = "REML"),se = FALSE,linewidth = 1.25) +
    scale_colour_manual(values = c("Regulon activity" = "#3B4CC0","SCN9A expression" = "#E64B35")) +
    labs(title = regulon.name,x = "Pseudotime",y = "Feature-wise scaled value") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
          axis.title = element_text(size = 10),axis.text = element_text(size = 9,colour = "black"),legend.position = "none")
  plot.list[[i]] <- p
}

## 手动建立统一图例
shared.legend <- ggplot() + 
  annotate("segment",x = 0.20,xend = 0.28,y = 0.5,yend = 0.5,colour = "#3B4CC0",linewidth = 1.5) +
  annotate("text",x = 0.30,y = 0.5,label = "Regulon activity",hjust = 0,size = 4.2) +
  annotate("segment",x = 0.58,xend = 0.66,y = 0.5,yend = 0.5,colour = "#E64B35",linewidth = 1.5) +
  annotate("text",x = 0.68,y = 0.5,label = "SCN9A expression",hjust = 0,size = 4.2) +
  xlim(0,1) + ylim(0,1) + theme_void()

## 组合为2列×5行
panel.body <- cowplot::plot_grid(plotlist = plot.list,labels = LETTERS[seq_along(plot.list)],ncol = 2,align = "hv",label_size = 12)
p.panel <- cowplot::plot_grid(shared.legend,panel.body,ncol = 1,rel_heights = c(0.035,1))

## 保存
ggsave(file.path(result_dir, "monocle2/SCN9A_HPAP_Top10_Regulon_Pseudotime_Panel.pdf"),p.panel,width = 8,height = 13,device = grDevices::cairo_pdf,limitsize = FALSE)



