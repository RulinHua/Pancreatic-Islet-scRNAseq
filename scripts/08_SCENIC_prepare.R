# SCENIC input preparation
# ========================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

########## SCENIC
cells.hpap <- which(sr$Cell_Type == "Beta" & sr$source == "HPAP" & sr$chemistry == "V3")
beta.hpap.v3 <- sr[, cells.hpap] ## HPAP V3 β细胞
cells.iidp <- which(sr$Cell_Type == "Beta" & sr$source == "IIDP" & sr$chemistry == "V3") 
beta.iidp.v3 <- sr[, cells.iidp] ## IIDP V3 β细胞

beta.hpap.v3 <- NormalizeData(beta.hpap.v3,normalization.method = "LogNormalize", scale.factor = 10000)
beta.iidp.v3 <- NormalizeData(beta.iidp.v3,normalization.method = "LogNormalize", scale.factor = 10000)

##### 根据HPAP V3进行基因过滤
gene.detected.cells <- Matrix::rowSums(beta.hpap.v3@assays$RNA@counts > 0)
min.detected.cells <- max(50,ceiling(0.005 * ncol(beta.hpap.v3@assays$RNA@counts)))
min.detected.cells
genes.kept <- names(gene.detected.cells[gene.detected.cells >= min.detected.cells])
genes.kept <- intersect(genes.kept,rownames(beta.iidp.v3@assays$RNA@counts))
length(genes.kept)
if (!"SCN9A" %in% genes.kept) {stop("SCN9A未通过基因过滤，请先查看SCN9A的检出细胞数")}

##### 按donor平衡抽样用于构建GRN
### 如果直接把全部 HPAP V3 细胞用于建网，细胞数特别多的 donor 会占更大权重。因此每个 donor 最多抽取300个β细胞。
meta.hpap <- as.data.table(beta.hpap.v3@meta.data,keep.rownames = "cell")
table(meta.hpap$donor_accession)
set.seed(123)
## 抽样
cells.hpap.grn <- meta.hpap[,.(cell = sample(cell,size = min(.N, 300))),by = donor_accession]$cell
length(cells.hpap.grn) ## 这部分细胞只用于 GRNBoost2建网 + motif剪枝，HPAP 和 IIDP 的全部细胞仍然用于 AUCell。

##### 准备表达矩阵
### HPAP建网矩阵
expression.hpap.grn <- beta.hpap.v3@assays$RNA@data[genes.kept,cells.hpap.grn]
### HPAP全量AUCell矩阵
expression.hpap.all <- beta.hpap.v3@assays$RNA@data[genes.kept,]
### IIDP全量AUCell矩阵
expression.iidp.all <- beta.iidp.v3@assays$RNA@data[genes.kept,]

##### 生成loom
loom.hpap.grn <- loomR::create(filename = file.path(result_dir, "SCENIC/HPAP_V3_GRN.loom"),data = expression.hpap.grn,do.transpose = TRUE)
loom.hpap.grn$close_all()
loom.hpap.all <- loomR::create(filename = file.path(result_dir, "SCENIC/HPAP_V3_ALL.loom"),data = expression.hpap.all,do.transpose = TRUE)
loom.hpap.all$close_all()
loom.iidp.all <- loomR::create(filename =  file.path(result_dir, "SCENIC/IIDP_V3_ALL.loom"),data = expression.iidp.all,do.transpose = TRUE)
loom.iidp.all$close_all()

