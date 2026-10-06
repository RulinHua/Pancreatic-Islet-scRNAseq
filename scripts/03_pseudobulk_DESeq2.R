# Donor-aware pseudobulk differential expression by cell type
# ===========================================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

########## 不同细胞类型中prediabetes vs Diabetes vs Normal以及T1D vs T2D vs Control Without Diabetes 的差异基因表格
minimum.cells <- 1
minimum.donors <- 3

meta.dt <- as.data.table(sr@meta.data,keep.rownames = "Cell")

## 检查每个metadata列在同一个donor内有多少个不同值
metadata.columns <- setdiff(colnames(meta.dt),c("Cell","donor_accession"))
metadata.check <- data.table(Metadata_column = character(),Maximum_unique_values_per_donor = integer())
for (v in metadata.columns) {
  column.check <- meta.dt[,.(n_unique = uniqueN(get(v),na.rm = F)),by = donor_accession]
  maximum.unique <- max(column.check$n_unique,na.rm = TRUE)
  metadata.check <- rbind(metadata.check,data.table(Metadata_column = v,Maximum_unique_values_per_donor = maximum.unique))
}
## 标记donor-level和cell/sample-level字段
metadata.check[,Metadata_type := ifelse(Maximum_unique_values_per_donor <= 1,"Donor-level","Cell/sample-varying")]
donor.metadata.columns <- metadata.check[Metadata_type == "Donor-level",Metadata_column]
donor.info <- unique(meta.dt[,c("donor_accession",donor.metadata.columns),with = FALSE],by = "donor_accession")

## 统计donor × Cell_Type细胞数量
cell.number <- meta.dt[,.(n_cells = .N),by = .(donor_accession,Cell_Type)]

## 每个donor在某个细胞类型中至少minimum.cells个细胞
cell.number <- cell.number[n_cells >= minimum.cells]
setorder(cell.number,Cell_Type,donor_accession)

## 为每个pseudobulk样本创建唯一ID
cell.number[,PB_ID := sprintf("PB%04d",seq_len(.N))]

## 将PB_ID对应回每个细胞
cell.metadata <- merge(meta.dt[,.(Cell,donor_accession,Cell_Type)],cell.number[,.(donor_accession,Cell_Type,n_cells,PB_ID)],by = c("donor_accession","Cell_Type"),all = FALSE,sort = FALSE)

## 建立pseudobulk metadata
pseudobulk.metadata <- merge(cell.number,donor.info,by = "donor_accession",all.x = TRUE,sort = FALSE)
setorder(pseudobulk.metadata,PB_ID)
pseudobulk.metadata <- as.data.frame(pseudobulk.metadata)
rownames(pseudobulk.metadata) <- pseudobulk.metadata$PB_ID

## 从sr的原始counts构建pseudobulk count矩阵
## 保证细胞顺序与sr一致
cells.use <- colnames(sr)[colnames(sr) %in% cell.metadata$Cell]
pb.vector <- setNames(cell.metadata$PB_ID,cell.metadata$Cell)
pb.vector <- pb.vector[cells.use]
any(is.na(pb.vector))

## 创建细胞 × PB_ID的稀疏分组矩阵
pb.factor <- factor(unname(pb.vector),levels = pseudobulk.metadata$PB_ID)
aggregation.matrix <- Matrix::sparse.model.matrix(~ 0 + pb.factor)
colnames(aggregation.matrix) <- levels(pb.factor)

## 基因 × 细胞 乘以 细胞 × PB_ID
raw.counts <- sr@assays$RNA@counts[,cells.use,drop = FALSE]
pseudobulk.counts <- raw.counts %*% aggregation.matrix
identical(colnames(pseudobulk.counts),pseudobulk.metadata$PB_ID)

## 检查汇总后是否仍然都是整数
maximum.decimal <- max(abs(pseudobulk.counts@x - round(pseudobulk.counts@x)),na.rm = TRUE)
maximum.decimal
message("Pseudobulk matrix: ",nrow(pseudobulk.counts)," genes × ",ncol(pseudobulk.counts)," pseudobulk samples")

## 创建Excel工作簿
wb <- createWorkbook()
celltypes.use <- sort(unique(as.character(pseudobulk.metadata$Cell_Type)))
## 每个细胞类型运行DESeq2
for (ct in celltypes.use) {
  message("Running cell type: ", ct)
  celltype.result.list <- list()
  ## 每个细胞类型运行两套分组
  for (grouping.name in c("HbA1c_group","Diabetes_status")) {
    ## HbA1c分组设置
    if (grouping.name == "HbA1c_group") {
      group.column <- "hba1c_group"
      group.levels <- c("Normal","Prediabetes","Diabetes")
      comparisons.use <- list(c("Prediabetes", "Normal"),c("Diabetes", "Normal"), c("Diabetes", "Prediabetes"))
    } else {
      ## Diabetes status分组设置
      group.column <- "diabetes_group"
      group.levels <- c("Control","T1D","T2D")
      comparisons.use <- list(c("T1D", "Control"),c("T2D", "Control"),c("T2D", "T1D"))
    }
    metadata.ct <- pseudobulk.metadata[pseudobulk.metadata$Cell_Type == ct & !is.na(pseudobulk.metadata[[group.column]]),,drop = FALSE]
    metadata.ct$Group <- factor(as.character(metadata.ct[[group.column]]),levels = group.levels)
    metadata.ct$source <- droplevels(factor(metadata.ct$source))
    
    ## 检查每组donor数
    donor.number <- table(factor(metadata.ct$Group,levels = group.levels))
    
    ## 保留donor数至少达到minimum.donors的组
    groups.keep <- names(donor.number[donor.number >= minimum.donors])
    if (length(groups.keep) < 2) {
      message("  Skipped ", grouping.name,": ",paste(names(donor.number),donor.number,sep = "=",collapse = ", "))
      next
    }
    
    ## 删除donor数不足的组
    metadata.ct <- metadata.ct[as.character(metadata.ct$Group) %in% groups.keep,,drop = FALSE]
    
    ## 重新设置因子，去掉空水平
    metadata.ct$Group <- factor(as.character(metadata.ct$Group),levels = group.levels[group.levels %in% groups.keep])
    metadata.ct$source <- droplevels( factor(metadata.ct$source))
    
    ## 只保留合格组之间的比较
    comparisons.available <- list()
    for (comparison in comparisons.use) {
      if (all(comparison %in%  groups.keep)) {comparisons.available[[length(comparisons.available) + 1]] <- comparison}
    }
    comparison.names <- character(0)
    for (comparison in comparisons.available) {
      comparison.names <- c(comparison.names,paste0(comparison[1],"_vs_",comparison[2]))
    }
    ## 8. 提取pseudobulk counts
    sample.ids <- as.character(metadata.ct$PB_ID)
    counts.ct <- pseudobulk.counts[,sample.ids,drop = FALSE]
    rownames(metadata.ct) <- metadata.ct$PB_ID
    
    counts.ct <- as.matrix(counts.ct)
    storage.mode(counts.ct) <- "integer"
    
    ## 构建设计矩阵
    source.removed <- FALSE
    source.status <- NA_character_
    
    ## source只有一个水平，不加入模型
    if (nlevels(metadata.ct$source) <= 1) {
      design.formula <- ~ Group
      source.status <- paste0("Source not included: only one source level")
    } else { 
      design.formula <- ~ source + Group
      design.matrix <- model.matrix(design.formula,data = metadata.ct)
      
      ## 如果不满秩，删除source
      if (qr(design.matrix)$rank < ncol(design.matrix)) {
        message("  ",ct," / ",grouping.name,": ~ source + Group is rank-deficient; source removed")
        design.formula <- ~ Group
        source.removed <- TRUE
        source.status <- paste0("Source removed because ~ source + Group was rank-deficient")
      }
    }
    
    ## 运行DESeq2
    dds <- DESeqDataSetFromMatrix(countData = counts.ct,colData = metadata.ct,design = design.formula)
    dds <- DESeq(dds,quiet = TRUE)
    for (comparison in comparisons.available) {
      numerator <- comparison[1]
      denominator <- comparison[2]
      res <- results(dds,contrast = c("Group",numerator,denominator),alpha = 0.05,pAdjustMethod = "BH",independentFiltering = TRUE)
      res.dt <- as.data.table(as.data.frame(res),keep.rownames = "Gene")
      numerator.number <- length(unique(metadata.ct$donor_accession[metadata.ct$Group ==numerator]))
      denominator.number <- length(unique(metadata.ct$donor_accession[metadata.ct$Group ==denominator]))
      res.dt[,`:=`(Cell_Type = ct,Grouping = grouping.name,
                   Contrast = paste0(numerator,"_vs_",denominator),
                   Numerator = numerator,Denominator = denominator,
                   N_numerator =numerator.number,N_denominator =denominator.number,
                   Significant_FDR_0.05 = !is.na(padj) &  padj < 0.05,
                   DEG_FDR_0.05_abs_log2FC_1 = !is.na(padj) &  padj < 0.05 & abs(log2FoldChange) >= 1)]
      setcolorder(res.dt,c("Gene","Cell_Type","Grouping","Contrast","Numerator","Denominator","N_numerator","N_denominator","log2FoldChange","lfcSE","stat","pvalue","padj","Significant_FDR_0.05","DEG_FDR_0.05_abs_log2FC_1"))
      result.name <- paste0(grouping.name,"__",numerator,"_vs_",denominator)
      res.dt <- res.dt[pvalue < 0.05]
      celltype.result.list[[result.name]] <- res.dt
    }
  }
  
  ## 当前细胞类型结果写入Excel sheet
  addWorksheet(wb,sheetName = ct)
  celltype.result <- rbindlist(celltype.result.list,use.names = TRUE, fill = TRUE)
  writeData(wb,sheet = ct,x = celltype.result,colNames = TRUE,rowNames = FALSE)
}
saveWorkbook(wb,file =file.path(result_dir, "PanKbase_islet_pseudobulk_DESeq2_DEGs_by_celltype.xlsx"),overwrite = TRUE)

