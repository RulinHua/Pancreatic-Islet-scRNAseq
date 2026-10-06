# SCN-family expression across disease and HbA1c groups
# =====================================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

##### SCN1A-SCN11A组间表达
scn_genes <- c("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A","SCN7A", "SCN8A", "SCN9A", "SCN10A", "SCN11A","TRPC4")
### diabetes_status分组
## bubble_plot
lst <- list()
stat_list <- list()
for (i in scn_genes) {
  dt <- DotPlot(sr,features = i,assay = "RNA",group.by = "Cell_Type",
                split.by = "diabetes_group",cols = "YlOrRd",
                dot.scale = 10,scale = FALSE)$data %>% as.data.table
  dt$id <- as.character(dt$id)
  dt$Group <- sub("^.*_(Control|T1D|T2D)$","\\1",dt$id)
  dt$Cell_Type <- sub("_(Control|T1D|T2D)$","",dt$id)
  dt$Group <- factor(dt$Group,levels = c("Control", "T1D", "T2D"))
  dt$Cell_Type <- factor(dt$Cell_Type,levels = rev(levels(sr$Cell_Type)))
  stat_list[[i]] <- dt[, .(Gene = i,Cell_Type = as.character(Cell_Type),Diabetes_Status = as.character(Group),Average_expression = avg.exp,Percent_expressing = pct.exp)]
  lst[[i]] <- ggplot(dt,aes(x = Group,y = Cell_Type,size = pct.exp,colour = avg.exp)) +
    geom_point() +
    coord_equal()+
    scale_size_continuous(name = "% expressing",range = c(0, 10),limits = c(0, max(dt$pct.exp, na.rm = TRUE)), breaks = scales::breaks_pretty(n = 4)) +
    scale_colour_gradient(name = "Average expression",low = "grey90",high = "#B2182B") +
    ggtitle(i) +xlab(NULL) +ylab(NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
          axis.text.x = element_text(angle = 30,hjust = 1),axis.text=element_text(color = "black"),
          text = element_text(size = 13),
          panel.grid = element_line(colour = "grey90"))
}
stat_table <- data.table::rbindlist(stat_list,use.names = TRUE,fill = TRUE)
stat_table <- stat_table[order(Gene,Cell_Type,Diabetes_Status)]
wb <- createWorkbook()
addWorksheet(wb,"diabetes_status")
writeData(wb, "diabetes_status", stat_table,colNames = T,rowNames =F)
p <- cowplot::plot_grid(plotlist = lst,ncol = 4,nrow = 3,align = "hv",axis = "tblr",byrow=T)
cowplot::save_plot(file.path(result_dir, "SCN_family_Control_T1D_T2D_bubble_plot.pdf"),p,ncol = 4,nrow =3,base_height =4.5,base_width = 4.5,limitsize = FALSE)

## 柱状图加扰动点图
meta.dt <- as.data.table(sr@meta.data,keep.rownames = "Cell")
cell.number <- meta.dt[,.(n_cells = .N),by = .(donor_accession,Cell_Type,diabetes_group)]
cell.number <- cell.number[n_cells >= 10] ## 一个donor在一种细胞类型中至少10个细胞
setorder(cell.number,Cell_Type,diabetes_group,donor_accession)
cell.number[,Average_Group := paste0("G",seq_len(.N))]
meta.use <- merge(meta.dt[,.(Cell,donor_accession,Cell_Type,diabetes_group)],cell.number,by = c("donor_accession","Cell_Type","diabetes_group"),all = FALSE,sort = FALSE)

sr.use <- subset(sr,cells = meta.use$Cell)
average.group.vector <- setNames(meta.use$Average_Group,meta.use$Cell)
sr.use$Average_Group <- unname(average.group.vector[colnames(sr.use)])

# 计算donor × Cell Type平均表达量
avg.matrix <- AverageExpression(object = sr.use,assays = "RNA",features = scn_genes,group.by = "Average_Group",layer = "data",return.seurat = FALSE,verbose = FALSE)$RNA
avg.dt <- as.data.table(t(as.matrix(avg.matrix)),keep.rownames = "Average_Group")
donor.mean <- melt(avg.dt,id.vars = "Average_Group",measure.vars = scn_genes, variable.name = "Gene",value.name = "Donor_mean_expression")
donor.mean <- merge(donor.mean,cell.number,by = "Average_Group",all.x = TRUE,sort = FALSE)

# 去掉 donor 数不足的疾病组
donor.count.table <- donor.mean[,.(n_donor = uniqueN(donor_accession)),by = .(Gene,Cell_Type,diabetes_group)]
donor.count.valid <- donor.count.table[n_donor >= 3] ## 每组至少3个donor
valid.combination <- donor.count.valid[,.(n_group = uniqueN(diabetes_group)),by = .(Gene,Cell_Type)][n_group >= 2] ## 每个细胞类型至少两个疾病组
donor.mean <- merge(donor.mean,donor.count.valid[,.(Gene,Cell_Type,diabetes_group)],by = c("Gene","Cell_Type","diabetes_group"),all = FALSE,sort = FALSE)
donor.mean <- merge(donor.mean,valid.combination[,.(Gene,Cell_Type)],by = c("Gene","Cell_Type"),all = FALSE,sort = FALSE)

# 计算柱高和 SEM
plot.summary <- donor.mean[,.(Group_mean_expression = mean(Donor_mean_expression),SD = sd(Donor_mean_expression),n_donor = uniqueN(donor_accession)),by = .(Gene,Cell_Type,diabetes_group)]
plot.summary[,SEM := SD / sqrt(n_donor)]
plot.summary[,Lower := pmax(Group_mean_expression - SEM,0)]
plot.summary[,Upper := Group_mean_expression + SEM]

# 画图
for (i in scn_genes) {
  celltype.plot.list <- list()
  for (ct in sort(unique(donor.mean$Cell_Type))) {
    plot.data <- donor.mean[Gene == i & Cell_Type == ct]
    summary.data <- plot.summary[Gene == i & Cell_Type == ct]
    groups.use <- c("Control","T1D","T2D")
    groups.use <- groups.use[groups.use %in% unique(as.character(plot.data$diabetes_group))]
    plot.data$diabetes_group <- factor(as.character(plot.data$diabetes_group),levels = groups.use)
    summary.data <- summary.data[as.character(diabetes_group) %in% groups.use]
    summary.data$diabetes_group <- factor(as.character(summary.data$diabetes_group),levels = groups.use)
    comparisons.use <- combn(groups.use,2,simplify = FALSE)    
    
    p.small <- ggplot(data = plot.data,aes(x = diabetes_group,y = Donor_mean_expression,group = diabetes_group)) +
      geom_col(data = summary.data,aes(x = diabetes_group,y = Group_mean_expression,fill = diabetes_group),inherit.aes = FALSE,width = 0.68,alpha = 0.75,colour = "black",linewidth = 0.35) +
      geom_errorbar(data = summary.data,aes(x = diabetes_group,ymin = Lower,ymax = Upper),inherit.aes = FALSE,width = 0.16,linewidth = 0.55,colour = "black") +
      geom_jitter(aes(fill = diabetes_group),shape = 21,colour = "black",stroke = 0.25,width = 0.12,height = 0,size = 1.9,alpha = 0.9) +
      ggpubr::stat_compare_means(comparisons = comparisons.use,method = "wilcox.test",method.args = list(exact = FALSE),label = "p.format",step.increase = 0.15,tip.length = 0.01,bracket.size = 0.25,size = 2.8,hide.ns = FALSE) +
      scale_fill_manual(values = c("Control" = "#4DBBD5","T1D" = "#E64B35","T2D" = "#00A087"),drop = TRUE) +
      labs(title = ct,x = NULL,y = "Donor mean expression") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
            axis.text.x = element_text(angle = 30,hjust = 1,colour = "black"),
            axis.text.y = element_text(colour = "black"),legend.position = "none",
            plot.margin = margin(t = 5,r = 5,b = 5,l = 5))
    celltype.plot.list[[ct]] <- p.small
  }
  ncol.plot <- 4
  nrow.plot <- ceiling(length(celltype.plot.list) / ncol.plot)
  p.grid <- cowplot::plot_grid(plotlist = celltype.plot.list, ncol = ncol.plot,nrow = nrow.plot,align = "hv",axis = "tblr",byrow = TRUE)
  p.title <- cowplot::ggdraw() +cowplot::draw_label(i,fontface = "bold",size = 18,x = 0.5,hjust = 0.5)
  p.final <- cowplot::plot_grid(p.title,p.grid,ncol = 1,rel_heights = c(0.06,1))
  fn <- paste0(file.path(result_dir, "SCN_donor_barplots_all_celltypes_diabetes_status/"),i,"_donor_mean_expression_by_celltype_Control_T1D_T2D.pdf")
  ggsave(filename = fn,plot = p.final,width = 16,height = 4.2 * nrow.plot + 0.6,units = "in",limitsize = FALSE)
}

# donor信息表
meta.all <- as.data.table(sr@meta.data,keep.rownames = "Cell")
meta.all <- meta.all[donor_accession %in% unique(donor.mean$donor_accession)] # 只保留实际用于画图的donor
# 检查每个metadata列在同一个donor内有多少个不同值
metadata.columns <- setdiff(colnames(meta.all),c("Cell","donor_accession"))
metadata.check <- data.table(Metadata_column = character(),Maximum_unique_values_per_donor = integer())

for (v in metadata.columns) {
  column.check <- meta.all[,.(n_unique = uniqueN(get(v),na.rm = TRUE)),by = donor_accession]
  maximum.unique <- max(column.check$n_unique,na.rm = TRUE)
  metadata.check <- rbind(metadata.check,data.table(Metadata_column = v,Maximum_unique_values_per_donor = maximum.unique))
}
# 标记donor-level和cell/sample-level字段
metadata.check[,Metadata_type := ifelse(Maximum_unique_values_per_donor <= 1,"Donor-level","Cell/sample-varying")]
donor.metadata.columns <- metadata.check[Metadata_type == "Donor-level",Metadata_column]
donor.metadata.columns <- setdiff(donor.metadata.columns,colnames(donor.mean))

donor.metadata <- unique(meta.all[,c("donor_accession",donor.metadata.columns),with = FALSE],by = "donor_accession")
donor.expression.metadata <- merge(donor.mean,donor.metadata,by = "donor_accession",all.x = TRUE,sort = FALSE)
donor.expression.metadata$Average_Group <- NULL
donor.expression.metadata$n_cells <- NULL
first.columns <- c("Gene","Cell_Type","diabetes_group","donor_accession","Donor_mean_expression")
other.columns <- setdiff(colnames(donor.expression.metadata),first.columns)
setcolorder(donor.expression.metadata,c(first.columns,other.columns))
#检查是否每个 donor–细胞类型–基因只有一行
duplicate.check <- donor.expression.metadata[,.N,by = .(donor_accession,Cell_Type,Gene)][N > 1]
duplicate.check

openxlsx::write.xlsx(x = list(Donor_expression = donor.expression.metadata),file = file.path(result_dir, "SCN_donor_barplots_all_celltypes_diabetes_status/SCN_diabetes_status_donor_expression_and_metadata.xlsx"),rowNames = FALSE,overwrite = TRUE)

## 
for (i in scn_genes) {
  p <- FeaturePlot(sr,features = i,reduction = "umap",split.by = "diabetes_group",
                   cols = c("grey90", "#B2182B"),pt.size = 0.1,order = TRUE,
                   keep.scale = "feature",raster = F,ncol = 3,combine = F)
  ## 从第三张图中准备色条
  legend.plot <- p[[3]] +theme(legend.position = "right",legend.title = element_text(size = 11),legend.text = element_text(size = 10)) +
    guides(colour = guide_colourbar(title = "Expression",title.position = "top",barheight = grid::unit(2.5, "cm"),barwidth = grid::unit(0.4, "cm")))
  ## 提取右侧色条
  # legend.grob <- ggplotGrob(legend.plot)
  # index <- which(legend.grob$layout$name == "guide-box-right")
  # ld <- legend.grob$grobs[[index]]
  ld <- cowplot::get_legend(legend.plot)
  
  umap.theme <- theme(aspect.ratio = 1,plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
                      axis.title.x = element_text(colour = "black"),axis.text.x = element_text(colour = "black"),axis.ticks.x = element_line(colour = "black"),axis.line.x = element_line(colour = "black"),
                      axis.title.y = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),axis.ticks.y = element_line(colour = "black"),
                      axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank())
  p1 <- p[[1]] +NoLegend() +labs(title = "Control",x = "UMAP_1",y = "UMAP_2") +umap.theme
  p2 <- p[[2]] +NoLegend() +labs(title = "T1D",x = "UMAP_1",y = "UMAP_2") + umap.theme
  p3 <- p[[3]] +NoLegend() +labs(title = "T2D",x = "UMAP_1",y = "UMAP_2") + umap.theme
  
  ## 只栅格化UMAP点
  p1 <- ggrastr::rasterise(p1,layers = "Point",dpi = 600)
  p2 <- ggrastr::rasterise(p2,layers = "Point",dpi = 600)
  p3 <- ggrastr::rasterise(p3,layers = "Point",dpi = 600)
  p.row <- cowplot::plot_grid(p1,p2,p3,ncol = 3,nrow = 1,align = "hv",axis = "tblr")
  p.content <- cowplot::plot_grid(p.row,ld,ncol = 2,rel_widths = c(1, 0.1))
  p.title <- cowplot::ggdraw() +cowplot::draw_label(i,fontface = "bold",size = 17,x = 0.5,hjust = 0.5)
  p.final <- cowplot::plot_grid(p.title,p.content,ncol = 1,rel_heights = c(0.08, 1))
  # fn <- paste0(file.path(result_dir, "SCN_UMAP_diabetes_status/"),i, "_Control_T1D_T2D_UMAP.png")
  fn <- paste0(file.path(result_dir, "SCN_UMAP_diabetes_status/"),i, "_Control_T1D_T2D_UMAP.pdf")
  ggsave(filename = fn,plot = p.final,device = cairo_pdf,width = 13,height = 4.8,units = "in",bg = "white")
  # ggsave(fn,p.final,device = ragg::agg_png,width = 13,height = 4.8,units = "in",dpi = 300,bg = "white")
}

### hba1c分组
## bubble plot
sr1 <- sr[,!is.na(sr$hba1c_group)]
lst <- list()
stat_list <- list()
for (i in scn_genes) {
  dt <- DotPlot(sr1,features = i,assay = "RNA",group.by = "Cell_Type",
                split.by = "hba1c_group",cols = "YlOrRd",
                dot.scale = 10,scale = FALSE)$data %>% as.data.table
  dt$id <- as.character(dt$id)
  dt$Group <- sub("^.*_(Normal|Prediabetes|Diabetes)$","\\1",dt$id)
  dt$Cell_Type <- sub("_(Normal|Prediabetes|Diabetes)$","",dt$id)
  dt$Group <- factor(dt$Group,levels = c("Normal", "Prediabetes", "Diabetes"))
  dt$Cell_Type <- factor(dt$Cell_Type,levels = rev(levels(sr1$Cell_Type)))
  stat_list[[i]] <- dt[, .(Gene = i,Cell_Type = as.character(Cell_Type),HbA1c_Group = as.character(Group),Average_expression = avg.exp,Percent_expressing = pct.exp)]
  lst[[i]] <- ggplot(dt,aes(x = Group,y = Cell_Type,size = pct.exp,colour = avg.exp)) +
    geom_point() +
    coord_equal()+
    scale_size_continuous(name = "% expressing",range = c(0, 10),limits = c(0, max(dt$pct.exp, na.rm = TRUE)), breaks = scales::breaks_pretty(n = 4)) +
    scale_colour_gradient(name = "Average expression",low = "grey90",high = "#B2182B") +
    ggtitle(i) +xlab(NULL) +ylab(NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),
          axis.text.x = element_text(angle = 30,hjust = 1),axis.text=element_text(color = "black"),
          text = element_text(size = 13),
          panel.grid = element_line(colour = "grey90"))
}
stat_table <- data.table::rbindlist(stat_list,use.names = TRUE,fill = TRUE)
stat_table <- stat_table[order(Gene,Cell_Type,HbA1c_Group)]
addWorksheet(wb,"HbA1c_Group")
writeData(wb, "HbA1c_Group", stat_table,colNames = T,rowNames =F)
saveWorkbook(wb, file.path(result_dir, "SCN_genes_TRPC4_CellType_expression_summary_DiabetesStatus_HbA1c.xlsx"), overwrite = TRUE)
p <- cowplot::plot_grid(plotlist = lst,ncol = 4,nrow = 3,align = "hv",axis = "tblr",byrow=T)
cowplot::save_plot(file.path(result_dir, "SCN_family_HbA1c_group_bubble_plot.pdf"),p,ncol = 4,nrow =3,base_height =4.5,base_width = 4.5,limitsize = FALSE)

## 柱状图加扰动点图
meta.dt <- as.data.table(sr1@meta.data,keep.rownames = "Cell")
cell.number <- meta.dt[,.(n_cells = .N),by = .(donor_accession,Cell_Type,hba1c_group)]
cell.number <- cell.number[n_cells >= 10] ## 一个donor在一种细胞类型中至少10个细胞
setorder(cell.number,Cell_Type,hba1c_group,donor_accession)
cell.number[,Average_Group := paste0("G",seq_len(.N))]
meta.use <- merge(meta.dt[,.(Cell,donor_accession,Cell_Type,hba1c_group)],cell.number,by = c("donor_accession","Cell_Type","hba1c_group"),all = FALSE,sort = FALSE)

sr.use <- subset(sr1,cells = meta.use$Cell)
average.group.vector <- setNames(meta.use$Average_Group,meta.use$Cell)
sr.use$Average_Group <- unname(average.group.vector[colnames(sr.use)])

# 计算donor × Cell Type平均表达量
avg.matrix <- AverageExpression(object = sr.use,assays = "RNA",features = scn_genes,group.by = "Average_Group",layer = "data",return.seurat = FALSE,verbose = FALSE)$RNA
avg.dt <- as.data.table(t(as.matrix(avg.matrix)),keep.rownames = "Average_Group")
donor.mean <- melt(avg.dt,id.vars = "Average_Group",measure.vars = scn_genes, variable.name = "Gene",value.name = "Donor_mean_expression")
donor.mean <- merge(donor.mean,cell.number,by = "Average_Group",all.x = TRUE,sort = FALSE)

# 去掉 donor 数不足的疾病组
donor.count.table <- donor.mean[,.(n_donor = uniqueN(donor_accession)),by = .(Gene,Cell_Type,hba1c_group)]
donor.count.valid <- donor.count.table[n_donor >= 3] ## 每组至少3个donor
valid.combination <- donor.count.valid[,.(n_group = uniqueN(hba1c_group)),by = .(Gene,Cell_Type)][n_group >= 2] ## 每个细胞类型至少两个疾病组
donor.mean <- merge(donor.mean,donor.count.valid[,.(Gene,Cell_Type,hba1c_group)],by = c("Gene","Cell_Type","hba1c_group"),all = FALSE,sort = FALSE)
donor.mean <- merge(donor.mean,valid.combination[,.(Gene,Cell_Type)],by = c("Gene","Cell_Type"),all = FALSE,sort = FALSE)

# 计算柱高和 SEM
plot.summary <- donor.mean[,.(Group_mean_expression = mean(Donor_mean_expression),SD = sd(Donor_mean_expression),n_donor = uniqueN(donor_accession)),by = .(Gene,Cell_Type,hba1c_group)]
plot.summary[,SEM := SD / sqrt(n_donor)]
plot.summary[,Lower := pmax(Group_mean_expression - SEM,0)]
plot.summary[,Upper := Group_mean_expression + SEM]

# 画图
for (i in scn_genes) {
  celltype.plot.list <- list()
  for (ct in sort(unique(donor.mean$Cell_Type))) {
    plot.data <- donor.mean[Gene == i & Cell_Type == ct]
    summary.data <- plot.summary[Gene == i & Cell_Type == ct]
    groups.use <- c("Normal","Prediabetes","Diabetes")
    groups.use <- groups.use[groups.use %in% unique(as.character(plot.data$hba1c_group))]
    plot.data$hba1c_group <- factor(as.character(plot.data$hba1c_group),levels = groups.use)
    summary.data <- summary.data[as.character(hba1c_group) %in% groups.use]
    summary.data$hba1c_group <- factor(as.character(summary.data$hba1c_group),levels = groups.use)
    comparisons.use <- combn(groups.use,2,simplify = FALSE)    
    
    p.small <- ggplot(data = plot.data,aes(x = hba1c_group,y = Donor_mean_expression,group = hba1c_group)) +
      geom_col(data = summary.data,aes(x = hba1c_group,y = Group_mean_expression,fill = hba1c_group),inherit.aes = FALSE,width = 0.68,alpha = 0.75,colour = "black",linewidth = 0.35) +
      geom_errorbar(data = summary.data,aes(x = hba1c_group,ymin = Lower,ymax = Upper),inherit.aes = FALSE,width = 0.16,linewidth = 0.55,colour = "black") +
      geom_jitter(aes(fill = hba1c_group),shape = 21,colour = "black",stroke = 0.25,width = 0.12,height = 0,size = 1.9,alpha = 0.9) +
      ggpubr::stat_compare_means(comparisons = comparisons.use,method = "wilcox.test",method.args = list(exact = FALSE),label = "p.format",step.increase = 0.15,tip.length = 0.01,bracket.size = 0.25,size = 2.8,hide.ns = FALSE) +
      scale_fill_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35"),drop = TRUE)+
      labs(title = ct,x = NULL,y = "Donor mean expression") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
            axis.text.x = element_text(angle = 30,hjust = 1,colour = "black"),
            axis.text.y = element_text(colour = "black"),legend.position = "none",
            plot.margin = margin(t = 5,r = 5,b = 5,l = 5))
    celltype.plot.list[[ct]] <- p.small
  }
  ncol.plot <- 4
  nrow.plot <- ceiling(length(celltype.plot.list) / ncol.plot)
  p.grid <- cowplot::plot_grid(plotlist = celltype.plot.list, ncol = ncol.plot,nrow = nrow.plot,align = "hv",axis = "tblr",byrow = TRUE)
  p.title <- cowplot::ggdraw() +cowplot::draw_label(i,fontface = "bold",size = 18,x = 0.5,hjust = 0.5)
  p.final <- cowplot::plot_grid(p.title,p.grid,ncol = 1,rel_heights = c(0.06,1))
  fn <- paste0(file.path(result_dir, "SCN_donor_barplots_all_celltypes_HbA1c_group/"),i,"_donor_mean_expression_by_celltype_Normal_Prediabetes_Diabetes.pdf")
  ggsave(filename = fn,plot = p.final,width = 16,height = 4.2 * nrow.plot + 0.6,units = "in",limitsize = FALSE)
}

# donor信息表
meta.all <- as.data.table(sr1@meta.data,keep.rownames = "Cell")
meta.all <- meta.all[donor_accession %in% unique(donor.mean$donor_accession)] # 只保留实际用于画图的donor
# 检查每个metadata列在同一个donor内有多少个不同值
metadata.columns <- setdiff(colnames(meta.all),c("Cell","donor_accession"))
metadata.check <- data.table(Metadata_column = character(),Maximum_unique_values_per_donor = integer())

for (v in metadata.columns) {
  column.check <- meta.all[,.(n_unique = uniqueN(get(v),na.rm = TRUE)),by = donor_accession]
  maximum.unique <- max(column.check$n_unique,na.rm = TRUE)
  metadata.check <- rbind(metadata.check,data.table(Metadata_column = v,Maximum_unique_values_per_donor = maximum.unique))
}
# 标记donor-level和cell/sample-level字段
metadata.check[,Metadata_type := ifelse(Maximum_unique_values_per_donor <= 1,"Donor-level","Cell/sample-varying")]
donor.metadata.columns <- metadata.check[Metadata_type == "Donor-level",Metadata_column]
donor.metadata.columns <- setdiff(donor.metadata.columns,colnames(donor.mean))

donor.metadata <- unique(meta.all[,c("donor_accession",donor.metadata.columns),with = FALSE],by = "donor_accession")
donor.expression.metadata <- merge(donor.mean,donor.metadata,by = "donor_accession",all.x = TRUE,sort = FALSE)
donor.expression.metadata$Average_Group <- NULL
donor.expression.metadata$n_cells <- NULL
first.columns <- c("Gene","Cell_Type","hba1c_group","donor_accession","Donor_mean_expression")
other.columns <- setdiff(colnames(donor.expression.metadata),first.columns)
setcolorder(donor.expression.metadata,c(first.columns,other.columns))
#检查是否每个 donor–细胞类型–基因只有一行
duplicate.check <- donor.expression.metadata[,.N,by = .(donor_accession,Cell_Type,Gene)][N > 1]
duplicate.check

openxlsx::write.xlsx(x = list(Donor_expression = donor.expression.metadata),file = file.path(result_dir, "SCN_donor_barplots_all_celltypes_HbA1c_group/SCN_HbA1c_donor_expression_and_metadata.xlsx"),rowNames = FALSE,overwrite = TRUE)

##
for (i in scn_genes) {
  p <- FeaturePlot(sr1,features = i,reduction = "umap",split.by = "hba1c_group",
                   cols = c("grey90", "#B2182B"),pt.size = 0.1,order = TRUE,
                   keep.scale = "feature",raster = F,ncol = 3,combine = F)
  ## 从第三张图中准备色条
  legend.plot <- p[[3]] +theme(legend.position = "right",legend.title = element_text(size = 11),legend.text = element_text(size = 10)) +
    guides(colour = guide_colourbar(title = "Expression",title.position = "top",barheight = grid::unit(2.5, "cm"),barwidth = grid::unit(0.4, "cm")))
  ## 提取右侧色条
  # legend.grob <- ggplotGrob(legend.plot)
  # index <- which(legend.grob$layout$name == "guide-box-right")
  # ld <- legend.grob$grobs[[index]]
  ld <- cowplot::get_legend(legend.plot)
  
  umap.theme <- theme(aspect.ratio = 1,plot.title = element_text(hjust = 0.5,size = 14,face = "bold"),
                      axis.title.x = element_text(colour = "black"),axis.text.x = element_text(colour = "black"),axis.ticks.x = element_line(colour = "black"),axis.line.x = element_line(colour = "black"),
                      axis.title.y = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),axis.ticks.y = element_line(colour = "black"),
                      axis.title.y.right = element_blank(),axis.text.y.right = element_blank(),axis.ticks.y.right = element_blank(),axis.line.y.right = element_blank())
  p1 <- p[[1]] +NoLegend() +labs(title = "Normal",x = "UMAP_1",y = "UMAP_2") +umap.theme
  p2 <- p[[2]] +NoLegend() +labs(title = "Prediabetes",x = "UMAP_1",y = "UMAP_2") + umap.theme
  p3 <- p[[3]] +NoLegend() +labs(title = "Diabetes",x = "UMAP_1",y = "UMAP_2") + umap.theme
  
  ## 只栅格化UMAP点
  p1 <- ggrastr::rasterise(p1,layers = "Point",dpi = 600)
  p2 <- ggrastr::rasterise(p2,layers = "Point",dpi = 600)
  p3 <- ggrastr::rasterise(p3,layers = "Point",dpi = 600)
  
  p.row <- cowplot::plot_grid(p1,p2,p3,ncol = 3,nrow = 1,align = "hv",axis = "tblr")
  p.content <- cowplot::plot_grid(p.row,ld,ncol = 2,rel_widths = c(1, 0.1))
  p.title <- cowplot::ggdraw() +cowplot::draw_label(i,fontface = "bold",size = 17,x = 0.5,hjust = 0.5)
  p.final <- cowplot::plot_grid(p.title,p.content,ncol = 1,rel_heights = c(0.08, 1))
  # fn <- paste0(file.path(result_dir, "SCN_UMAP_HbA1c_group/"),i, "_Normal_Prediabetes_Diabetes_UMAP.png")
  # ggsave(fn,p.final,device = ragg::agg_png,width = 13,height = 4.8,units = "in",dpi = 300,bg = "white")
  fn <- paste0(file.path(result_dir, "SCN_UMAP_HbA1c_group/"),i, "_Normal_Prediabetes_Diabetes_UMAP.pdf")
  ggsave(filename = fn,plot = p.final,device = cairo_pdf,width = 13,height = 4.8,units = "in",bg = "white")
}

