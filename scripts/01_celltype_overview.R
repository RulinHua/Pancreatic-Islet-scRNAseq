# Cell-type overview and marker visualization
# ===========================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

##### umap图谱
## diabetes_status分组细胞类型
colors <- scCancer::getDefaultColors(n = uniqueN(sr$Cell_Type),type = 2)
names(colors) <- sort(unique(sr$Cell_Type))
p1 <- DimPlot(sr, reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("total")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p2 <- DimPlot(sr[,sr$description_of_diabetes_status=="non-diabetic"], reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("Control")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p3 <- DimPlot(sr[,sr$description_of_diabetes_status=="type 1 diabetes"], reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("T1D")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p4 <- DimPlot(sr[,sr$description_of_diabetes_status=="type 2 diabetes"], reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("T2D")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p <- cowplot::plot_grid(p1[[1]],p2[[1]],p3[[1]],p4[[1]],ncol = 4,nrow = 1,align = "hv",axis = "tblr",byrow=T)
fn <- file.path(result_dir, "celltype_diabetes_status_umapplot.pdf")
cowplot::save_plot(fn,p,ncol = 4,nrow =1,base_height = 5,base_width = 6)

## hba1c分组细胞类型
colors <- scCancer::getDefaultColors(n = uniqueN(sr$Cell_Type),type = 2)
names(colors) <- sort(unique(sr$Cell_Type))
p1 <- DimPlot(sr[,!is.na(sr$hba1c_group)], reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("total")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p2 <- DimPlot(sr[,!is.na(sr$hba1c_group) & sr$hba1c_group=="Normal"], reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("Normal")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p3 <- DimPlot(sr[,!is.na(sr$hba1c_group) & sr$hba1c_group=="Prediabetes"], reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("Prediabetes")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p4 <- DimPlot(sr[,!is.na(sr$hba1c_group) & sr$hba1c_group=="Diabetes"], reduction = "umap", group.by = "Cell_Type",cols=colors,label=F,pt.size = 0.1,raster=FALSE)+ggtitle("Diabetes")+theme(plot.title = element_text(hjust = 0.5),aspect.ratio=1,text = element_text(size = 16))
p <- cowplot::plot_grid(p1[[1]],p2[[1]],p3[[1]],p4[[1]],ncol = 4,nrow = 1,align = "hv",axis = "tblr",byrow=T)
fn <- file.path(result_dir, "celltype_hba1c_group_umapplot.pdf")
cowplot::save_plot(fn,p,ncol = 4,nrow =1,base_height = 5,base_width = 6)

##### 分组细胞类型比例条形图
## diabetes_status分组
dt <- data.table(celltype=sr$Cell_Type,group=sr$description_of_diabetes_status)
dt <- dt[,.(count=.N),by=.(celltype,group)]
dt$pct <- 0
dt$group[dt$group=="non-diabetic"] <- "Control"
dt$group[dt$group=="type 1 diabetes"] <- "T1D"
dt$group[dt$group=="type 2 diabetes"] <- "T2D"
dt$pct[dt$group=="Control"] <- dt$count[dt$group=="Control"]/sum(dt$count[dt$group=="Control"])
dt$pct[dt$group=="T1D"] <- dt$count[dt$group=="T1D"]/sum(dt$count[dt$group=="T1D"])
dt$pct[dt$group=="T2D"] <- dt$count[dt$group=="T2D"]/sum(dt$count[dt$group=="T2D"])

dt$group <- factor(dt$group,levels = c("Control","T1D","T2D"))
colors <- scCancer::getDefaultColors(n = uniqueN(sr$Cell_Type),type = 2)
names(colors) <- sort(unique(sr$Cell_Type))
p <- ggplot(dt,aes(x=group,y=pct, fill=celltype))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = colors)+theme_minimal()+
  theme(legend.text = element_text(size = 16),legend.title = element_text(size = 16),axis.title.x = element_blank(),text = element_text(size = 16),axis.text = element_text(color = "black"),axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave(file.path(result_dir, "celltype_diabetes_status_pencentage_barplot.pdf"),p,height = 3.5,width = 5)

## hba1c分组
dt <- data.table(celltype=sr$Cell_Type,group=sr$hba1c_group)
dt <- dt[!is.na(group)]
dt <- dt[,.(count=.N),by=.(celltype,group)]
dt$pct <- 0
dt$pct[dt$group=="Normal"] <- dt$count[dt$group=="Normal"]/sum(dt$count[dt$group=="Normal"])
dt$pct[dt$group=="Prediabetes"] <- dt$count[dt$group=="Prediabetes"]/sum(dt$count[dt$group=="Prediabetes"])
dt$pct[dt$group=="Diabetes"] <- dt$count[dt$group=="Diabetes"]/sum(dt$count[dt$group=="Diabetes"])

dt$group <- factor(dt$group,levels = c("Normal","Prediabetes","Diabetes"))
colors <- scCancer::getDefaultColors(n = uniqueN(sr$Cell_Type),type = 2)
names(colors) <- sort(unique(sr$Cell_Type))
p <- ggplot(dt,aes(x=group,y=pct, fill=celltype))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = colors)+theme_minimal()+
  theme(legend.text = element_text(size = 16),legend.title = element_text(size = 16),axis.title.x = element_blank(),text = element_text(size = 16),axis.text = element_text(color = "black"),axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave(file.path(result_dir, "celltype_hba1c_group_pencentage_barplot.pdf"),p,height = 3.5,width = 5)

### 细胞类型markers
cells <- lapply(sort(unique(sr$Cell_Type)), function(x){
  cells <- colnames(sr)[sr$Cell_Type==x]
  if(length(cells)>300){cells <- sample(cells,300)}
  return(cells)
}) %>% unlist

Idents(sr) <- sr$Cell_Type
diff.expr.genes <- FindAllMarkers(sr[,cells], assay = "RNA",only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
diff.expr.genes <- as.data.table(diff.expr.genes)
dt <- diff.expr.genes[p_val_adj<0.05,.SD[order(avg_log2FC,decreasing = T)][1],by=gene]

genes <- c("PTPRC", "TYROBP", "C1QC","PECAM1", "VWF", "PLVAP", 
           "PRSS1", "PRSS2", "CPA1","COL1A1", "COL3A1", "FAP",
           "RGS5", "NDUFA4L2", "FABP4","CDK1", "MKI67", "TOP2A",
           "GCG", "TTR", "MAFB","PPY", "GHRL","CHRM3", 
           "MUC5B", "KRT19", "TFF1","SST", "HHEX", "RBP4",
           "INS", "IAPP", "NKX6-1","SPP1", "MMP7", "CFTR")
dt <- DotPlot(sr,features=genes, group.by ="Cell_Type",assay = "RNA")$data %>% as.data.table
dt$features.plot <- factor(dt$features.plot,levels = genes)
ct <- sort(unique(sr$Cell_Type))
dt$id <- factor(dt$id,levels = ct)
p <- ggplot(dt)+
  geom_tile(aes(x=features.plot,y=id),fill="#F5F5F5",colour = "white",size=0.2)+
  geom_point(aes(x=features.plot,y=id,color=avg.exp.scaled,size=pct.exp))+
  coord_equal()+
  scale_size_continuous(range = c(0,5))+
  scale_color_gradientn(colors = c("grey","blue"))+
  theme(axis.text.x = element_text(size = 16,angle = 90,hjust = 1,vjust = 0.4,color="black"),axis.text.y = element_text(size=16,color="black"))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.text = element_text(size = 16),legend.title = element_text(size = 16),
        axis.title = element_blank(),axis.ticks = element_blank(),axis.ticks.length = unit(0, "cm"),
        panel.grid =element_blank(),legend.position = "top",panel.background = element_blank())
ggsave(file.path(result_dir, "celltype_markers_dotplot.pdf"),p,height = 4.5,width = 12)

genes <- c("REG1A","CTRB2","PRSS1","PRSS2","CPA1","COL6A1","PDGFRB","RGS5","GCG",
           "PPY","IAPP","INS","MKI67","CDK1","PTPRC","SST","GHRL","KRT19","PLVAP",
           "ESAM","VWF","PECAM1","MUC5B")
dt <- DotPlot(sr,features=genes, group.by ="Cell_Type",assay = "RNA")$data %>% as.data.table
dt$features.plot <- factor(dt$features.plot,levels = genes)
ct <- sort(unique(sr$Cell_Type))
dt$id <- factor(dt$id,levels = rev(sort(unique(as.character(dt$id)))))
p <- ggplot(dt)+
  geom_tile(aes(x=features.plot,y=id),fill="#F5F5F5",colour = "white",size=0.2)+
  geom_point(aes(x=features.plot,y=id,color=avg.exp.scaled,size=pct.exp))+
  coord_equal()+
  scale_size_continuous(range = c(0,5))+
  scale_color_gradientn(colors = c("lightgrey", "red"))+
  theme(axis.text.x = element_text(size = 16,angle = 90,hjust = 1,vjust = 0.4,color="black"),axis.text.y = element_text(size=16,color="black"))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),legend.text = element_text(size = 16),legend.title = element_text(size = 16),
        axis.title = element_blank(),axis.ticks = element_blank(),axis.ticks.length = unit(0, "cm"),
        panel.grid =element_blank(),legend.position = "top",panel.background = element_blank())
ggsave(file.path(result_dir, "Islet_celltype_DotPlot_PanKbase_v3.4_reference_markers.pdf"),p,height = 4.5,width = 12)

