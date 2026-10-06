# GSEA and predefined functional-module analyses
# ==============================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

########## GSEA
## 使用DESeq2 Wald stat排序
rank.dt <- res.primary[Gene != "SCN9A" & !is.na(stat),.(SYMBOL = Gene,stat,log2FoldChange,pvalue,padj)]

## SYMBOL转换成ENTREZID
gene.map <- clusterProfiler::bitr(rank.dt$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
rank.entrez <- merge(rank.dt,gene.map,by = "SYMBOL",all = FALSE)

## 同一个ENTREZID对应多个SYMBOL时，保留|stat|最大的记录
rank.entrez[, abs_stat := abs(stat)]
setorder(rank.entrez,ENTREZID,-abs_stat)
rank.entrez <- rank.entrez[,.SD[1],by = ENTREZID]

## 创建clusterProfiler需要的命名数值向量
gene.list <- rank.entrez$stat
names(gene.list) <- rank.entrez$ENTREZID

## 从大到小排序
gene.list <- sort(gene.list,decreasing = TRUE)

## 检查
length(gene.list)
anyDuplicated(names(gene.list)) # 应0
all(diff(gene.list) <= 0) # 应TRUE

## 运行GO Biological Process GSEA
set.seed(123)
go.gsea <- clusterProfiler::gseGO(geneList = gene.list,OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",minGSSize = 10,maxGSSize = 500,exponent = 1,eps = 0,pvalueCutoff = 1,pAdjustMethod = "BH",verbose = FALSE,seed = TRUE,by = "fgsea")
go.gsea <- clusterProfiler::setReadable(go.gsea,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
go.table <- as.data.table(as.data.frame(go.gsea))
go.table[,Direction := ifelse(NES > 0,"Positive association with SCN9A","Negative association with SCN9A")]
go.table <- go.table[order(p.adjust, -abs(NES))]

## 运行KEGG GSEA
set.seed(123)
kegg.gsea <- clusterProfiler::gseKEGG(geneList = gene.list,organism = "hsa",keyType = "ncbi-geneid",minGSSize = 10,maxGSSize = 500,exponent = 1,eps = 0,pvalueCutoff = 1,pAdjustMethod = "BH",verbose = FALSE,seed = TRUE,by = "fgsea",use_internal_data = FALSE)
kegg.gsea <- clusterProfiler::setReadable(kegg.gsea,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
kegg.table <- as.data.table(as.data.frame(kegg.gsea))
kegg.table[,Direction := ifelse(NES > 0,"Positive association with SCN9A","Negative association with SCN9A")]
kegg.table <- kegg.table[order(p.adjust, -abs(NES))]

## GO Top通路
go.significant <- go.table[!is.na(p.adjust) & p.adjust < 0.05]
go.positive <- go.significant[NES > 0][order(p.adjust, -NES)]
go.negative <- go.significant[NES < 0][order(p.adjust, NES)]
go.positive <- head(go.positive,5)
go.negative <- head(go.negative,5)
go.top <- rbindlist(list(go.positive,go.negative),use.names = TRUE,fill = TRUE)

## KEGG Top通路
kegg.significant <- kegg.table[!is.na(p.adjust) & p.adjust < 0.05]
kegg.positive <- kegg.significant[NES > 0][order(p.adjust, -NES)]
kegg.negative <- kegg.significant[NES < 0][order(p.adjust, NES)]
kegg.positive <- head(kegg.positive,5)
kegg.negative <- head(kegg.negative,5)
kegg.top <- rbindlist(list(kegg.positive,kegg.negative),use.names = TRUE,fill = TRUE)

## GO Ridge Plot
go.top.ids <- unique(as.character(go.top$ID)) ## 提取前面选出的GO通路ID
go.gsea.top <- go.gsea
go.gsea.top@result <- go.gsea@result[go.gsea@result$ID %in% go.top.ids,,drop = FALSE] ## 只保留正向Top 5和负向Top 5
go.gsea.top@result <- go.gsea.top@result[match(go.top.ids,go.gsea.top@result$ID),,drop = FALSE] ## 按照go.top中的顺序排列
p.go.ridge <- enrichplot::ridgeplot(go.gsea.top,showCategory = nrow(go.gsea.top@result),fill = "p.adjust") +
  labs(title = "Top GO biological processes associated with SCN9A",
       subtitle = "Top positive and negative associations",
       x = "DESeq2 Wald statistic",y = NULL,fill = "FDR") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9,colour = "black"),
        axis.text.x = element_text(colour = "black"))
ggsave(file.path(result_dir, "SCN9A_continuous_GSEA/SCN9A_GO_BP_top5_positive_top5_negative_GSEA_ridgeplot.pdf"),p.go.ridge,width = 10,height = 7)

## KEGG ridge plot
kegg.top.ids <- unique(as.character(kegg.top$ID))
kegg.gsea.top <- kegg.gsea
kegg.gsea.top@result <- kegg.gsea@result[kegg.gsea@result$ID %in% kegg.top.ids,,drop = FALSE]
kegg.gsea.top@result <- kegg.gsea.top@result[match(kegg.top.ids,kegg.gsea.top@result$ID),,drop = FALSE] ## 按照kegg.top中的顺序排列
p.kegg.ridge <- enrichplot::ridgeplot(kegg.gsea.top,showCategory = nrow(kegg.gsea.top@result),fill = "p.adjust") +
  labs(title = "Top KEGG pathways associated with SCN9A",subtitle = "Top positive and negative associations",
       x = "DESeq2 Wald statistic",y = NULL,fill = "FDR") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 9,colour = "black"),
        axis.text.x = element_text(colour = "black"))
ggsave(file.path(result_dir, "SCN9A_continuous_GSEA/SCN9A_KEGG_top5_positive_top5_negative_GSEA_ridgeplot.pdf"),p.kegg.ridge,width = 10,height = 7)

##### 选定的GO KEGG通路山峰图
## 选定的GO-BP通路
go.selected.ids <- c(
  "GO:0007156",  # homophilic cell adhesion...
  "GO:0007416",  # synapse assembly
  "GO:0001702",  # gastrulation with mouth forming second
  "GO:0035249",  # synaptic transmission, glutamatergic
  "GO:0008038",  # neuron recognition
  "GO:0006119",  # oxidative phosphorylation
  "GO:0042773",  # ATP synthesis coupled electron transport
  "GO:0042775",  # mitochondrial ATP synthesis...
  "GO:0019646",  # aerobic electron transport chain
  "GO:0002181"   # cytoplasmic translation
)

## 选定的KEGG通路
kegg.selected.ids <- c(
  "hsa04930",    # Type II diabetes mellitus
  "hsa04724",    # Glutamatergic synapse
  "hsa04070",    # Phosphatidylinositol signaling system
  "hsa04919",    # Thyroid hormone signaling pathway
  "hsa00562",    # Inositol phosphate metabolism
  "hsa03010",    # Ribosome
  "hsa00190",    # Oxidative phosphorylation
  "hsa03050",    # Proteasome
  "hsa05012",    # Parkinson disease
  "hsa05020"     # Prion disease
)

## GO-BP山峰图
go.selected.result <- go.gsea@result[match(go.selected.ids,go.gsea@result$ID),,drop = FALSE]
go.gsea.selected <- go.gsea ## 建立只包含选定通路的gseaResult对象
go.gsea.selected@result <- go.selected.result

p <- enrichplot::ridgeplot(go.gsea.selected,showCategory = nrow(go.gsea.selected@result),
                           fill = "p.adjust",core_enrichment = TRUE,label_format = 45,
                           orderBy = "NES",decreasing = FALSE) +
  labs(title = "Selected GO biological processes associated with SCN9A",
       subtitle = "Positive and negative GSEA enrichment",
       x = "DESeq2 Wald statistic",y = NULL,fill = "FDR") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 11),
        axis.text.y = element_text(size = 9,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black"),
        legend.title = element_text(size = 10),legend.text = element_text(size = 9))
ggsave(file.path(result_dir, "SCN9A_continuous_GSEA/SCN9A_selected_GO_BP_GSEA_ridgeplot.pdf"),p,width = 11,height = 8)

## KEGG山峰图
kegg.selected.result <- kegg.gsea@result[match(kegg.selected.ids,kegg.gsea@result$ID),,drop = FALSE]
kegg.gsea.selected <- kegg.gsea
kegg.gsea.selected@result <- kegg.selected.result

p <- enrichplot::ridgeplot(kegg.gsea.selected,showCategory = nrow(kegg.gsea.selected@result),
                           fill = "p.adjust",core_enrichment = TRUE,label_format = 42,
                           orderBy = "NES",decreasing = FALSE) +
  labs(title = "Selected KEGG pathways associated with SCN9A",
       subtitle = "Positive and negative GSEA enrichment",
       x = "DESeq2 Wald statistic",y = NULL,fill = "FDR") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,size = 15,face = "bold"),
        plot.subtitle = element_text(hjust = 0.5,size = 11),
        axis.text.y = element_text(size = 9,colour = "black"),
        axis.text.x = element_text(size = 10,colour = "black"),
        legend.title = element_text(size = 10),legend.text = element_text(size = 9))

ggsave(file.path(result_dir, "SCN9A_continuous_GSEA/SCN9A_selected_KEGG_GSEA_ridgeplot.pdf"),p,width = 10,height = 8)

## 输出Excel
openxlsx::write.xlsx(list(GO_BP_all = go.table, KEGG_all = kegg.table),file =file.path(result_dir, "SCN9A_continuous_GSEA/SCN9A_continuous_association_GO_KEGG_GSEA_results.xlsx"),overwrite = TRUE)

##### SCN9A expression vs predefined Ca2+ / regulated-exocytosis modules
### Define predefined GO modules
## Ca2+ module:
## GO:0098703 calcium ion import across plasma membrane
GO_CA <- "GO:0098703"

## Stimulus-secretion module:
## GO:0017156 calcium-ion regulated exocytosis
GO_EXOCYTOSIS <- "GO:0017156"

### Extract GO gene sets from org.Hs.eg.db
## GOALL而不是GO：
## GOALL包括gene对该GO term的直接和间接(parent term)注释，
## 更适合构建完整的GO module。
required.go.columns <- c("GOALL","ONTOLOGYALL")
go.annotation <- AnnotationDbi::select(org.Hs.eg.db,keys = rownames(vst.matrix),keytype = "SYMBOL",columns = c("GOALL","ONTOLOGYALL"))
go.annotation <- as.data.table(go.annotation)
go.annotation <- unique(go.annotation[!is.na(GOALL) & ONTOLOGYALL == "BP"])

### Get module genes
ca.genes <- unique(go.annotation[GOALL == GO_CA,SYMBOL])
exo.genes <- unique(go.annotation[GOALL == GO_EXOCYTOSIS,SYMBOL])

## 只保留当前beta-cell VST矩阵中实际检测到的genes
ca.genes.detected <- intersect(ca.genes,rownames(vst.matrix))
exo.genes.detected <- intersect(exo.genes,rownames(vst.matrix))

## 防止SCN9A进入module，避免self-correlation
ca.genes.detected <- setdiff(ca.genes.detected,"SCN9A")
exo.genes.detected <- setdiff(exo.genes.detected,"SCN9A")

## 两个module的overlap
module.overlap <- intersect(ca.genes.detected,exo.genes.detected)
cat("Overlap between modules:",length(module.overlap),"genes\n")

## Export predefined module gene lists
module.gene.table <- rbindlist(list(data.table(Module = "calcium ion import across plasma membrane",GO_ID = GO_CA,Gene = ca.genes.detected),
                                    data.table(Module = "calcium-ion regulated exocytosis",GO_ID = GO_EXOCYTOSIS,Gene = exo.genes.detected)))
module.list <- list(Ca2_import_across_plasma_membrane = ca.genes.detected,Ca2_regulated_exocytosis =exo.genes.detected)

### PRIMARY ANALYSIS: GSVA
gsva.param <- GSVA::gsvaParam(exprData = vst.matrix,geneSets = module.list,kcdf = "Gaussian",minSize = 10,maxSize = Inf)
gsva.score <- GSVA::gsva(gsva.param,verbose = FALSE)
gsva.score <- as.matrix(gsva.score)
gsva.score <- gsva.score[,colnames(vst.matrix),drop = FALSE]

### 8. Construct donor-level analysis dataframe
module.df <- data.frame(PB_ID = colnames(vst.matrix),
                        SCN9A_expression = pseudobulk.metadata[colnames(vst.matrix),"SCN9A_PB_log2"],
                        Ca2_GSVA =as.numeric(gsva.score["Ca2_import_across_plasma_membrane",colnames(vst.matrix)]),
                        Exocytosis_GSVA = as.numeric(gsva.score["Ca2_regulated_exocytosis",colnames(vst.matrix)]),
                        donor_accession = pseudobulk.metadata[colnames(vst.matrix),"donor_accession"], 
                        stringsAsFactors = FALSE)

## 附加其他metadata
for(current.variable in c("source","chemistry","sex","age","diabetes_group","hba1c_group","SCN9A_percent_detected")){
  module.df[[current.variable]] <- pseudobulk.metadata[module.df$PB_ID,current.variable]
}

### 9. Primary Spearman correlations
cor.ca <- cor.test(module.df$SCN9A_expression,module.df$Ca2_GSVA,method = "spearman",exact = FALSE)
cor.exo <- cor.test(module.df$SCN9A_expression,module.df$Exocytosis_GSVA,method = "spearman",exact = FALSE)
print(cor.ca)
print(cor.exo)

## Summarize correlation statistics
correlation.result <- data.table(Module = c("Ca2+ influx","Ca2+-regulated exocytosis"),
                                 GO_ID = c(GO_CA,GO_EXOCYTOSIS),
                                 N = c(sum(complete.cases(module.df$SCN9A_expression,module.df$Ca2_GSVA)),
                                       sum(complete.cases(module.df$SCN9A_expression,module.df$Exocytosis_GSVA))),
                                 Spearman_rho = c(unname(cor.ca$estimate),unname(cor.exo$estimate)),
                                 P_value = c(cor.ca$p.value,cor.exo$p.value))

## 因为两个primary hypotheses，同时输出BH-adjusted P
correlation.result[,FDR := p.adjust(P_value,method = "BH")]

### 相关系数散点图
## SCN9A vs Ca2+ GSVA
rho.ca <- unname(cor.ca$estimate)
pvalue.ca <- ifelse(cor.ca$p.value < 0.001,format(cor.ca$p.value,scientific = TRUE,digits = 3),sprintf("%.3f",cor.ca$p.value))

p.ca <- ggplot(module.df,aes(x = SCN9A_expression,y = Ca2_GSVA,colour = hba1c_group)) +
  geom_point(size = 3,alpha = 0.85) +
  geom_smooth(aes(group = 1),method = "lm",formula = y ~ x,se = TRUE,colour = "black",linewidth = 0.7) +
  scale_colour_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")) +
  # annotate("text",x = -Inf,y = Inf,label = label.ca,hjust = -0.08,vjust = 1.15,size = 3.5)+
  labs(title = "SCN9A and Ca²⁺ influx",
       subtitle = paste0("Donor-level Spearman \u03c1 = ",round(rho.ca, 3),"; P = ",pvalue.ca),
       x = paste0("\u03b2-cell SCN9A expression\n","log2(DESeq2 normalized count + 1)"),
       y = "Ca²⁺ influx GSVA score",
       colour = "HbA1c group") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"))
ggsave(file.path(result_dir, "SCN9A_vs_Ca2_module_GSVA.pdf"),p.ca,width = 6,height = 5,device = cairo_pdf)

## SCN9A vs regulated exocytosis GSVA
rho.exo <- unname(cor.exo$estimate)
pvalue.exo <- ifelse(cor.exo$p.value < 0.001,format(cor.exo$p.value,scientific = TRUE,digits = 3),sprintf("%.3f",cor.exo$p.value))

p.exo <- ggplot(module.df,aes(x = SCN9A_expression,y = Exocytosis_GSVA,colour = hba1c_group)) +
  geom_point(size = 3,alpha = 0.85) +
  geom_smooth(aes(group = 1),method = "lm",formula = y ~ x,se = TRUE,colour = "black",linewidth = 0.7) +
  scale_colour_manual(values = c("Normal" = "#4DBBD5","Prediabetes" = "#EFC000","Diabetes" = "#E64B35")) +
  labs(title = "SCN9A and Ca²⁺-regulated exocytosis",
       subtitle = paste0("Donor-level Spearman \u03c1 = ",round(rho.exo, 3),"; P = ",pvalue.exo),
       x = paste0("\u03b2-cell SCN9A expression\n","log2(DESeq2 normalized count + 1)"),
       y = "Ca²⁺-regulated exocytosis\nGSVA score",
       colour = "HbA1c group") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(colour = "black"))
ggsave(file.path(result_dir, "SCN9A_vs_regulated_exocytosis_module_GSVA.pdf"),p.exo,width = 6,height = 5,device = cairo_pdf)

### SENSITIVITY ANALYSIS:
## mean gene-wise Z-score
## 目的：检查结论是否依赖GSVA算法
## 每个gene在donor之间Z-score
gene.z <- t(scale(t(vst.matrix)))

## 删除zero variance等导致的non-finite genes
valid.z.genes <- rownames(gene.z)[apply(gene.z,1,function(x){all(is.finite(x))})]
gene.z <- gene.z[valid.z.genes,,drop = FALSE]

ca.z.genes <- intersect(ca.genes.detected,rownames(gene.z))
exo.z.genes <- intersect(exo.genes.detected,rownames(gene.z))

module.df$Ca2_meanZ <- colMeans(gene.z[ca.z.genes,module.df$PB_ID,drop = FALSE],na.rm = TRUE)
module.df$Exocytosis_meanZ <- colMeans(gene.z[exo.z.genes,module.df$PB_ID,drop = FALSE],na.rm = TRUE)

## mean-Z correlations
cor.ca.z <- cor.test(module.df$SCN9A_expression,module.df$Ca2_meanZ,method = "spearman",exact = FALSE)
cor.exo.z <- cor.test(module.df$SCN9A_expression,module.df$Exocytosis_meanZ,method = "spearman",exact = FALSE)

sensitivity.result <- data.table(Module = c("Ca2+ influx","Ca2+-regulated exocytosis"),
                                 Method = "Mean gene-wise Z-score",
                                 Spearman_rho = c(unname(cor.ca.z$estimate),unname(cor.exo.z$estimate)),
                                 P_value = c(cor.ca.z$p.value,cor.exo.z$p.value))
sensitivity.result[,FDR :=p.adjust(P_value,method = "BH")]
cat("\n==============================\n","SENSITIVITY RESULTS: MEAN-Z\n","==============================\n")
print(sensitivity.result)

### 输出module基因文件
openxlsx::write.xlsx(list(Module_gene_sets =module.gene.table),
                     file = file.path(result_dir, "SCN9A_functional_GO_module_gene_mapping.xlsx"),
                     overwrite = TRUE,rowNames = FALSE)

