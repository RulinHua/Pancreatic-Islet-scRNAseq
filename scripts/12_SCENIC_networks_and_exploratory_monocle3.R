# SCENIC direct-target/network analysis and exploratory Monocle3 code
# ===================================================================
# Run scripts/00_setup.R first in the same R session.
# This file was reorganized from the original analysis script without intentionally
# changing the statistical analysis logic. Dataset-specific results may evolve.

##### 回到 R中读取 regulon 靶基因
regulon.targets <- fread(file.path(result_dir, "SCENIC/HPAP_V3_regulon_targets.tsv"))
regulon.size <- regulon.targets[,.(regulon_size = uniqueN(target)),by = .(regulon,TF)]
setorder(regulon.size,-regulon_size)

##### 从 adjacency 提取 TF–SCN9A 关系
adjacency <- fread(file.path(result_dir, "SCENIC/HPAP_V3_adjacencies.tsv"))
setnames(adjacency,names(adjacency),tolower(names(adjacency)))

### GRNBoost2 共表达候选
SCN9A.adjacency <- adjacency[target == "SCN9A"]
setorder(SCN9A.adjacency,-importance)
fwrite(SCN9A.adjacency,file.path(result_dir, "SCN9A_all_GRNBoost2_upstream_TFs.tsv"),sep = "\t")

##### motif 剪枝后仍包含 SCN9A 的 TF
SCN9A.regulons <- regulon.targets[target == "SCN9A"]
SCN9A.regulons
### 合并 adjacency importance
SCN9A.direct <- merge(SCN9A.regulons[,.(regulon,TF,target_weight = weight)],SCN9A.adjacency[,.(TF = tf,adjacency_importance = importance)],by = "TF",all.x = TRUE)
### 再合并差异 regulon 结果
SCN9A.direct <- merge(SCN9A.direct,result.combined,by = "regulon",all.x = TRUE)
### 加入 regulon 大小
SCN9A.direct <- merge(SCN9A.direct,regulon.size,by = c("regulon","TF"),all.x = TRUE)
### 分类
SCN9A.direct[,evidence_class := "motif-supported candidate"]
SCN9A.direct[!is.na(FDR_HPAP) & FDR_HPAP < 0.05,evidence_class := "HPAP-supported direct candidate"]
SCN9A.direct[replicated == TRUE,evidence_class := "HPAP-IIDP replicated direct candidate"]

setorder(SCN9A.direct,FDR_meta,FDR_HPAP,-adjacency_importance)
fwrite(SCN9A.direct,file.path(result_dir, "SCN9A_motif_supported_upstream_TFs.tsv"),sep = "\t")

##### regulon的Cytoscape网络图
### Top10 regulon Cytoscape网络数据
network.output.dir <- file.path(result_dir, "Cytoscape_Top10_Regulon_Networks")

## 根据GO注释定义黏附相关基因
adhesion.root.go <- "GO:0007155"
all.go.annotation <- AnnotationDbi::select(org.Hs.eg.db,
                                           keys = AnnotationDbi::keys(org.Hs.eg.db,keytype = "ENTREZID"),
                                           keytype = "ENTREZID",
                                           columns = c("SYMBOL","GOALL","ONTOLOGYALL"))
all.go.annotation <- as.data.table(all.go.annotation)
table(all.go.annotation$GOALL == adhesion.root.go,useNA = "ifany") ## 检查GO:0007155是否存在

## 提取cell adhesion相关基因
adhesion.annotation <- all.go.annotation[GOALL == adhesion.root.go & ONTOLOGYALL == "BP"]
## 提取gene symbol
adhesion.genes <- unique(adhesion.annotation$SYMBOL[!is.na(adhesion.annotation$SYMBOL) & adhesion.annotation$SYMBOL != ""])

## 保存每个网络的数据
network.data.list <- list()
network.summary <- list()

## 每个网络显示的普通靶基因数量
number.total.targets <- 30
number.adhesion.targets <- 10
for(regulon.name in top10.regulons.hpap){
  ## 当前regulon的全部SCENIC靶基因
  current.targets <- copy(regulon.targets[regulon == regulon.name & !is.na(target) &  target != ""])
  if(nrow(current.targets) == 0){
    warning(paste0("No targets found for ",regulon.name))
    next
  }
  
  ## TF名称
  tf.name <- unique(current.targets$TF)
  
  ## 删除重复靶基因和TF自身
  current.targets <- current.targets[target != tf.name]
  current.targets <- current.targets[order(-weight,na.last = TRUE)]
  current.targets <- unique(current.targets,by = "target")
  
  ## 判断SCN9A和黏附相关靶基因
  current.targets[,is_SCN9A := target == "SCN9A"]
  current.targets[,is_adhesion := target %in% adhesion.genes]
  current.adhesion.targets <- current.targets[is_adhesion == TRUE,target]
  SCN9A.present <- any(current.targets$is_SCN9A)
  
  ## 按SCENIC weight选择黏附相关靶基因
  selected.adhesion.targets <- head(current.targets[is_adhesion == TRUE][order(-weight,na.last = TRUE)]$target,number.adhesion.targets)
  
  ## SCN9A如果属于当前regulon，则强制保留
  forced.targets <- unique(c(if(SCN9A.present) "SCN9A",selected.adhesion.targets))
  
  ## 计算还可以加入多少普通靶基因
  number.other.targets <- max(0,number.total.targets - length(forced.targets))
  
  ## 按weight选择剩余的普通靶基因
  selected.other.targets <- head(current.targets[!target %in% forced.targets][order(-weight,na.last = TRUE)]$target,number.other.targets)
  
  ## 最终用于绘图的靶基因最多30个
  selected.target.names <- unique(c(forced.targets,selected.other.targets))
  selected.targets <- current.targets[target %in% selected.target.names]
  
  ## 靶基因分类
  selected.targets[,node_class := fcase(is_SCN9A,"SCN9A",is_adhesion,"Adhesion_related",default = "Other_target")]
  
  ## Edge表
  edge.table <- selected.targets[,.(source = tf.name,target,interaction = "motif_supported",weight,edge_class = node_class)]
  
  ## 处理缺失权重
  edge.table[,edge_weight := fifelse(is.finite(weight),weight,0)]
  
  ## 将边宽缩放至1—5
  if(length(unique(edge.table$edge_weight)) > 1){
    edge.table[,edge_width := 1 + 4 * (edge_weight - min(edge_weight)) / (max(edge_weight) - min(edge_weight))]
  } else {
    edge.table[,edge_width := 2.5]
  }
  
  ## Node表
  node.table <- rbind(data.table(id = tf.name,display_label = tf.name,node_class = "TF",node_size = 55),
                      selected.targets[,.(id = target,
                                          display_label = target,
                                          node_class,
                                          node_size = fcase(node_class == "SCN9A",45,node_class == "Adhesion_related",36,default = 24))]
  )
  node.table <- unique(node.table,by = "id")
  
  ## 保存Cytoscape导入文件
  safe.regulon.name <- gsub("[^A-Za-z0-9_-]+","_",regulon.name)
  fwrite(node.table,file.path(network.output.dir,paste0(safe.regulon.name,"_nodes.tsv")),sep = "\t")
  fwrite(edge.table,file.path(network.output.dir,paste0(safe.regulon.name,"_edges.tsv")),sep = "\t")
  
  ## 保存GraphML文件
  graph.object <- graph_from_data_frame(d = as.data.frame(edge.table[,.(source,target,interaction,weight = edge_weight,edge_width,edge_class)]),directed = TRUE,vertices = as.data.frame(node.table))
  write_graph(graph.object,file.path(network.output.dir,paste0(safe.regulon.name,"_network.graphml")),format = "graphml")
  
  ## 保留对象，后续RCy3直接画图
  network.data.list[[regulon.name]] <- list(nodes = as.data.frame(node.table),edges = as.data.frame(edge.table))
  network.summary[[regulon.name]] <- data.table(regulon = regulon.name,TF = tf.name,
                                                total_SCENIC_targets = nrow(current.targets),
                                                displayed_targets = nrow(selected.targets),
                                                SCN9A_is_target = SCN9A.present,
                                                total_adhesion_targets = length(current.adhesion.targets),
                                                displayed_adhesion_targets = length(selected.adhesion.targets),
                                                displayed_adhesion_genes = paste(selected.adhesion.targets,collapse = ";"))
}

### 汇总结果
network.summary <- rbindlist(network.summary,fill = TRUE)
network.summary
fwrite(network.summary,file.path(network.output.dir,"Top10_regulon_network_summary.tsv"),sep = "\t")

### 绘制Top10 regulon网络图
for(regulon.name in names(network.data.list)){
  ## 当前网络的node和edge表
  node.table <- as.data.table(network.data.list[[regulon.name]]$nodes)
  edge.table <- as.data.table(network.data.list[[regulon.name]]$edges)
  
  ## 所有节点都显示基因名称
  node.table[,display_label := id]
  
  ## 当前regulon的TF
  tf.name <- node.table[node_class == "TF",id][1]
  
  ## 构建igraph网络
  graph.object <- graph_from_data_frame(d = as.data.frame(edge.table[,.(source,target,edge_width,edge_class)]),directed = TRUE,vertices = as.data.frame(node.table))
  
  ## 找到TF节点
  tf.index <- which(V(graph.object)$name == tf.name)
  
  ## 绘图
  p <- ggraph(graph.object,layout = "star",center = tf.index) +
    ## TF → target的边
    geom_edge_link(aes(width = edge_width,colour = edge_class),alpha = 0.70,
                   arrow = arrow(length = unit(2.5,"mm"),type = "closed"),
                   end_cap = circle(4,"mm")) +
    ## 节点
    geom_node_point(aes(size = node_size,fill = node_class),shape = 21,colour = "black",stroke = 0.6) +
    ## 所有基因名称
    geom_node_text(aes(label = display_label),repel = TRUE,size = 4.5,fontface = "bold") +
    ## 节点大小
    # scale_size_manual(values = c(TF = 12,SCN9A = 10,Adhesion_related = 8,Other_target = 6)) +
    scale_size_identity() +
    ## 节点颜色
    scale_fill_manual(values = c(TF = "#E64B35",SCN9A = "#DC0000",Adhesion_related = "#4DBBD5",Other_target = "grey80")) +
    ## 边颜色
    scale_edge_colour_manual(values = c(SCN9A = "#DC0000",Adhesion_related = "#4DBBD5",Other_target = "grey65")) +
    ## edge_width原来是1-5，这里映射成合适的实际绘图线宽
    scale_edge_width_continuous(range = c(0.5,2.5)) +
    ggtitle(regulon.name) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),legend.position = "right")
  
  ## 安全文件名
  safe.regulon.name <- gsub("[^A-Za-z0-9_-]+","_",regulon.name)
  
  ## 保存PDF
  ggsave(file.path(network.output.dir,paste0(safe.regulon.name,"_network.pdf")),p,width = 11,height = 10,units = "in")
}

# ########## monocle3
# beta.hpap.v3 <- sr[,sr$Cell_Type=="Beta" & sr$source == "HPAP" & sr$chemistry == "V3"]
# beta.hpap.v3 <- NormalizeData(beta.hpap.v3,assay = "RNA",normalization.method = "LogNormalize",scale.factor = 10000,verbose = FALSE)
# beta.hpap.v3 <- FindVariableFeatures(beta.hpap.v3,selection.method = "vst",nfeatures = 2000,verbose = FALSE)
# trajectory.genes <- VariableFeatures(beta.hpap.v3)
# excluded.genes <- unique(c("MALAT1",grep("^MT-|^RPL|^RPS",trajectory.genes,value = TRUE)))
# hemoglobin.genes <- c("HBA1","HBA2","HBB","HBD","HBG1","HBG2")
# excluded.genes <- c(excluded.genes,hemoglobin.genes)
# trajectory.genes <- setdiff(trajectory.genes,excluded.genes)
# 
# counts.hpap.v3 <- GetAssayData(beta.hpap.v3,assay = "RNA",slot = "counts")
# keep.genes <- Matrix::rowSums(counts.hpap.v3) > 0
# counts.hpap.v3 <- counts.hpap.v3[keep.genes,,drop = FALSE]
# 
# cell.metadata <- beta.hpap.v3@meta.data[colnames(counts.hpap.v3),,drop = FALSE]
# gene.metadata <- data.frame(gene_short_name = rownames(counts.hpap.v3),row.names = rownames(counts.hpap.v3))
# 
# cds.hpap.v3 <- new_cell_data_set(expression_data = counts.hpap.v3,cell_metadata = cell.metadata,gene_metadata = gene.metadata)
# cds.hpap.v3 <- estimate_size_factors(cds.hpap.v3)
# cds.hpap.v3 <- detect_genes(cds.hpap.v3,min_expr = 0.1)
# trajectory.genes <- intersect(trajectory.genes,rownames(cds.hpap.v3))
# length(trajectory.genes)
# 
# set.seed(123)
# cds.hpap.v3 <- preprocess_cds(
#   cds.hpap.v3,
#   method = "PCA",
#   num_dim = 50,
#   norm_method = "log",
#   use_genes = trajectory.genes,
#   scaling = TRUE,
#   verbose = TRUE
# )
# plot_pc_variance_explained(cds.hpap.v3)
# cds.hpap.v3 <- preprocess_cds(
#   cds.hpap.v3,
#   method = "PCA",
#   num_dim = 30,
#   norm_method = "log",
#   use_genes = trajectory.genes,
#   scaling = TRUE,
#   verbose = TRUE
# )
# set.seed(123)
# 
# cds.hpap.v3 <- reduce_dimension(
#   cds.hpap.v3,
#   reduction_method = "UMAP",
#   preprocess_method = "PCA",
#   umap.metric = "cosine",
#   umap.n_neighbors = 30,
#   umap.min_dist = 0.1,
#   umap.fast_sgd = FALSE,
#   cores = 20,
#   verbose = TRUE
# )
# cds.hpap.v3 <- cluster_cells(
#   cds.hpap.v3,
#   reduction_method = "UMAP",
#   cluster_method = "leiden",
#   k = 20,
#   partition_qval = 0.05,
#   random_seed = 123,
#   verbose = TRUE
# )
# 
# colData(cds.hpap.v3)$monocle3_cluster <- as.character(
#   cluster.vector
# )
# 
# colData(cds.hpap.v3)$monocle3_partition <- as.character(
#   partition.vector
# )
# 
# p.cluster <- plot_cells(
#   cds.hpap.v3,
#   color_cells_by = "monocle3_cluster",
#   show_trajectory_graph = FALSE,
#   label_cell_groups = TRUE,
#   label_leaves = FALSE,
#   label_branch_points = FALSE,
#   cell_size = 0.08
# )
# ggsave(file.path(result_dir, "Beta_cell_Monocle3/HPAP_V3_Beta_UMAP_Cluster.png")
#   ,
#   p.cluster,
#   width = 7,
#   height = 6
# )
# p.partition <- plot_cells(
#   cds.hpap.v3,
#   color_cells_by = "monocle3_partition",
#   show_trajectory_graph = FALSE,
#   label_cell_groups = TRUE,
#   label_leaves = FALSE,
#   label_branch_points = FALSE,
#   cell_size = 0.08
# )
# ggsave(file.path(result_dir, "Beta_cell_Monocle3/HPAP_V3_Beta_UMAP_Partition.png"),
#   p.partition,
#   width = 5,
#   height = 4)
# 
# p.donor <- plot_cells(
#   cds.hpap.v3,
#   color_cells_by = "donor_accession",
#   show_trajectory_graph = FALSE,
#   label_cell_groups = FALSE,
#   label_leaves = FALSE,
#   label_branch_points = FALSE,
#   cell_size = 0.05
# ) +
#   theme(aspect.ratio=1,legend.position = "none")
# ggsave(file.path(result_dir, "Beta_cell_Monocle3/HPAP_V3_Beta_UMAP_Donor.png"),
#   p.donor,
#   width = 5,
#   height = 4
# )
# 
# p.hba1c <- plot_cells(
#   cds.hpap.v3,
#   color_cells_by = "hba1c_group",
#   show_trajectory_graph = FALSE,
#   label_cell_groups = FALSE,
#   label_leaves = FALSE,
#   label_branch_points = FALSE,
#   cell_size = 0.08
# ) +
#   scale_color_manual(
#     values = c(
#       "Normal" = "#4DBBD5",
#       "Prediabetes" = "#EFC000",
#       "Diabetes" = "#E64B35"
#     ),
#     na.value = "grey80"
#   )
# ggsave(file.path(result_dir, "Beta_cell_Monocle3/HPAP_V3_Beta_UMAP_HbA1cGroup.png"),p.hba1c,width = 5,height = 4)       

