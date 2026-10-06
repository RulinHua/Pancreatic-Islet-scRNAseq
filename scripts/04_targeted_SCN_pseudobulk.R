# Targeted SCN-family pseudobulk analysis
# =======================================
# Run scripts/00_setup.R first in the same R session.
#
# This module complements the global cell-type pseudobulk analysis in script 03.
# It focuses on SCN-family genes and TRPC4, uses donor-level biological replicates,
# and adds covariate-aware models for diabetes-status comparisons.

minimum.cells <- 1
minimum.donors <- 3
scn_genes <- c(
  "SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A",
  "SCN7A", "SCN8A", "SCN9A", "SCN10A", "SCN11A", "TRPC4"
)

outdir <- file.path(result_dir, "SCN_pseudobulk_DESeq2_diabetes_group")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Donor metadata and donor x cell-type pseudobulk counts
# -----------------------------------------------------------------------------
meta.dt <- as.data.table(sr@meta.data, keep.rownames = "Cell")
meta.dt[, `:=`(
  Cell = as.character(Cell),
  donor_accession = as.character(donor_accession),
  Cell_Type = as.character(Cell_Type),
  diabetes_group = as.character(diabetes_group),
  hba1c_group = as.character(hba1c_group),
  source = as.character(source),
  chemistry = as.character(chemistry),
  sex = as.character(sex),
  age = as.numeric(as.character(`age_(years)`))
)]

for (v in c("donor_accession", "Cell_Type", "diabetes_group", "hba1c_group", "source", "chemistry", "sex")) {
  set(meta.dt, which(!is.na(meta.dt[[v]]) & trimws(meta.dt[[v]]) == ""), v, NA_character_)
}

donor.info <- unique(meta.dt[, .(
  donor_accession, diabetes_group, hba1c_group,
  source, chemistry, sex, age
)])

duplicated.donor <- donor.info[, .N, by = donor_accession][N > 1]
if (nrow(duplicated.donor) > 0) {
  warning("Some donors have more than one donor-level metadata record; inspect duplicated.donor.")
}

# Complete covariate information is required for the covariate-adjusted diabetes model.
donor.info <- donor.info[
  complete.cases(donor.info[, .(diabetes_group, source, chemistry, sex, age)])
]

donor.info[, diabetes_group := factor(diabetes_group, levels = c("Control", "T1D", "T2D"))]
donor.info[, hba1c_group := factor(hba1c_group, levels = c("Normal", "Prediabetes", "Diabetes"))]
donor.info[, source := factor(source)]
donor.info[, chemistry := factor(chemistry)]
donor.info[, sex := factor(sex)]

cell.number <- meta.dt[
  donor_accession %in% donor.info$donor_accession,
  .(n_cells = .N),
  by = .(donor_accession, Cell_Type)
]
cell.number <- merge(cell.number, donor.info, by = "donor_accession", all.x = TRUE, sort = FALSE)
cell.number <- cell.number[n_cells >= minimum.cells]

# Require at least three donors per disease group and at least two groups per cell type.
donor.number.table <- cell.number[, .(n_donor = uniqueN(donor_accession)), by = .(Cell_Type, diabetes_group)]
valid.group <- donor.number.table[n_donor >= minimum.donors]
cell.number <- merge(
  cell.number,
  valid.group[, .(Cell_Type, diabetes_group)],
  by = c("Cell_Type", "diabetes_group"),
  all = FALSE,
  sort = FALSE
)
valid.celltype <- cell.number[, .(n_group = uniqueN(diabetes_group)), by = Cell_Type][n_group >= 2]
cell.number <- merge(cell.number, valid.celltype[, .(Cell_Type)], by = "Cell_Type", all = FALSE, sort = FALSE)

setorder(cell.number, Cell_Type, diabetes_group, donor_accession)
cell.number[, PB_ID := sprintf("PB%04d", seq_len(.N))]

meta.use <- merge(
  meta.dt[, .(Cell, donor_accession, Cell_Type)],
  cell.number[, .(donor_accession, Cell_Type, PB_ID)],
  by = c("donor_accession", "Cell_Type"),
  all = FALSE,
  sort = FALSE
)

sr.use <- subset(sr, cells = meta.use$Cell)
pb.id.vector <- setNames(meta.use$PB_ID, meta.use$Cell)
sr.use$PB_ID <- unname(pb.id.vector[colnames(sr.use)])
stopifnot(!any(is.na(sr.use$PB_ID)))

pseudobulk.counts <- AggregateExpression(
  object = sr.use,
  assays = "RNA",
  group.by = "PB_ID",
  slot = "counts",
  return.seurat = FALSE,
  verbose = FALSE
)$RNA

pseudobulk.metadata <- as.data.frame(cell.number)
pseudobulk.metadata <- pseudobulk.metadata[
  match(colnames(pseudobulk.counts), pseudobulk.metadata$PB_ID),
  , drop = FALSE
]
rownames(pseudobulk.metadata) <- pseudobulk.metadata$PB_ID
stopifnot(identical(colnames(pseudobulk.counts), rownames(pseudobulk.metadata)))

# -----------------------------------------------------------------------------
# 2. Covariate-aware SCN/TRPC4 DESeq2 models by cell type
# -----------------------------------------------------------------------------
scn.result.list <- list()
scn.expression.list <- list()
celltypes.use <- sort(unique(as.character(pseudobulk.metadata$Cell_Type)))
comparisons.all <- list(
  c("T1D", "Control"),
  c("T2D", "Control"),
  c("T2D", "T1D")
)

for (ct in celltypes.use) {
  message("Running targeted pseudobulk model: ", ct)

  sample.ids <- rownames(pseudobulk.metadata)[pseudobulk.metadata$Cell_Type == ct]
  metadata.ct <- pseudobulk.metadata[sample.ids, , drop = FALSE]
  counts.ct <- as.matrix(pseudobulk.counts[, sample.ids, drop = FALSE])
  storage.mode(counts.ct) <- "integer"

  group.order <- c("Control", "T1D", "T2D")
  group.order <- group.order[group.order %in% unique(as.character(metadata.ct$diabetes_group))]
  metadata.ct$diabetes_group <- factor(as.character(metadata.ct$diabetes_group), levels = group.order)
  metadata.ct$source <- droplevels(factor(metadata.ct$source))
  metadata.ct$chemistry <- droplevels(factor(metadata.ct$chemistry))
  metadata.ct$sex <- droplevels(factor(metadata.ct$sex))
  metadata.ct$age_scaled <- as.numeric(scale(metadata.ct$age))

  # Retain categorical covariates only when they contribute independent information.
  candidate.variables <- c("source", "chemistry", "sex")
  variables.use <- candidate.variables[vapply(
    candidate.variables,
    function(v) nlevels(metadata.ct[[v]]) > 1,
    logical(1)
  )]

  confounded.variables <- character(0)
  for (v in variables.use) {
    pair.matrix <- model.matrix(reformulate(c(v, "diabetes_group")), data = metadata.ct)
    if (qr(pair.matrix)$rank < ncol(pair.matrix)) {
      confounded.variables <- c(confounded.variables, v)
    }
  }

  if (length(confounded.variables) > 0) {
    message(
      "Skipped ", ct, ": diabetes_group is completely confounded with ",
      paste(confounded.variables, collapse = ", ")
    )
    next
  }

  selected.variables <- character(0)
  current.matrix <- model.matrix(reformulate("diabetes_group"), data = metadata.ct)

  for (v in variables.use) {
    trial.formula <- reformulate(c(selected.variables, v, "diabetes_group"))
    trial.matrix <- model.matrix(trial.formula, data = metadata.ct)
    expected.df <- nlevels(metadata.ct[[v]]) - 1
    added.df <- qr(trial.matrix)$rank - qr(current.matrix)$rank

    if (added.df == expected.df) {
      selected.variables <- c(selected.variables, v)
      current.matrix <- trial.matrix
    } else {
      message(ct, ": removed redundant covariate ", v)
    }
  }

  design.terms <- selected.variables
  if (sum(is.finite(metadata.ct$age_scaled)) > 1 && sd(metadata.ct$age_scaled, na.rm = TRUE) > 0) {
    design.terms <- c(design.terms, "age_scaled")
  }
  design.terms <- c(design.terms, "diabetes_group")
  design.formula <- reformulate(design.terms)
  design.matrix <- model.matrix(design.formula, data = metadata.ct)

  if (qr(design.matrix)$rank < ncol(design.matrix)) {
    stop(
      "Rank-deficient final design for cell type ", ct,
      ": ", paste(deparse(design.formula), collapse = "")
    )
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts.ct,
    colData = metadata.ct,
    design = design.formula
  )
  dds <- DESeq(dds, quiet = TRUE)

  scn.present <- intersect(scn_genes, rownames(counts.ct))
  scn.normalized <- counts(dds, normalized = TRUE)[scn.present, , drop = FALSE]
  scn.log2.normalized <- log2(scn.normalized + 1)

  expression.dt <- as.data.table(t(scn.log2.normalized), keep.rownames = "PB_ID")
  expression.dt <- melt(
    expression.dt,
    id.vars = "PB_ID",
    measure.vars = scn.present,
    variable.name = "Gene",
    value.name = "log2_normalized_count",
    variable.factor = FALSE
  )
  metadata.ct.dt <- as.data.table(metadata.ct, keep.rownames = "PB_ID")
  expression.dt <- merge(expression.dt, metadata.ct.dt, by = "PB_ID", all.x = TRUE, sort = FALSE)
  scn.expression.list[[ct]] <- expression.dt

  for (comparison in comparisons.all) {
    numerator <- comparison[1]
    denominator <- comparison[2]
    if (!all(c(numerator, denominator) %in% group.order)) next

    res <- results(
      dds,
      contrast = c("diabetes_group", numerator, denominator),
      alpha = 0.05,
      pAdjustMethod = "BH"
    )
    res.dt <- as.data.table(as.data.frame(res), keep.rownames = "Gene")
    res.dt[, `:=`(
      Cell_Type = ct,
      Contrast = paste0(numerator, "_vs_", denominator),
      Numerator = numerator,
      Denominator = denominator,
      Design = paste(deparse(design.formula), collapse = "")
    )]
    setcolorder(
      res.dt,
      c(
        "Gene", "Cell_Type", "Contrast", "Numerator", "Denominator", "Design",
        "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
      )
    )

    safe.celltype <- gsub("[^A-Za-z0-9]+", "_", ct)
    fwrite(
      res.dt,
      file.path(outdir, paste0(safe.celltype, "_", numerator, "_vs_", denominator, "_DESeq2_all_genes.csv"))
    )

    scn.res.dt <- merge(data.table(Gene = scn_genes), res.dt, by = "Gene", all.x = TRUE, sort = FALSE)
    result.name <- paste(ct, numerator, denominator, sep = "__")
    scn.result.list[[result.name]] <- scn.res.dt
  }
}

scn.result.table <- rbindlist(scn.result.list, use.names = TRUE, fill = TRUE, idcol = "Result_ID")
fwrite(
  scn.result.table,
  file.path(outdir, "SCN_TRPC4_DESeq2_pseudobulk_all_celltypes_all_contrasts.csv")
)

scn.expression.table <- rbindlist(scn.expression.list, use.names = TRUE, fill = TRUE)
fwrite(
  scn.expression.table,
  file.path(outdir, "SCN_TRPC4_donor_log2_normalized_counts_all_celltypes.csv")
)

# -----------------------------------------------------------------------------
# 3. Donor-level beta-cell SCN9A and alpha-cell TRPC4 expression tables
# -----------------------------------------------------------------------------
pseudobulk.metadata.all <- pseudobulk.metadata

beta.ids <- rownames(pseudobulk.metadata.all)[pseudobulk.metadata.all$Cell_Type == "Beta"]
beta.metadata <- pseudobulk.metadata.all[beta.ids, , drop = FALSE]
beta.counts <- as.matrix(pseudobulk.counts[, beta.ids, drop = FALSE])
stopifnot(identical(colnames(beta.counts), rownames(beta.metadata)))

dds.beta <- DESeqDataSetFromMatrix(beta.counts, beta.metadata, design = ~ 1)
dds.beta <- estimateSizeFactors(dds.beta, type = "ratio")
beta.normalized <- counts(dds.beta, normalized = TRUE)

beta.metadata$Beta_total_pseudobulk_UMI <- colSums(beta.counts)
beta.metadata$Beta_size_factor <- sizeFactors(dds.beta)
beta.metadata$SCN9A_Beta_raw_count <- beta.counts["SCN9A", ]
beta.metadata$SCN9A_Beta_normalized_count <- beta.normalized["SCN9A", ]
beta.metadata$SCN9A_Beta_log2_normalized_count <- log2(beta.normalized["SCN9A", ] + 1)
beta.metadata$SCN9A_Beta_detected <- beta.counts["SCN9A", ] > 0

alpha.ids <- rownames(pseudobulk.metadata.all)[pseudobulk.metadata.all$Cell_Type == "Alpha"]
alpha.metadata <- pseudobulk.metadata.all[alpha.ids, , drop = FALSE]
alpha.counts <- as.matrix(pseudobulk.counts[, alpha.ids, drop = FALSE])
stopifnot(identical(colnames(alpha.counts), rownames(alpha.metadata)))

dds.alpha <- DESeqDataSetFromMatrix(alpha.counts, alpha.metadata, design = ~ 1)
dds.alpha <- estimateSizeFactors(dds.alpha, type = "ratio")
alpha.normalized <- counts(dds.alpha, normalized = TRUE)

alpha.metadata$Alpha_total_pseudobulk_UMI <- colSums(alpha.counts)
alpha.metadata$Alpha_size_factor <- sizeFactors(dds.alpha)
alpha.metadata$TRPC4_Alpha_raw_count <- alpha.counts["TRPC4", ]
alpha.metadata$TRPC4_Alpha_normalized_count <- alpha.normalized["TRPC4", ]
alpha.metadata$TRPC4_Alpha_log2_normalized_count <- log2(alpha.normalized["TRPC4", ] + 1)
alpha.metadata$TRPC4_Alpha_detected <- alpha.counts["TRPC4", ] > 0

openxlsx::write.xlsx(
  list(Beta_SCN9A = beta.metadata, Alpha_TRPC4 = alpha.metadata),
  file = file.path(result_dir, "PanKbase_donor_metadata_Beta_SCN9A_Alpha_TRPC4_pseudobulk.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

# -----------------------------------------------------------------------------
# 4. Focused beta-cell SCN9A comparisons and source-adjusted visualizations
# -----------------------------------------------------------------------------
# Diabetes-status model
beta.ids <- rownames(pseudobulk.metadata.all)[pseudobulk.metadata.all$Cell_Type == "Beta"]
beta.metadata <- pseudobulk.metadata.all[beta.ids, , drop = FALSE]
beta.counts <- as.matrix(pseudobulk.counts[, beta.ids, drop = FALSE])

dds.beta.diabetes <- DESeqDataSetFromMatrix(
  beta.counts,
  beta.metadata,
  design = ~ source + diabetes_group
)
dds.beta.diabetes <- DESeq(dds.beta.diabetes, quiet = TRUE)

beta.diabetes.results <- list()
for (comparison in comparisons.all) {
  numerator <- comparison[1]
  denominator <- comparison[2]
  res <- results(
    dds.beta.diabetes,
    contrast = c("diabetes_group", numerator, denominator),
    alpha = 0.05,
    pAdjustMethod = "BH"
  )
  res.dt <- as.data.table(as.data.frame(res), keep.rownames = "Gene")
  res.dt[, Contrast := paste0(numerator, "_vs_", denominator)]
  beta.diabetes.results[[paste(comparison, collapse = "_")]] <- res.dt
}
res.beta.diabetes <- rbindlist(beta.diabetes.results)
fwrite(
  res.beta.diabetes[Gene == "SCN9A"],
  file.path(outdir, "Beta_SCN9A_diabetes_status_DESeq2.csv")
)

vsd <- varianceStabilizingTransformation(dds.beta.diabetes, blind = FALSE)
vst.matrix <- assay(vsd)
group.design <- model.matrix(~ diabetes_group, data = as.data.frame(colData(dds.beta.diabetes)))
source.adjusted.matrix <- limma::removeBatchEffect(
  vst.matrix,
  batch = colData(dds.beta.diabetes)$source,
  design = group.design
)

plot.data <- data.table(
  PB_ID = rownames(beta.metadata),
  donor_accession = beta.metadata$donor_accession,
  source = beta.metadata$source,
  diabetes_group = beta.metadata$diabetes_group,
  n_cells = beta.metadata$n_cells,
  Expression = as.numeric(source.adjusted.matrix["SCN9A", rownames(beta.metadata)])
)
plot.data[, diabetes_group := factor(diabetes_group, levels = c("Control", "T1D", "T2D"))]

summary.data <- plot.data[, .(
  Mean = mean(Expression),
  SD = sd(Expression),
  n_donor = .N
), by = diabetes_group]
summary.data[, SEM := SD / sqrt(n_donor)]
summary.data[, `:=`(Lower = Mean - SEM, Upper = Mean + SEM)]

pvalue.data <- res.beta.diabetes[
  Gene == "SCN9A",
  .(Contrast, log2FoldChange, pvalue, padj)
]
pvalue.data[, c("group1", "group2") := tstrsplit(Contrast, "_vs_", fixed = TRUE)]
pvalue.data[, label := fifelse(
  is.na(padj), "FDR = NA",
  fifelse(padj < 0.001, "FDR < 0.001", paste0("FDR = ", formatC(padj, format = "f", digits = 3)))
)]
maximum.y <- max(c(plot.data$Expression, summary.data$Upper), na.rm = TRUE)
step.y <- max(diff(range(plot.data$Expression, na.rm = TRUE)) * 0.18, 0.25)
pvalue.data[, y.position := maximum.y + step.y * seq_len(.N)]
x.labels <- setNames(
  paste0(summary.data$diabetes_group, "\n(n=", summary.data$n_donor, ")"),
  summary.data$diabetes_group
)

p <- ggplot(plot.data, aes(x = diabetes_group, y = Expression)) +
  geom_violin(aes(fill = diabetes_group), width = 0.82, trim = TRUE, alpha = 0.35, colour = "black", linewidth = 0.45) +
  geom_errorbar(
    data = summary.data,
    aes(x = diabetes_group, ymin = Lower, ymax = Upper),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.9,
    colour = "black"
  ) +
  geom_point(
    data = summary.data,
    aes(x = diabetes_group, y = Mean),
    inherit.aes = FALSE,
    shape = 23,
    size = 4,
    stroke = 0.9,
    fill = "white",
    colour = "black"
  ) +
  geom_jitter(
    aes(fill = diabetes_group),
    shape = 21,
    colour = "black",
    stroke = 0.3,
    width = 0.10,
    size = 2.5,
    alpha = 0.8
  ) +
  ggpubr::stat_pvalue_manual(
    pvalue.data,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.35,
    size = 3.5
  ) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "T1D" = "#E64B35", "T2D" = "#00A087"), drop = FALSE) +
  scale_x_discrete(labels = x.labels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(
    title = "SCN9A expression in Beta cells",
    subtitle = "Diabetes status",
    x = NULL,
    y = "Source-adjusted VST expression"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

ggsave(
  file.path(result_dir, "SCN9A_Beta_DiabetesStatus_VST_Violin.pdf"),
  p,
  width = 6,
  height = 6.5
)

# HbA1c-group model
beta.ids <- rownames(pseudobulk.metadata.all)[
  pseudobulk.metadata.all$Cell_Type == "Beta" & !is.na(pseudobulk.metadata.all$hba1c_group)
]
beta.metadata <- pseudobulk.metadata.all[beta.ids, , drop = FALSE]
beta.counts <- as.matrix(pseudobulk.counts[, beta.ids, drop = FALSE])

dds.beta.hba1c <- DESeqDataSetFromMatrix(
  beta.counts,
  beta.metadata,
  design = ~ source + hba1c_group
)
dds.beta.hba1c <- DESeq(dds.beta.hba1c, quiet = TRUE)

hba1c.comparisons <- list(
  c("Prediabetes", "Normal"),
  c("Diabetes", "Normal"),
  c("Diabetes", "Prediabetes")
)

beta.hba1c.results <- list()
for (comparison in hba1c.comparisons) {
  numerator <- comparison[1]
  denominator <- comparison[2]
  res <- results(
    dds.beta.hba1c,
    contrast = c("hba1c_group", numerator, denominator),
    alpha = 0.05,
    pAdjustMethod = "BH"
  )
  res.dt <- as.data.table(as.data.frame(res), keep.rownames = "Gene")
  res.dt[, Contrast := paste0(numerator, "_vs_", denominator)]
  beta.hba1c.results[[paste(comparison, collapse = "_")]] <- res.dt
}
res.beta.hba1c <- rbindlist(beta.hba1c.results)
fwrite(
  res.beta.hba1c[Gene == "SCN9A"],
  file.path(outdir, "Beta_SCN9A_HbA1c_group_DESeq2.csv")
)

vsd <- varianceStabilizingTransformation(dds.beta.hba1c, blind = FALSE)
vst.matrix <- assay(vsd)
group.design <- model.matrix(~ hba1c_group, data = as.data.frame(colData(dds.beta.hba1c)))
source.adjusted.matrix <- limma::removeBatchEffect(
  vst.matrix,
  batch = colData(dds.beta.hba1c)$source,
  design = group.design
)

plot.data <- data.table(
  PB_ID = rownames(beta.metadata),
  donor_accession = beta.metadata$donor_accession,
  source = beta.metadata$source,
  hba1c_group = beta.metadata$hba1c_group,
  n_cells = beta.metadata$n_cells,
  Expression = as.numeric(source.adjusted.matrix["SCN9A", rownames(beta.metadata)])
)
plot.data[, hba1c_group := factor(hba1c_group, levels = c("Normal", "Prediabetes", "Diabetes"))]

summary.data <- plot.data[, .(
  Mean = mean(Expression),
  SD = sd(Expression),
  n_donor = .N
), by = hba1c_group]
summary.data[, SEM := SD / sqrt(n_donor)]
summary.data[, `:=`(Lower = Mean - SEM, Upper = Mean + SEM)]

pvalue.data <- res.beta.hba1c[
  Gene == "SCN9A",
  .(Contrast, log2FoldChange, pvalue, padj)
]
pvalue.data[, c("group1", "group2") := tstrsplit(Contrast, "_vs_", fixed = TRUE)]
pvalue.data[, label := fifelse(
  is.na(padj), "FDR = NA",
  fifelse(padj < 0.001, "FDR < 0.001", paste0("FDR = ", formatC(padj, format = "f", digits = 3)))
)]
maximum.y <- max(c(plot.data$Expression, summary.data$Upper), na.rm = TRUE)
step.y <- max(diff(range(plot.data$Expression, na.rm = TRUE)) * 0.18, 0.25)
pvalue.data[, y.position := maximum.y + step.y * seq_len(.N)]
x.labels <- setNames(
  paste0(summary.data$hba1c_group, "\n(n=", summary.data$n_donor, ")"),
  summary.data$hba1c_group
)

p <- ggplot(plot.data, aes(x = hba1c_group, y = Expression)) +
  geom_violin(aes(fill = hba1c_group), width = 0.82, trim = TRUE, alpha = 0.35, colour = "black", linewidth = 0.45) +
  geom_errorbar(
    data = summary.data,
    aes(x = hba1c_group, ymin = Lower, ymax = Upper),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.9,
    colour = "black"
  ) +
  geom_point(
    data = summary.data,
    aes(x = hba1c_group, y = Mean),
    inherit.aes = FALSE,
    shape = 23,
    size = 4,
    stroke = 0.9,
    fill = "white",
    colour = "black"
  ) +
  geom_jitter(
    aes(fill = hba1c_group),
    shape = 21,
    colour = "black",
    stroke = 0.3,
    width = 0.10,
    size = 2.5,
    alpha = 0.8
  ) +
  ggpubr::stat_pvalue_manual(
    pvalue.data,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.35,
    size = 3.5
  ) +
  scale_fill_manual(values = c("Normal" = "#4DBBD5", "Prediabetes" = "#EFC000", "Diabetes" = "#E64B35"), drop = FALSE) +
  scale_x_discrete(labels = x.labels, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(
    title = "SCN9A expression in Beta cells",
    subtitle = "HbA1c group",
    x = NULL,
    y = "Source-adjusted VST expression"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

ggsave(
  file.path(result_dir, "SCN9A_Beta_HbA1cGroup_VST_Violin.pdf"),
  p,
  width = 6,
  height = 6.5
)

# -----------------------------------------------------------------------------
# 5. Focused alpha-cell TRPC4 pairwise comparisons
# -----------------------------------------------------------------------------
alpha.ids <- rownames(pseudobulk.metadata.all)[pseudobulk.metadata.all$Cell_Type == "Alpha"]
alpha.metadata <- pseudobulk.metadata.all[alpha.ids, , drop = FALSE]
alpha.counts <- as.matrix(pseudobulk.counts[, alpha.ids, drop = FALSE])

dds.alpha.diabetes <- DESeqDataSetFromMatrix(
  alpha.counts,
  alpha.metadata,
  design = ~ source + diabetes_group
)
dds.alpha.diabetes <- DESeq(dds.alpha.diabetes, quiet = TRUE)

alpha.diabetes.results <- list()
for (comparison in comparisons.all) {
  numerator <- comparison[1]
  denominator <- comparison[2]
  res <- results(
    dds.alpha.diabetes,
    contrast = c("diabetes_group", numerator, denominator),
    alpha = 0.05,
    pAdjustMethod = "BH"
  )
  res.dt <- as.data.table(as.data.frame(res), keep.rownames = "Gene")
  res.dt[, Contrast := paste0(numerator, "_vs_", denominator)]
  alpha.diabetes.results[[paste(comparison, collapse = "_")]] <- res.dt
}
res.alpha.diabetes <- rbindlist(alpha.diabetes.results)
fwrite(
  res.alpha.diabetes[Gene == "TRPC4"],
  file.path(outdir, "Alpha_TRPC4_diabetes_status_DESeq2.csv")
)

alpha.ids <- rownames(pseudobulk.metadata.all)[
  pseudobulk.metadata.all$Cell_Type == "Alpha" & !is.na(pseudobulk.metadata.all$hba1c_group)
]
alpha.metadata <- pseudobulk.metadata.all[alpha.ids, , drop = FALSE]
alpha.counts <- as.matrix(pseudobulk.counts[, alpha.ids, drop = FALSE])

dds.alpha.hba1c <- DESeqDataSetFromMatrix(
  alpha.counts,
  alpha.metadata,
  design = ~ source + hba1c_group
)
dds.alpha.hba1c <- DESeq(dds.alpha.hba1c, quiet = TRUE)

alpha.hba1c.results <- list()
for (comparison in hba1c.comparisons) {
  numerator <- comparison[1]
  denominator <- comparison[2]
  res <- results(
    dds.alpha.hba1c,
    contrast = c("hba1c_group", numerator, denominator),
    alpha = 0.05,
    pAdjustMethod = "BH"
  )
  res.dt <- as.data.table(as.data.frame(res), keep.rownames = "Gene")
  res.dt[, Contrast := paste0(numerator, "_vs_", denominator)]
  alpha.hba1c.results[[paste(comparison, collapse = "_")]] <- res.dt
}
res.alpha.hba1c <- rbindlist(alpha.hba1c.results)
fwrite(
  res.alpha.hba1c[Gene == "TRPC4"],
  file.path(outdir, "Alpha_TRPC4_HbA1c_group_DESeq2.csv")
)
