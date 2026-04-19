# =========================================================
# 文件名: 02_run_nichenet_analysis.R
# 作用:
#   这是项目的主分析脚本。
#
# 本脚本做的事情:
#   1. 读取处理好的 Seurat 对象和网络文件
#   2. 定义 receiver / sender
#   3. 提取 receiver 表达基因和 receptor
#   4. 找候选 ligands
#   5. 在 receiver 内部做差异表达分析
#   6. 定义 gene set of interest
#   7. 计算 ligand activity
#   8. 取 top ligands
#   9. 预测 top ligands 的 target genes 和 receptors
#   10. 保存关键结果
#
# 这是你整个学习型复现仓库中最核心的分析脚本。
# =========================================================

# -----------------------------
# 第 0 步：加载包
# -----------------------------
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(tibble)

# -----------------------------
# 第 1 步：创建结果目录
# -----------------------------
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 第 2 步：读取处理好的输入对象
# -----------------------------
seuratObj <- readRDS("data/processed/seuratObj_prepared.rds")
lr_network <- readRDS("data/processed/lr_network_prepared.rds")
ligand_target_matrix <- readRDS("data/processed/ligand_target_matrix_prepared.rds")
weighted_networks <- readRDS("data/processed/weighted_networks_prepared.rds")

# -----------------------------
# 第 3 步：定义 receiver 和 sender
# -----------------------------
# receiver = 被信号影响的细胞类型
receiver <- "CD8 T"

# sender_celltypes = 可能发出信号的细胞类型
sender_celltypes <- c("CD4 T", "Treg", "Mono", "NK", "B", "DC")

cat("Receiver: ", receiver, "\n", sep = "")
cat("Sender cell types: ", paste(sender_celltypes, collapse = ", "), "\n", sep = "")

# -----------------------------
# 第 4 步：找 receiver 中表达的基因
# -----------------------------
# get_expressed_genes() 的作用：
# 找某类细胞中“表达出来”的基因
#
# pct = 0.05 的意思：
# 至少 5% 的该类细胞表达这个基因，才认为“表达”
expressed_genes_receiver <- get_expressed_genes(
  receiver,
  seuratObj,
  pct = 0.05
)

cat("Number of expressed genes in receiver: ", length(expressed_genes_receiver), "\n", sep = "")
print(head(expressed_genes_receiver))

# -----------------------------
# 第 5 步：找 receiver 中表达的 receptor
# -----------------------------
# lr_network$to 表示配对网络中的 receptor 一列
all_receptors <- unique(lr_network$to)

# intersect() 的作用：
# 求交集
# 这里的意思是：哪些 receptor 同时也在 receiver 中表达
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

cat("Number of expressed receptors in receiver: ", length(expressed_receptors), "\n", sep = "")
print(head(expressed_receptors))

# -----------------------------
# 第 6 步：得到 sender-agnostic 候选 ligand
# -----------------------------
# 逻辑：
# 只要某个 ligand 对应的 receptor 在 receiver 中表达，
# 它就先进入候选列表
potential_ligands <- lr_network %>%
  filter(to %in% expressed_receptors) %>%
  pull(from) %>%
  unique()

cat("Number of sender-agnostic ligands: ", length(potential_ligands), "\n", sep = "")
print(head(potential_ligands))

# -----------------------------
# 第 7 步：得到 sender-focused 候选 ligand
# -----------------------------
# 对每一种 sender 细胞，分别找表达基因
list_expressed_genes_sender <- sender_celltypes %>%
  unique() %>%
  lapply(get_expressed_genes, seuratObj, 0.05)

# 把多个 sender 的表达基因合并成一个长向量
expressed_genes_sender <- list_expressed_genes_sender %>%
  unlist() %>%
  unique()

# sender-focused 候选 ligand =
# receiver 能接收 + sender 真的表达
potential_ligands_focused <- intersect(
  potential_ligands,
  expressed_genes_sender
)

cat("Number of sender-focused ligands: ", length(potential_ligands_focused), "\n", sep = "")
print(head(potential_ligands_focused))

# -----------------------------
# 第 8 步：定义分析条件
# -----------------------------
condition_oi <- "LCMV"
condition_reference <- "SS"

cat("Condition of interest: ", condition_oi, "\n", sep = "")
cat("Reference condition: ", condition_reference, "\n", sep = "")

# -----------------------------
# 第 9 步：取出 receiver 细胞
# -----------------------------
# subset() 的作用：
# 从 Seurat 对象中取子集
seurat_obj_receiver <- subset(seuratObj, idents = receiver)

cat("Receiver subset object:\n")
print(seurat_obj_receiver)

# -----------------------------
# 第 10 步：在 receiver 内部做差异表达分析
# -----------------------------
# FindMarkers() 是 Seurat 中做差异分析的经典函数
#
# group.by = "aggregate" 表示：
# 按 meta.data 里的 aggregate 列进行条件分组
DE_table_receiver <- FindMarkers(
  object = seurat_obj_receiver,
  ident.1 = condition_oi,
  ident.2 = condition_reference,
  group.by = "aggregate",
  min.pct = 0.05
) %>%
  rownames_to_column("gene")

cat("DE result head:\n")
print(head(DE_table_receiver))

cat("DE result column names:\n")
print(colnames(DE_table_receiver))

# -----------------------------
# 第 11 步：定义 gene set of interest
# -----------------------------
# 这里我们按你的原始思路，使用 avg_log2FC 这一列
fc_column <- "avg_log2FC"

# 筛选标准：
# 1. p_val_adj <= 0.05
# 2. abs(log2FC) >= 0.25
geneset_oi <- DE_table_receiver %>%
  filter(
    p_val_adj <= 0.05,
    abs(.data[[fc_column]]) >= 0.25
  ) %>%
  pull(gene)

# 只保留 ligand_target_matrix 中存在的基因
geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]

cat("Number of genes in gene set of interest: ", length(geneset_oi), "\n", sep = "")
print(head(geneset_oi))

# -----------------------------
# 第 12 步：定义 background gene set
# -----------------------------
# 背景基因 = receiver 中表达，且存在于 ligand_target_matrix 中的基因
background_expressed_genes <- expressed_genes_receiver[
  expressed_genes_receiver %in% rownames(ligand_target_matrix)
]

cat("Number of background genes: ", length(background_expressed_genes), "\n", sep = "")
print(head(background_expressed_genes))

# -----------------------------
# 第 13 步：做 ligand activity analysis
# -----------------------------
# predict_ligand_activities() 的作用：
# 给每个候选 ligand 打分，看看谁最能解释 gene set of interest
#
# 注意：
# 这里保留你原始脚本的思路，使用 potential_ligands
# 而不是 potential_ligands_focused
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# 按 aupr_corrected 从高到低排序，并给出排名
ligand_activities <- ligand_activities %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = rank(desc(aupr_corrected)))

cat("Top ligand activities:\n")
print(head(ligand_activities, 10))

# -----------------------------
# 第 14 步：提取 top ligands
# -----------------------------
# 这里取前 30 个高分 ligand
#
# slice_max() 比 top_n() 更新，也更推荐
best_upstream_ligands <- ligand_activities %>%
  slice_max(order_by = aupr_corrected, n = 30) %>%
  arrange(desc(aupr_corrected)) %>%
  pull(test_ligand)

cat("Number of top ligands: ", length(best_upstream_ligands), "\n", sep = "")
print(head(best_upstream_ligands, 10))

# -----------------------------
# 第 15 步：检查关键结果
# -----------------------------
cat("geneset_oi length: ", length(geneset_oi), "\n", sep = "")
cat("background_expressed_genes length: ", length(background_expressed_genes), "\n", sep = "")
cat("Top 10 ligand activities:\n")
print(head(ligand_activities, 10))
cat("Top 10 best ligands:\n")
print(head(best_upstream_ligands, 10))

# -----------------------------
# 第 16 步：预测 top ligands 的 target genes
# -----------------------------
# 这里对每个 top ligand 找它最可能影响的 target genes
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 100
  ) %>%
  bind_rows() %>%
  drop_na()

cat("Active ligand-target links head:\n")
print(head(active_ligand_target_links_df))
cat("Number of active ligand-target links rows: ", nrow(active_ligand_target_links_df), "\n", sep = "")

# -----------------------------
# 第 17 步：预测 top ligands 的 receptors
# -----------------------------
# 根据：
# 1. top ligands
# 2. receiver 中表达的 receptors
# 3. ligand-receptor 网络
# 4. weighted_networks$lr_sig
# 来推测优先的 ligand-receptor 配对
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands,
  expressed_receptors,
  lr_network,
  weighted_networks$lr_sig
)

cat("Ligand-receptor links head:\n")
print(head(ligand_receptor_links_df))
cat("Number of ligand-receptor links rows: ", nrow(ligand_receptor_links_df), "\n", sep = "")

# -----------------------------
# 第 18 步：保存表格结果
# -----------------------------
write.csv(
  DE_table_receiver,
  "results/tables/DE_table_receiver.csv",
  row.names = FALSE
)

write.csv(
  ligand_activities,
  "results/tables/ligand_activities.csv",
  row.names = FALSE
)

write.csv(
  active_ligand_target_links_df,
  "results/tables/active_ligand_target_links.csv",
  row.names = FALSE
)

write.csv(
  ligand_receptor_links_df,
  "results/tables/ligand_receptor_links.csv",
  row.names = FALSE
)

write.csv(
  data.frame(best_upstream_ligands = best_upstream_ligands),
  "results/tables/best_upstream_ligands.csv",
  row.names = FALSE
)

# -----------------------------
# 第 19 步：把关键对象保存成一个总的 RDS
# -----------------------------
# 这样后面的画图脚本就可以直接读取，不用重新跑分析
analysis_results <- list(
  receiver = receiver,
  sender_celltypes = sender_celltypes,
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  expressed_genes_receiver = expressed_genes_receiver,
  expressed_receptors = expressed_receptors,
  potential_ligands = potential_ligands,
  potential_ligands_focused = potential_ligands_focused,
  DE_table_receiver = DE_table_receiver,
  geneset_oi = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_activities = ligand_activities,
  best_upstream_ligands = best_upstream_ligands,
  active_ligand_target_links_df = active_ligand_target_links_df,
  ligand_receptor_links_df = ligand_receptor_links_df,
  ligand_target_matrix = ligand_target_matrix
)

saveRDS(analysis_results, "results/rds/analysis_results.rds")

# -----------------------------
# 第 20 步：打印完成提示
# -----------------------------
cat("02_run_nichenet_analysis.R 运行完成：主分析结果已保存。\n")
