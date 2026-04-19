# =========================================================
# 这个文件是主分析流程
#
# 1. 读入已经准备好的对象
# 2. 定义 receiver 和 sender
# 3. 找 receiver 里表达的基因和 receptor
# 4. 找候选 ligands
# 5. 在 receiver 内部做差异分析
# 6. 定义 gene set of interest
# 7. 给 ligand 打分并排序
# 8. 取 top ligands
# 9. 往下预测 target genes 和 receptors
# 10. 保存结果，供后面画图使用
# =========================================================

library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(tibble)

# -----------------------------
# 先确保结果目录存在
# -----------------------------
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 读入前一步准备好的对象
# -----------------------------
# readRDS():
# read = 读
# RDS = R 单对象保存格式
# 看到 readRDS，先理解成“把保存好的 R 对象读回来”
seuratObj <- readRDS("data/processed/seuratObj_prepared.rds")
lr_network <- readRDS("data/processed/lr_network_prepared.rds")
ligand_target_matrix <- readRDS("data/processed/ligand_target_matrix_prepared.rds")
weighted_networks <- readRDS("data/processed/weighted_networks_prepared.rds")

# -----------------------------
# 定义 receiver 和 sender
# -----------------------------
# receiver:
# 被信号影响、重点观察的细胞类型
#
# sender_celltypes:
# 可能发出信号的细胞类型
receiver <- "CD8 T"
sender_celltypes <- c("CD4 T", "Treg", "Mono", "NK", "B", "DC")

cat("Receiver: ", receiver, "\n", sep = "")
cat("Sender cell types: ", paste(sender_celltypes, collapse = ", "), "\n", sep = "")

# -----------------------------
# 找 receiver 中表达的基因
# -----------------------------
# get_expressed_genes():
# get = 获取
# expressed genes = 表达基因
#
# pct = 0.05:
# 至少 5% 的这类细胞表达该基因，才算“表达”
expressed_genes_receiver <- get_expressed_genes(
  receiver,
  seuratObj,
  pct = 0.05
)

cat("Number of expressed genes in receiver: ", length(expressed_genes_receiver), "\n", sep = "")
print(head(expressed_genes_receiver))

# -----------------------------
# 找 receiver 中表达的 receptor
# -----------------------------
# unique():
# 去重
# 看到 unique，先记成“只保留唯一值”
all_receptors <- unique(lr_network$to)

# intersect():
# 交集
# inter = 交叉
# sect = 部分
# 看到 intersect，直接理解成“两边共同拥有的部分”
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

cat("Number of expressed receptors in receiver: ", length(expressed_receptors), "\n", sep = "")
print(head(expressed_receptors))

# -----------------------------
# 找 sender-agnostic 候选 ligands
# -----------------------------
# 这里先不管 sender 是否真的表达
# 只要 receiver 有对应 receptor，就先把 ligand 放进候选列表
potential_ligands <- lr_network %>%
  # filter():
  # 筛选行
  # filter 这个词本身就是“过滤”
  filter(to %in% expressed_receptors) %>%

  # pull():
  # 从数据框里“抽出”某一列
  pull(from) %>%

  # unique():
  # 去重
  unique()

cat("Number of sender-agnostic ligands: ", length(potential_ligands), "\n", sep = "")
print(head(potential_ligands))

# -----------------------------
# 找 sender-focused 候选 ligands
# -----------------------------
# 对每一种 sender 细胞，分别找表达基因
#
# lapply():
# l = list
# apply = 应用
# 先记成“把同一个函数作用到列表里的每个元素上”
list_expressed_genes_sender <- sender_celltypes %>%
  unique() %>%
  lapply(get_expressed_genes, seuratObj, 0.05)

# unlist():
# 把列表拆开，变成一个长向量
expressed_genes_sender <- list_expressed_genes_sender %>%
  unlist() %>%
  unique()

# focused 版本更严格：
# 不仅 receiver 能接到，还要求 sender 里真的表达这个 ligand
potential_ligands_focused <- intersect(
  potential_ligands,
  expressed_genes_sender
)

cat("Number of sender-focused ligands: ", length(potential_ligands_focused), "\n", sep = "")
print(head(potential_ligands_focused))

# -----------------------------
# 定义实验条件
# -----------------------------
condition_oi <- "LCMV"
condition_reference <- "SS"

cat("Condition of interest: ", condition_oi, "\n", sep = "")
cat("Reference condition: ", condition_reference, "\n", sep = "")

# -----------------------------
# 先从整个对象里取出 receiver 细胞
# -----------------------------
# subset():
# sub = 子 / 下一级
# set = 集合
# 看到 subset，先理解成“从大对象里切一个子集出来”
seurat_obj_receiver <- subset(seuratObj, idents = receiver)

cat("Receiver subset object:\n")
print(seurat_obj_receiver)

# -----------------------------
# 在 receiver 内部做差异表达分析
# -----------------------------
# FindMarkers():
# marker 可以理解成“区分两组的特征”
# 这里用它找出 LCMV 和 SS 条件下，CD8 T 内部差异表达的基因
DE_table_receiver <- FindMarkers(
  object = seurat_obj_receiver,
  ident.1 = condition_oi,
  ident.2 = condition_reference,
  group.by = "aggregate",
  min.pct = 0.05
) %>%
  # rownames_to_column():
  # 把原来藏在行名里的基因名变成正式的一列
  rownames_to_column("gene")

cat("DE result head:\n")
print(head(DE_table_receiver))

cat("DE result column names:\n")
print(colnames(DE_table_receiver))

# -----------------------------
# 定义 gene set of interest
# -----------------------------
# 这里沿用差异分析结果里的 avg_log2FC
fc_column <- "avg_log2FC"

# 先筛显著、再筛变化幅度
#
# p_val_adj <= 0.05:
# 校正后的 p 值显著
#
# abs(... ) >= 0.25:
# abs = absolute value（绝对值）
# 这里表示无论上调还是下调，只看变化幅度是否够大
geneset_oi <- DE_table_receiver %>%
  filter(
    p_val_adj <= 0.05,
    abs(.data[[fc_column]]) >= 0.25
  ) %>%
  pull(gene)

# 只保留 ligand_target_matrix 里存在的基因
# 否则后面 NicheNet 无法对这些基因打分
geneset_oi <- geneset_oi[geneset_oi %in% rownames(ligand_target_matrix)]

cat("Number of genes in gene set of interest: ", length(geneset_oi), "\n", sep = "")
print(head(geneset_oi))

# -----------------------------
# 定义 background gene set
# -----------------------------
# 背景基因：
# receiver 中表达，并且也出现在 ligand_target_matrix 里的基因
background_expressed_genes <- expressed_genes_receiver[
  expressed_genes_receiver %in% rownames(ligand_target_matrix)
]

cat("Number of background genes: ", length(background_expressed_genes), "\n", sep = "")
print(head(background_expressed_genes))

# -----------------------------
# 计算 ligand activity
# -----------------------------
# predict_ligand_activities():
# predict = 预测
# ligand activities = ligand 活性
#
# 这里的核心问题可以直接记成：
# “哪个 ligand 最能解释 geneset_oi 这组基因变化？”
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# arrange():
# 排序
#
# desc():
# descending = 降序
#
# mutate():
# 增加或修改列
# 先记成“给表补一列/改一列”
ligand_activities <- ligand_activities %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = rank(desc(aupr_corrected)))

cat("Top ligand activities:\n")
print(head(ligand_activities, 10))

# -----------------------------
# 取前 30 个 top ligands
# -----------------------------
# slice_max():
# 取某一列数值最大的前 n 个
best_upstream_ligands <- ligand_activities %>%
  slice_max(order_by = aupr_corrected, n = 30) %>%
  arrange(desc(aupr_corrected)) %>%
  pull(test_ligand)

cat("Number of top ligands: ", length(best_upstream_ligands), "\n", sep = "")
print(head(best_upstream_ligands, 10))

# -----------------------------
# 顺手检查几个关键结果
# -----------------------------
cat("geneset_oi length: ", length(geneset_oi), "\n", sep = "")
cat("background_expressed_genes length: ", length(background_expressed_genes), "\n", sep = "")
cat("Top 10 ligand activities:\n")
print(head(ligand_activities, 10))
cat("Top 10 best ligands:\n")
print(head(best_upstream_ligands, 10))

# -----------------------------
# 预测 top ligands 的 target genes
# -----------------------------
# 这里的逻辑是：
# 对每个 top ligand，分别往下找最可能影响的 target genes
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 100
  ) %>%

  # bind_rows():
  # bind = 绑定
  # rows = 行
  # 把多个结果表按行拼起来
  bind_rows() %>%

  # drop_na():
  # drop = 去掉
  # na = 缺失值
  drop_na()

cat("Active ligand-target links head:\n")
print(head(active_ligand_target_links_df))
cat("Number of active ligand-target links rows: ", nrow(active_ligand_target_links_df), "\n", sep = "")

# -----------------------------
# 预测 top ligands 的 receptors
# -----------------------------
# 这里结合：
# 1. top ligands
# 2. receiver 中表达的 receptors
# 3. ligand-receptor 网络
# 4. 带权重的 lr_sig 网络
# 来评估更优先的 ligand-receptor 配对
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
# 把表格结果存下来
# -----------------------------
# write.csv():
# write = 写出
# csv = 逗号分隔表格格式
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
# 把关键对象打包存成一个 RDS
# -----------------------------
# list():
# 把多个对象打包成一个列表
#
# 这样做的好处：
# 后面的画图脚本可以直接读取 analysis_results.rds
# 不用每次都把主分析重跑一遍
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

cat("02_run_nichenet_analysis.R 跑完：主分析结果已经保存好。\n")
