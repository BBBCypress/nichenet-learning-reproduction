# =========================================================
# 文件名: 01_prepare_data.R
# 作用:
#   这个脚本负责准备分析所需的数据输入。
#
# 本脚本做的事情:
#   1. 创建项目目录
#   2. 下载官方示例 Seurat 对象
#   3. 下载 NicheNet 网络文件
#   4. 读取并预处理 Seurat 对象
#   5. 检查对象和网络文件是否正常
#   6. 把处理好的对象保存到 data/processed/
#
# 为什么要单独拆出来:
#   因为“数据准备”和“正式分析”是两件事。
#   单独拆开后，主分析脚本会更干净，也更容易重跑。
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
# 第 1 步：创建项目目录
# -----------------------------
# dir.create() 的作用：
# 创建文件夹
#
# recursive = TRUE 表示：
# 如果上层目录还不存在，也一起创建
#
# showWarnings = FALSE 表示：
# 如果目录已经存在，不要重复弹警告
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 第 2 步：设置下载超时时间
# -----------------------------
# options(timeout = 600) 的意思是：
# 把下载等待时间设成 600 秒（10 分钟）
# 避免 Zenodo 速度慢时中途断掉
options(timeout = 600)

# -----------------------------
# 第 3 步：定义要下载的文件路径
# -----------------------------
# 这样做的好处是：
# 1. 以后如果要改文件名或链接，只改这里
# 2. 后面写代码时更简洁
seurat_url <- "https://zenodo.org/record/3531889/files/seuratObj.rds"
seurat_file <- "data/raw/seuratObj.rds"

lr_url <- "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"
lr_file <- "data/raw/lr_network_mouse_21122021.rds"

ltm_url <- "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
ltm_file <- "data/raw/ligand_target_matrix_mouse.rds"

wn_url <- "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"
wn_file <- "data/raw/weighted_networks_mouse.rds"

# -----------------------------
# 第 4 步：下载 Seurat 对象
# -----------------------------
# file.exists() 的作用：
# 检查某个文件是否已经存在
#
# 这里的逻辑是：
# 如果文件还没有下载，就下载；
# 如果已经有了，就直接跳过，避免重复下载。
if (!file.exists(seurat_file)) {
  download.file(
    url = seurat_url,
    destfile = seurat_file,
    mode = "wb"
  )
}

# -----------------------------
# 第 5 步：下载 NicheNet 网络文件
# -----------------------------
# 这三个文件分别代表：
# 1. ligand-receptor 配对网络
# 2. ligand-target 先验矩阵
# 3. 带权重的网络信息
if (!file.exists(lr_file)) {
  download.file(
    url = lr_url,
    destfile = lr_file,
    mode = "wb"
  )
}

if (!file.exists(ltm_file)) {
  download.file(
    url = ltm_url,
    destfile = ltm_file,
    mode = "wb"
  )
}

if (!file.exists(wn_file)) {
  download.file(
    url = wn_url,
    destfile = wn_file,
    mode = "wb"
  )
}

# -----------------------------
# 第 6 步：检查原始文件是否下载成功
# -----------------------------
cat("Seurat object exists: ", file.exists(seurat_file), "\n", sep = "")
cat("LR network exists: ", file.exists(lr_file), "\n", sep = "")
cat("Ligand-target matrix exists: ", file.exists(ltm_file), "\n", sep = "")
cat("Weighted networks exists: ", file.exists(wn_file), "\n", sep = "")

# -----------------------------
# 第 7 步：读取 Seurat 对象
# -----------------------------
# readRDS() 用来读取 .rds 文件
seuratObj <- readRDS(seurat_file)

# -----------------------------
# 第 8 步：更新旧版 Seurat 对象结构
# -----------------------------
# 有些官方示例对象可能来自较老版本 Seurat
# UpdateSeuratObject() 会把它升级成当前版本可用的结构
seuratObj <- UpdateSeuratObject(seuratObj)

# -----------------------------
# 第 9 步：把基因别名转成标准 symbol
# -----------------------------
# alias_to_symbol_seurat() 是 nichenetr 提供的函数
# 这样做能减少后面和 NicheNet 网络匹配不上基因名的问题
#
# "mouse" 表示当前数据是小鼠数据
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

# -----------------------------
# 第 10 步：设置当前 identity
# -----------------------------
# Idents() 是 Seurat 里很重要的概念
# 这里把当前细胞分组方式设为 celltype
Idents(seuratObj) <- seuratObj$celltype

# -----------------------------
# 第 11 步：读取 NicheNet 网络对象
# -----------------------------
lr_network <- readRDS(lr_file)
ligand_target_matrix <- readRDS(ltm_file)
weighted_networks <- readRDS(wn_file)

# 去掉重复的 ligand-receptor 配对
# distinct(from, to) 表示：
# 只保留唯一的 ligand-receptor 组合
lr_network <- lr_network %>% distinct(from, to)

# -----------------------------
# 第 12 步：做基础数据检查
# -----------------------------
# 这些检查很适合学习性项目：
# 既能帮助你确认对象没读错，也能体现你有“先检查再分析”的习惯
cat("\n========== Basic checks ==========\n")
cat("Seurat object summary:\n")
print(seuratObj)

cat("\nMeta data head:\n")
print(head(seuratObj@meta.data))

cat("\nCelltype distribution:\n")
print(table(seuratObj$celltype))

cat("\nCondition distribution (aggregate):\n")
print(table(seuratObj$aggregate))

cat("\nLR network head:\n")
print(head(lr_network))

cat("\nLigand-target matrix dimension:\n")
print(dim(ligand_target_matrix))

cat("\nWeighted networks names:\n")
print(names(weighted_networks))
cat("=================================\n")

# -----------------------------
# 第 13 步：保存处理好的对象
# -----------------------------
# saveRDS() 的作用：
# 把对象保存成 .rds 文件，方便后续脚本直接读取
saveRDS(seuratObj, "data/processed/seuratObj_prepared.rds")
saveRDS(lr_network, "data/processed/lr_network_prepared.rds")
saveRDS(ligand_target_matrix, "data/processed/ligand_target_matrix_prepared.rds")
saveRDS(weighted_networks, "data/processed/weighted_networks_prepared.rds")

# -----------------------------
# 第 14 步：打印完成提示
# -----------------------------
cat("01_prepare_data.R 运行完成：数据和网络文件已准备好。\n")
