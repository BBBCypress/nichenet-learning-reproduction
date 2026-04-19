# =========================================================
# 这个文件只负责“数据准备”
#
# 先做这些事：
# 1. 建目录
# 2. 下载示例数据
# 3. 读取 Seurat 对象
# 4. 读取 NicheNet 网络文件
# 5. 做基础检查
# 6. 把处理好的对象存下来
#
# 这一部分先跑通，后面主分析才有东西可用
# 如果这里出问题，先查下载、读入和对象格式，不要急着改分析参数
# =========================================================

library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(tibble)

# -----------------------------
# 先把项目目录建出来
# -----------------------------
# dir.create()：
# dir = directory（目录/文件夹）
# create = 创建
#
# recursive = TRUE：
# 上层目录不存在时，一起建出来
#
# showWarnings = FALSE：
# 如果目录已经存在，就别重复弹警告
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 把下载超时时间设长一点
# -----------------------------
# options()：
# 用来设置 R 的全局选项
#
# timeout = 超时时间
# 下载站点慢的时候，这一步能减少中途断掉的情况
options(timeout = 600)

# -----------------------------
# 先把下载地址和本地文件名写清楚
# -----------------------------
# 这样后面读代码时更容易看出：
# 文件从哪里来，最后存到哪里
seurat_url <- "https://zenodo.org/record/3531889/files/seuratObj.rds"
seurat_file <- "data/raw/seuratObj.rds"

lr_url <- "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"
lr_file <- "data/raw/lr_network_mouse_21122021.rds"

ltm_url <- "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
ltm_file <- "data/raw/ligand_target_matrix_mouse.rds"

wn_url <- "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"
wn_file <- "data/raw/weighted_networks_mouse.rds"

# -----------------------------
# 没有下载过的文件，再下载
# -----------------------------
# file.exists()：
# exists = 存在
# 一看到这个函数，就先理解成“文件在不在”
#
# download.file()：
# download = 下载
# file = 文件
#
# mode = "wb"：
# 二进制文件用 wb 更稳，.rds 这里就按这个写
if (!file.exists(seurat_file)) {
  download.file(
    url = seurat_url,
    destfile = seurat_file,
    mode = "wb"
  )
}

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
# 简单确认一下文件有没有真的下到本地
# -----------------------------
# cat()：
# 常用来把文字直接打印到控制台
# 先记成“直接输出一行提示”
cat("Seurat object exists: ", file.exists(seurat_file), "\n", sep = "")
cat("LR network exists: ", file.exists(lr_file), "\n", sep = "")
cat("Ligand-target matrix exists: ", file.exists(ltm_file), "\n", sep = "")
cat("Weighted networks exists: ", file.exists(wn_file), "\n", sep = "")

# -----------------------------
# 读取 Seurat 对象
# -----------------------------
# readRDS()：
# read = 读取
# RDS = R 的单对象保存格式
#
# 看到 readRDS，直接记成：
# “把 .rds 文件里的 R 对象读回来”
seuratObj <- readRDS(seurat_file)

# -----------------------------
# 更新旧版 Seurat 对象格式
# -----------------------------
# UpdateSeuratObject()：
# 名字很直白，看到就知道是“更新 Seurat 对象”
#
# 官方示例对象有时来自旧版本 Seurat
# 这里先转成当前版本更容易继续往下用的结构
seuratObj <- UpdateSeuratObject(seuratObj)

# -----------------------------
# 把基因别名换成标准 gene symbol
# -----------------------------
# alias = 别名
# symbol = 标准名
#
# 这一步主要是为了后面和 NicheNet 网络对基因名时更稳
# 不然同一个基因可能因为命名不统一匹配不上
#
# "mouse" 表示这里按小鼠数据处理
seuratObj <- alias_to_symbol_seurat(seuratObj, "mouse")

# -----------------------------
# 把当前 identity 设成 celltype
# -----------------------------
# Idents()：
# Seurat 里很常见
# 先简单记成“当前按什么标签给细胞分组”
#
# 这里设成 celltype，后面按细胞类型取子集会更直接
Idents(seuratObj) <- seuratObj$celltype

# -----------------------------
# 读取 NicheNet 网络对象
# -----------------------------
lr_network <- readRDS(lr_file)
ligand_target_matrix <- readRDS(ltm_file)
weighted_networks <- readRDS(wn_file)

# -----------------------------
# 去掉重复的 ligand-receptor 配对
# -----------------------------
# distinct()：
# 来自 dplyr
# 先记成“去重”
#
# 这里按 from 和 to 两列一起去重
# 只保留唯一的 ligand-receptor 组合
lr_network <- lr_network %>% distinct(from, to)

# -----------------------------
# 做几个最基础的检查
# -----------------------------
# 先检查对象本身、meta.data、细胞类型分布、条件分布
# 再检查网络对象的基本结构
# 这一部分虽然简单，但很有用：
# 先确认输入没问题，再往后做分析
cat("\n========== Basic checks ==========\n")

cat("Seurat object summary:\n")
print(seuratObj)

cat("\nMeta data head:\n")
# head()：
# head = 头部
# 看前几行，快速确认结构
print(head(seuratObj@meta.data))

cat("\nCelltype distribution:\n")
# table()：
# 把每个类别出现多少次统计出来
print(table(seuratObj$celltype))

cat("\nCondition distribution (aggregate):\n")
print(table(seuratObj$aggregate))

cat("\nLR network head:\n")
print(head(lr_network))

cat("\nLigand-target matrix dimension:\n")
# dim()：
# dimension = 维度
# 经常用来看矩阵或数据框的行列数
print(dim(ligand_target_matrix))

cat("\nWeighted networks names:\n")
# names()：
# 看列表或对象里有哪些元素名
print(names(weighted_networks))

cat("=================================\n")

# -----------------------------
# 把处理好的对象存下来
# -----------------------------
# saveRDS() 和 readRDS() 是一对
# save = 保存
#
# 这样做的好处：
# 后面的主分析脚本可以直接读取处理好的对象
# 不用每次都从头下载、更新、转换
saveRDS(seuratObj, "data/processed/seuratObj_prepared.rds")
saveRDS(lr_network, "data/processed/lr_network_prepared.rds")
saveRDS(ligand_target_matrix, "data/processed/ligand_target_matrix_prepared.rds")
saveRDS(weighted_networks, "data/processed/weighted_networks_prepared.rds")

cat("01_prepare_data.R 跑完：数据和网络文件已经准备好。\n")
