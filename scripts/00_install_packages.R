# =========================================================
# 文件名: 00_install_packages.R
# 作用:
#   这个脚本专门负责安装并加载本项目所需的 R 包。
#
# 为什么要单独写成一个文件:
#   因为“安装环境”和“正式分析”最好分开。
#   这样主分析脚本会更干净，别人看也更清楚。
#
# 本脚本做的事情:
#   1. 检查 devtools 是否安装
#   2. 安装常用 CRAN 包：Seurat / SeuratObject / tidyverse
#   3. 从 GitHub 安装 nichenetr
#   4. 加载所有分析会用到的包
#
# 使用方法:
#   source("scripts/00_install_packages.R")
# 或者在 RStudio 里直接运行整个脚本
# =========================================================

# -----------------------------
# 第 0 步：定义一个“小工具函数”
# -----------------------------
# 这个函数的作用是：
# 如果某个包还没有安装，就自动帮你安装。
#
# 这样写的好处：
# 1. 代码更简洁
# 2. 不用每个包都手动写 if
# 3. 初学时更容易理解“先检查，再安装”的思路
install_if_missing <- function(pkg) {
  # requireNamespace() 的作用：
  # 检查某个包是否已经存在于当前 R 环境中
  #
  # quietly = TRUE 的意思是：
  # 检查时尽量少打印无关提示
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# -----------------------------
# 第 1 步：安装 devtools
# -----------------------------
# devtools 主要用来从 GitHub 安装包
# 因为 nichenetr 常用 GitHub 安装，所以这里先确保 devtools 存在
install_if_missing("devtools")

# -----------------------------
# 第 2 步：安装 CRAN 基础包
# -----------------------------
# 这里把项目需要的 CRAN 包先列成一个向量
cran_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "tidyverse",
  "tibble"
)

# 循环安装这些包
# for (...) 是最基础的循环结构之一
# 这里的意思是：对 cran_pkgs 里的每个包名，都执行一次安装检查
for (pkg in cran_pkgs) {
  install_if_missing(pkg)
}

# -----------------------------
# 第 3 步：安装 nichenetr
# -----------------------------
# nichenetr 不一定总在 CRAN 上，因此这里用 GitHub 安装
#
# 说明：
# 1. 只有当 nichenetr 没有安装时才安装
# 2. 如果你的本地已经装好了，就不会重复安装
if (!requireNamespace("nichenetr", quietly = TRUE)) {
  devtools::install_github("saeyslab/nichenetr")
}

# -----------------------------
# 第 4 步：加载项目需要的包
# -----------------------------
# library() 的作用是：
# 把已经安装好的包加载到当前 R 会话里
#
# 为什么要显式写出来：
# 别人一看就知道这个项目依赖哪些包
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(tibble)

# -----------------------------
# 第 5 步：打印版本信息
# -----------------------------
# packageVersion() 的作用是查看包版本
# 这一步不是必须的，但很适合学习型仓库：
# 你自己和师姐都能看到你用的是哪个版本
cat("========== Package versions ==========\n")
cat("Seurat: ", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("SeuratObject: ", as.character(packageVersion("SeuratObject")), "\n", sep = "")
cat("nichenetr: ", as.character(packageVersion("nichenetr")), "\n", sep = "")
cat("tidyverse: ", as.character(packageVersion("tidyverse")), "\n", sep = "")
cat("=====================================\n")

# -----------------------------
# 第 6 步：打印完成提示
# -----------------------------
cat("00_install_packages.R 运行完成：依赖包已准备好。\n")
