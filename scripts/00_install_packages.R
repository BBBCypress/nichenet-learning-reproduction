# =========================================================
# 先把环境准备好
# 先安装需要的包，再加载需要的包
# 把“环境准备”和“正式分析”分开写，后面排错会更清楚
# 如果这里出问题，先看包有没有装好，不要急着怀疑分析流程
# 定义一个小工具函数：缺哪个包，就装哪个包
# function = 定义一个函数
# 可以先把 function 理解成“自定义一个可重复使用的小工具”
# 这里把函数名写成 install_if_missing
# install = 安装
# missing = 缺少
# 这个名字很好记：缺什么，就安装什么
install_if_missing <- function(pkg) {

  # requireNamespace()：
  # 先检查某个包在不在
  #
  # require = 需要
  # namespace = 命名空间
  #
  # 现在不用死抠“命名空间”的严格定义，
  # 先记住这个函数在这里的实际用途：
  # “检查这个包是否已经可用”
  #
  # quietly = TRUE：
  # quietly = 安静地
  # 表示检查时尽量不要打印太多提示
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# -----------------------------
# 先确保 devtools 已经安装
# -----------------------------
# devtools 主要用来从 GitHub 安装包
# 后面安装 nichenetr 会用到它
install_if_missing("devtools")

# -----------------------------
# 列出这次分析需要的 CRAN 包
# -----------------------------
# c() 的作用是把多个值合成一个向量
# 可以先把 c 记成 collect / combine 的感觉：
# “把几个东西收成一串”
cran_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "tidyverse",
  "tibble"
)

# -----------------------------
# 用 for 循环逐个检查、逐个安装
# -----------------------------
# for = 对这一组对象一个一个地做同样的事
# 这里就是：对 cran_pkgs 里的每个包名，都跑一次 install_if_missing()
for (pkg in cran_pkgs) {
  install_if_missing(pkg)
}

# -----------------------------
# 安装 nichenetr
# -----------------------------
# 这个包通常直接从 GitHub 安装
# 所以这里先检查：如果没装，再安装
if (!requireNamespace("nichenetr", quietly = TRUE)) {

  # devtools::install_github()
  #
  # :: 的意思：
  # 明确调用“某个包里的某个函数”
  # 这里表示：调用 devtools 包里的 install_github()
  #
  # install_github 这个名字很好记：
  # install = 安装
  # github = 从 GitHub 来
  devtools::install_github("saeyslab/nichenetr")
}

# -----------------------------
# 正式加载后面分析会用到的包
# -----------------------------
# library()：
# 把包加载到当前 R 会话里
#
# library 原义是“图书馆”
# 这里可以记成：
# “把这本工具书拿到手边开始用”
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(tibble)

# -----------------------------
# 顺手打印版本信息，方便以后查环境
# -----------------------------
# packageVersion()：
# package = 包
# version = 版本
# 一看到这个函数，就知道它是“查看包版本”
#
# as.character()：
# 把版本对象转成字符，打印时更稳妥
#
# cat()：
# cat 在 R 里常用来把文字直接输出到控制台
# 可以先记成“把内容连续打印出来”
cat("========== Package versions ==========\n")
cat("Seurat: ", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("SeuratObject: ", as.character(packageVersion("SeuratObject")), "\n", sep = "")
cat("nichenetr: ", as.character(packageVersion("nichenetr")), "\n", sep = "")
cat("tidyverse: ", as.character(packageVersion("tidyverse")), "\n", sep = "")
cat("=====================================\n")

# -----------------------------
# 最后打印一个完成提示
# -----------------------------
cat("00_install_packages.R 跑完：依赖包已经准备好。\n")
