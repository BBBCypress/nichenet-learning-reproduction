# =========================================================
# 这个文件只负责画图
#
# 主分析那一步已经把关键结果存成了 RDS
# 这里直接把结果读回来，再整理成适合画图的矩阵
#
# 这一部分先看清楚两件事：
# 1. 图 3a 画的是什么
# 2. 图 3b 画的是什么
#
# 图 3a：
# 看 top ligands 的 activity 分数
#
# 图 3b：
# 看 top ligands 和 predicted target genes 之间的关系
# =========================================================

library(nichenetr)
library(tidyverse)
library(tibble)
library(ggplot2)

# -----------------------------
# 先确保画图输出目录存在
# -----------------------------
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 读取主分析结果
# -----------------------------
# readRDS():
# read = 读取
# RDS = R 的单对象保存格式
#
# 这里直接把上一步保存好的分析结果整个读回来
analysis_results <- readRDS("results/rds/analysis_results.rds")

# 从 analysis_results 这个列表里取出后面画图要用的对象
# $ 的作用可以先记成：
# “从一个列表/对象里取出某个命名元素”
ligand_activities <- analysis_results$ligand_activities
best_upstream_ligands <- analysis_results$best_upstream_ligands
active_ligand_target_links_df <- analysis_results$active_ligand_target_links_df
ligand_target_matrix <- analysis_results$ligand_target_matrix

# -----------------------------
# 准备图 3a 用的矩阵
# -----------------------------
# 这一步的目标很简单：
# 从 ligand_activities 里只保留 top ligands，
# 再把 aupr_corrected 整理成一个热图能直接使用的矩阵
vis_ligand_aupr <- ligand_activities %>%
  # filter():
  # 过滤出想保留的行
  # filter 这个词本身就是“筛选”
  filter(test_ligand %in% best_upstream_ligands) %>%

  # column_to_rownames():
  # 把某一列拿去做行名
  #
  # 这里把 test_ligand 变成行名，
  # 后面画图时每一行就对应一个 ligand
  column_to_rownames("test_ligand") %>%

  # select():
  # 只保留想看的列
  select(aupr_corrected) %>%

  # arrange():
  # 排序
  # 这里按 aupr_corrected 从低到高排
  arrange(aupr_corrected) %>%

  # as.matrix():
  # as = 转成
  # matrix = 矩阵
  # 这里把数据框转成矩阵，方便热图函数使用
  as.matrix(ncol = 1)

cat("vis_ligand_aupr head:\n")
print(head(vis_ligand_aupr))

cat("vis_ligand_aupr dim:\n")
# dim():
# dimension = 维度
# 看矩阵有多少行、多少列
print(dim(vis_ligand_aupr))

# -----------------------------
# 画图 3a：Ligand activity heatmap
# -----------------------------
# make_heatmap_ggplot():
# make = 生成
# heatmap = 热图
# ggplot = 返回 ggplot 风格图对象
#
# 这里画的是：
# 每个 top ligand 的 activity 分数
p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Prioritized ligands",
  "Ligand activity",
  legend_title = "AUPR",
  color = "darkorange"
) +
  theme(axis.text.x.top = element_blank())

# print():
# 在当前会话里把图显示出来
print(p_ligand_aupr)

# -----------------------------
# 把图 3a 保存成 PDF 和 PNG
# -----------------------------
# ggsave():
# gg = ggplot
# save = 保存
#
# 这是保存 ggplot 图最常用的函数之一
ggsave(
  filename = "results/plots/ligand_activity_heatmap.pdf",
  plot = p_ligand_aupr,
  width = 5,
  height = 8
)

ggsave(
  filename = "results/plots/ligand_activity_heatmap.png",
  plot = p_ligand_aupr,
  width = 5,
  height = 8,
  dpi = 300
)

# -----------------------------
# 准备图 3b 用的矩阵
# -----------------------------
# prepare_ligand_target_visualization():
# 先把 ligand-target 结果整理成适合热图展示的格式
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33
)

# -----------------------------
# 整理 ligand 的显示顺序
# -----------------------------
# intersect():
# 交集
# 看到 intersect，先理解成“两边共有的那些元素”
#
# 这里只保留：
# 既在 best_upstream_ligands 里，
# 也在 active_ligand_target_links 列名里出现的 ligand
order_ligands <- intersect(
  best_upstream_ligands,
  colnames(active_ligand_target_links)
) %>%
  # rev():
  # reverse = 反转顺序
  # 这里反转一下，是为了让图里的排列更符合想要的展示方向
  rev()

# -----------------------------
# 整理 target gene 的显示顺序
# -----------------------------
order_targets <- active_ligand_target_links_df$target %>%
  unique() %>%
  intersect(rownames(active_ligand_target_links))

# -----------------------------
# 构建最终画图矩阵
# -----------------------------
# 先按 order_targets 和 order_ligands 取出矩阵子集
# 再用 t() 做转置
#
# t():
# transpose = 转置
# 这个函数很常用，名字短但一定要记住
# 看到 t() 就先想“行列对调”
vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

cat("vis_ligand_target head:\n")
print(head(vis_ligand_target[, 1:min(5, ncol(vis_ligand_target))]))

cat("vis_ligand_target dim:\n")
print(dim(vis_ligand_target))

# -----------------------------
# 画图 3b：Predicted target genes heatmap
# -----------------------------
# 这里画的是：
# top ligands 和 predicted target genes 之间的关系强弱
p_ligand_target <- make_heatmap_ggplot(
  vis_ligand_target,
  "Prioritized ligands",
  "Predicted target genes",
  color = "purple",
  legend_title = "Regulatory potential"
) +
  scale_fill_gradient2(
    low = "whitesmoke",
    high = "purple"
  )

print(p_ligand_target)

# -----------------------------
# 把图 3b 保存成 PDF 和 PNG
# -----------------------------
ggsave(
  filename = "results/plots/predicted_target_genes_heatmap.pdf",
  plot = p_ligand_target,
  width = 10,
  height = 8
)

ggsave(
  filename = "results/plots/predicted_target_genes_heatmap.png",
  plot = p_ligand_target,
  width = 10,
  height = 8,
  dpi = 300
)

# -----------------------------
# 顺手把画图矩阵也存下来
# -----------------------------
# 这样做有两个好处：
# 1. 后面想复查画图输入时更方便
# 2. 如果以后想换别的画图方法，不用重新整理矩阵
write.csv(
  vis_ligand_aupr,
  "results/tables/vis_ligand_aupr_matrix.csv"
)

write.csv(
  vis_ligand_target,
  "results/tables/vis_ligand_target_matrix.csv"
)

cat("03_plot_results.R 跑完：图和画图矩阵都已经保存好。\n")
