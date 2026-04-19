# =========================================================
# 文件名: 03_plot_results.R
# 作用:
#   这个脚本负责把主分析结果整理成可视化图。
#
# 本脚本做的事情:
#   1. 读取主分析结果
#   2. 准备图 3a 所需的 ligand activity 矩阵
#   3. 绘制 ligand activity heatmap
#   4. 准备图 3b 所需的 ligand-target 矩阵
#   5. 绘制 predicted target genes heatmap
#   6. 保存图像到 results/plots/
#
# 说明:
#   这里保留了你原来代码的主要画图思路，
#   但去掉了控制台输出痕迹和重复代码。
# =========================================================

# -----------------------------
# 第 0 步：加载包
# -----------------------------
library(nichenetr)
library(tidyverse)
library(tibble)
library(ggplot2)

# -----------------------------
# 第 1 步：创建图像输出目录
# -----------------------------
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 第 2 步：读取主分析结果
# -----------------------------
analysis_results <- readRDS("results/rds/analysis_results.rds")

# 从列表中取出后面画图要用的对象
ligand_activities <- analysis_results$ligand_activities
best_upstream_ligands <- analysis_results$best_upstream_ligands
active_ligand_target_links_df <- analysis_results$active_ligand_target_links_df
ligand_target_matrix <- analysis_results$ligand_target_matrix

# -----------------------------
# 第 3 步：准备图 3a 的矩阵
# -----------------------------
# 逻辑：
# 从 ligand_activities 表里只保留 top ligands，
# 然后把 aupr_corrected 整理成一个 30 x 1 的矩阵
vis_ligand_aupr <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>%
  select(aupr_corrected) %>%
  arrange(aupr_corrected) %>%
  as.matrix(ncol = 1)

cat("vis_ligand_aupr head:\n")
print(head(vis_ligand_aupr))
cat("vis_ligand_aupr dim:\n")
print(dim(vis_ligand_aupr))

# -----------------------------
# 第 4 步：绘制图 3a
# -----------------------------
# make_heatmap_ggplot() 是 nichenetr 提供的辅助函数
# 用来快速画热图
p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  "Prioritized ligands",
  "Ligand activity",
  legend_title = "AUPR",
  color = "darkorange"
) +
  theme(axis.text.x.top = element_blank())

# 在当前会话中显示图
print(p_ligand_aupr)

# 保存图 3a
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
# 第 5 步：准备图 3b 的矩阵
# -----------------------------
# prepare_ligand_target_visualization() 的作用：
# 把 ligand-target 的长表整理成适合画热图的矩阵形式
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33
)

# order_ligands:
# 让图中的 ligand 顺序和高分 ligand 的顺序一致，再反转一下方便显示
order_ligands <- intersect(
  best_upstream_ligands,
  colnames(active_ligand_target_links)
) %>%
  rev()

# order_targets:
# 让 target gene 顺序和 active_ligand_target_links_df 中出现的顺序一致
order_targets <- active_ligand_target_links_df$target %>%
  unique() %>%
  intersect(rownames(active_ligand_target_links))

# 构建最终画图矩阵
# t() 的作用是转置矩阵
# 让 ligand 放在行，target gene 放在列，更符合原始作图习惯
vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

cat("vis_ligand_target head:\n")
print(head(vis_ligand_target[, 1:min(5, ncol(vis_ligand_target))]))
cat("vis_ligand_target dim:\n")
print(dim(vis_ligand_target))

# -----------------------------
# 第 6 步：绘制图 3b
# -----------------------------
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

# 在当前会话中显示图
print(p_ligand_target)

# 保存图 3b
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
# 第 7 步：保存画图矩阵
# -----------------------------
# 这样以后如果想在别的软件里继续画图，也方便导出
write.csv(
  vis_ligand_aupr,
  "results/tables/vis_ligand_aupr_matrix.csv"
)

write.csv(
  vis_ligand_target,
  "results/tables/vis_ligand_target_matrix.csv"
)

# -----------------------------
# 第 8 步：打印完成提示
# -----------------------------
cat("03_plot_results.R 运行完成：图像和画图矩阵已保存。\n")
