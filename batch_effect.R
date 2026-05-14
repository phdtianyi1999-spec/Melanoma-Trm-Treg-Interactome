table(sce_11_downstream$orig.ident)
table(sce_11_downstream$cell_type)
# If not installed, please run the following line
# install.packages(c("Seurat", "ggplot2", "patchwork", "gridExtra", "pals"))
getwd()

library(Seurat)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(pals)

# 1. Clean up potential spaces and fix the factor order
clean_names <- trimws(as.character(sce_11_downstream$orig.ident))
exact_samples <- sort(unique(clean_names))
sce_11_downstream$orig.ident <- factor(clean_names, levels = exact_samples)

# 2. Generate a [pure, unnamed] color vector (unname is key; no names, just the colors themselves)
n_samples <- length(exact_samples)
pure_colors <- unname(as.vector(polychrome(36))[1:n_samples])

# 3. Plotting: Remove 'cols' in Seurat and use ggplot2's scale_color_manual to force overwrite
p_combined <- DimPlot(sce_11_downstream, 
                      reduction = "umap", 
                      group.by = "orig.ident", 
                      pt.size = 0.5,
                      raster = FALSE) + 
  scale_color_manual(values = pure_colors) +  # <--- Use ggplot2 directly to force coloring
  ggtitle("UMAP: All Samples Combined") +
  theme_bw() + 
  theme(legend.position = "right", 
        legend.text = element_text(size = 9))

p_split <- DimPlot(sce_11_downstream, 
                   reduction = "umap", 
                   split.by = "orig.ident", 
                   ncol = 5, 
                   pt.size = 0.5,
                   raster = FALSE) + 
  scale_color_manual(values = pure_colors) +  # <--- Also force coloring
  ggtitle("UMAP: Split by Sample") +
  theme_bw() +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "gray90"),
        strip.text = element_text(size = 10, face = "bold")) 

# 4. Table extraction and final layout assembly
cell_counts <- as.data.frame(table(sce_11_downstream$orig.ident))
colnames(cell_counts) <- c("Sample", "Cell_Count")

table_plot <- tableGrob(cell_counts, 
                        rows = NULL, 
                        theme = ttheme_minimal(base_size = 10,
                                               core = list(bg_params = list(fill = c("#f0f0f0", "white"))),
                                               colhead = list(bg_params = list(fill = "#525252"),
                                                              fg_params = list(col = "white", fontface = "bold"))))

top_row <- p_combined + wrap_elements(panel = table_plot) + plot_layout(widths = c(3, 1))
final_layout <- top_row / p_split + plot_layout(heights = c(1.2, 2))

# 5. Export high-quality image
ggsave(filename = "Batch_Check_Final.png", 
       plot = final_layout, 
       width = 18, 
       height = 16, 
       dpi = 300,
       bg = "white")

# 5. Export high-quality PDF vector graphic
ggsave(filename = "Batch_Check_Final.pdf", 
       plot = final_layout, 
       width = 18, 
       height = 16,
       bg = "white") # Keep white background to prevent turning black in some readers