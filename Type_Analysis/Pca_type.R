library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(FactoMineR)
library(factoextra)

# Read the data
typeflo <- read.csv("typeflo.csv", row.names = 1)

# PCA for typeflo.csv

typeflo_numeric <- as.data.frame(apply(typeflo, 2, function(x) as.numeric(gsub(",", ".", x))))
rownames(typeflo_numeric) <- rownames(typeflo)


pca_typeflo <- rda(typeflo_numeric)
pca_scores_typeflo <- scores(pca_typeflo, display = c("sites", "species"))


pca_df_typeflo <- data.frame(
  PC1 = pca_scores_typeflo$sites[, 1],
  PC2 = pca_scores_typeflo$sites[, 2],
  Group = substr(rownames(typeflo_numeric), 1, 2)
)
main_plot <- ggplot(pca_df_typeflo, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pca_df_typeflo)), vjust = 1.5, hjust = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_fixed() +
  theme_minimal() +
  labs(title = "PCA of Plant Type Data",
       x = paste0("PC1 (", round(summary(pca_typeflo)$cont$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_typeflo)$cont$importance[2, 2] * 100, 1), "%)"))

pca_result_type <- prcomp(typeflo_numeric, scale = TRUE)

circle_plot <- fviz_pca_var(pca_result_type, col.var = "contrib", repel = TRUE,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             pointsize = 2, labelsize = 2) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

# save main plot in a png of size 20x20 cm
ggsave("main_plot.png", main_plot, width = 20, height = 20, units = "cm")
# save main plot in a png of size 20x20 cm
ggsave("circle_plot.png", circle_plot, width = 20, height = 20, units = "cm")


