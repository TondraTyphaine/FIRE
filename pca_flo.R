library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(FactoMineR)
library(factoextra)

# Read the data
flo <- read.csv("flo.csv", row.names = 1)

flo_numeric <- as.data.frame(apply(flo[, -1], 2, as.numeric))
rownames(flo_numeric) <- rownames(flo)

pca_flo <- rda(flo_numeric)
pca_scores <- scores(pca_flo, display = c("sites", "species"))


pca_df <- data.frame(
  PC1 = pca_scores$sites[, 1],
  PC2 = pca_scores$sites[, 2],
  Group = substr(rownames(flo_numeric), 1, 2)
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pca_df)), vjust = 1.5, hjust = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_equal() +
  theme_minimal() +
  labs(title = "PCA of Flora Data",
       x = paste0("PC1 (", round(summary(pca_flo)$cont$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_flo)$cont$importance[2, 2] * 100, 1), "%)"))

# save ggplot in a png file with max size
ggsave("pca_flo.png", width = 20, height = 20, units = "cm")

