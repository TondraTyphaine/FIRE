library(vegan)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

data <- read.csv("data_flo.csv")

# merge park and quadrat
pa_data <- data %>%
  mutate(park_quadrat = paste(Au, Quadrat, sep = "_"))

# presence/absence table
presence_absence <- pa_data %>%
  select(park_quadrat, Espece) %>%
  distinct() %>%  # rm duplicates
  mutate(presence = 1) %>%
  pivot_wider(names_from = Espece, values_from = presence, values_fill = list(presence = 0))

# change the type of presence_absnce from tibble
presence_absence <- as.data.frame(presence_absence)
rownames(presence_absence) <- presence_absence$park_quadrat
presence_absence <- presence_absence[, -1]

head(presence_absence)

# (forgot why i did that, im clearly not using it)
data_ji <- read.csv("presence_absence.csv", row.names = 1)

# Compute Jaccard distance and convert to similarity matrix
jaccard_dist <- vegdist(presence_absence, method = "jaccard")
jaccard_sim <- 1 - as.matrix(jaccard_dist)
jaccard_clust <- hclust(jaccard_dist, method = "average")

pheatmap(jaccard_sim, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         # labels_row = rownames(data_ji),
         # labels_col = colnames(data_ji),
         main = "Jaccard Index Heatmap for Park Quadrants",
         display_numbers = TRUE,  # optional
         color = colorRampPalette(c("white", "red"))(50))

plot(jaccard_clust,
     main = "Jaccard Dendrogram of Species per Park",
     xlab = "Parks",
     sub = "",
     cex = 0.9)

nmds <- metaMDS(presence_absence, distance = "jaccard")
plot(nmds, main = "NMDS of Park Quadrants")

## follow up with some more diversity data, shanon, richness, simpson and accumulation below too
## im using the vegan package for most of it

richness <- specnumber(presence_absence)
richness_df <- data.frame(Site = names(richness), Richness = richness)

ggplot(richness_df, aes(x = Site, y = Richness)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Species Richness by Site", x = "Site", y = "Number of Species")

shannon <- diversity(presence_absence, index = "shannon")
shannon_df <- data.frame(Site = names(shannon), Shannon = shannon)

ggplot(shannon_df, aes(x = Site, y = Shannon)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Shannon Diversity by Site", x = "Site", y = "Shannon Diversity")

simpson <- diversity(presence_absence, index = "simpson")
simpson_df <- data.frame(Site = names(simpson), Simpson = simpson)

ggplot(simpson_df, aes(x = Site, y = Simpson)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Simpson Diversity by Site", x = "Site", y = "Simpson Diversity")

# if i decide to use bray curtis instead of jaccard, i have that 
bray_curtis <- vegdist(presence_absence, method = "bray")
hclust_result <- hclust(bray_curtis)

plot(hclust_result, main = "Cluster Dendrogram of Sites", xlab = "Sites")


## overcomplicated code just to have accumulation curves per park
park_names <- unique(sub("_.*", "", rownames(presence_absence)))
accum_curves <- list()

# a loop to get the data per park
for (park in park_names) {
  park_data <- presence_absence[grep(paste0("^", park), rownames(presence_absence)), ]
  accum_curves[[park]] <- specaccum(park_data, method = "random")
}
# prep for plot
plot_data <- data.frame()
for (park in names(accum_curves)) {
  curve_data <- data.frame(
    Sites = accum_curves[[park]]$sites,
    Richness = accum_curves[[park]]$richness,
    SD = accum_curves[[park]]$sd,
    Park = park
  )
  plot_data <- rbind(plot_data, curve_data)
}

# plot de toto
ggplot(plot_data, aes(x = Sites, y = Richness, color = Park)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Species Accumulation Curves by Park",
       x = "Number of Sites",
       y = "Species Richness") +
  theme(legend.position = "right")

# im an idiot, the accumulation curve should be done with the abundance, that's the TODO for tomorrow, gotta sleep now