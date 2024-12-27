###########PROJET FIRE PEDOLOGIE######

#########Charger les bibliothèques nécessaires
library(vegan)       # Pour PERMANOVA et analyses multivariées
library(cluster)     # Pour la CAH (Classification Hiérarchique Ascendante)
library(ape)         # Pour la PCoA
library(ggplot2)     # Pour les visualisations
library(FD)          # Pour calculer la distance de Gower
library(reshape2)    # Pour la mise en forme des données longues
library(FD)

## chemin vers le fichier CSV
setwd("C:/Users/bemma/OneDrive/Documents/Projet FIRE")

data <- read.csv("PEDO1.csv", sep = ";", stringsAsFactors = FALSE)

# Ajouter une colonne pour le parc
# Extraire le nom du parc (première partie de la colonne Sample)
data$Parc <- sub("_.*", "", data$Sample)

###ÉTAPE 1 : Classification Hiérarchique Ascendante
# Objectif : Visualiser la dissimilarité entre les carottes

# Calcule de la matrice de dissimilarité avec les colonnes d'épaisseur
dist_matrix <- gowdis(data[, grep("_thickness", colnames(data))])

# Appliquer l'algorithme UPGMA
cah_upgma <- hclust(dist_matrix, method = "average")

# Visualiser le dendrogramme
plot(cah_upgma, 
     main = "Dendrogramme - CAH avec UPGMA", 
     xlab = "Carottes", 
     ylab = "Dissimilarité")


#CAH par quadrat au sein d'un parc
# Pour le Domaine de Meric 
parc_select <- "DM"
data_parc <- subset(data, Parc == parc_select)

# Vérifier si les quadrats existent dans les données
if (!"Quadrat" %in% colnames(data_parc)) {
  stop("Les données doivent inclure une colonne 'Quadrat'.")
}

# Regrouper les données par quadrat et calculer la moyenne des épaisseurs
quadrat_summary <- aggregate(data_parc[, grep("_thickness", colnames(data_parc))],
                             by = list(Quadrat = data_parc$Quadrat),
                             mean)

# Calcule de la matrice de dissimilarité entre quadrats
dist_matrix_quadrat <- gowdis(quadrat_summary[, -1])  # Exclure la colonne 'Quadrat'

# Appliquer la CAH
cah_quadrat <- hclust(dist_matrix_quadrat, method = "average")

# Visualisation du dendrogramme pour les quadrats dans le parc sélectionné
plot(cah_quadrat, 
     main = paste("Dendrogramme - CAH entre Quadrats (Parc :", parc_select, ")"), 
     xlab = "Quadrats", 
     ylab = "Dissimilarité")


# Boucle pour analyser chaque parc
for (parc in unique(data$Parc)) {
  # Filtrer les données pour le parc
  data_parc <- subset(data, Parc == parc)
  
  # Vérifier si les quadrats existent dans les données
  if (!"Quadrat" %in% colnames(data_parc)) {
    stop("Les données doivent inclure une colonne 'Quadrat'.")
  }
  
  # Calcule des moyennes par quadrat
  quadrat_summary <- aggregate(data_parc[, grep("_thickness", colnames(data_parc))],
                               by = list(Quadrat = data_parc$Quadrat),
                               mean)
  
  # Calcule de la matrice de dissimilarité entre quadrats
  dist_matrix_quadrat <- gowdis(quadrat_summary[, -1])  # Exclure la colonne 'Quadrat'
  
  # Appliquation de la CAH
  cah_quadrat <- hclust(dist_matrix_quadrat, method = "average")
  
  # Visualiser le dendrogramme
  plot(cah_quadrat, 
       main = paste("Dendrogramme - CAH entre Quadrats (Parc :", parc, ")"), 
       xlab = "Quadrats", 
       ylab = "Dissimilarité")
}



###"#ÉTAPE 2 : Analyse en Coordonnées Principales (PCoA) 
# Objectif : Réduire la dimensionnalité et observer la répartition des carottes

# Calcule de la matrice de dissimilarité transformée (racine carrée s'il le faut)
dist_sqrt <- sqrt(as.matrix(dist_matrix))

# Effectuer la PCoA
pcoa_result <- cmdscale(dist_sqrt, k = 2, eig = TRUE)

# Extraire les coordonnées des axes
pcoa_axes <- as.data.frame(pcoa_result$points)
colnames(pcoa_axes) <- c("Axe1", "Axe2")

# Ajout des informations de parc aux résultats
pcoa_axes$Parc <- data$Parc

# Calcule de la proportion de variance expliquée
variance_expliquee <- pcoa_result$eig / sum(pcoa_result$eig)
cat("Variance expliquée par les axes:\n")
print(variance_expliquee)

# Visualisation des résultats de la PCoA avec ellipses
pcoa_plot <- ggplot(pcoa_axes, aes(x = Axe1, y = Axe2, color = Parc)) +
  geom_point(size = 3) +
  stat_ellipse() +
  labs(title = "PCoA avec Distance de Gower", x = "Axe 1", y = "Axe 2") +
  theme_minimal()
print(pcoa_plot)

# Corrélation avec les axes PCoA
envfit_result <- envfit(pcoa_axes[, c("Axe1", "Axe2")], data[, grep("_thickness", colnames(data))], permutations = 999)
print(envfit_result)

# Extraire les vecteurs significatifs et convertir en data.frame
vectors_df <- as.data.frame(envfit_result$vectors$arrows)
colnames(vectors_df) <- c("Axis1", "Axis2")  # Pour renommer les colonnes pour correspondre aux axes

# Ajouter les vecteurs au graphique PCoA
pcoa_plot_vectors <- pcoa_plot +
  geom_segment(data = vectors_df,
               aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
               arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  labs(title = "PCoA avec vecteurs de corrélation")
print(pcoa_plot_vectors)



#ÉTAPE 3 : PERMANOVA 
# Objectif : Tester l'effet des facteurs (Parc)

# Tester l'effet du parc sur les dissimilarités
permanova_result <- adonis2(dist_matrix ~ Parc, data = data, permutations = 999)
cat("Résultats de la PERMANOVA:\n")
print(permanova_result)

# ----- DISPERSION BÊTA -----
# Tester également la dispersion entre groupes
betadisper_result <- betadisper(dist_matrix, data$Parc)

# Tester la significativité des dispersions
anova_result <- anova(betadisper_result)
cat("Résultats de l'ANOVA sur la dispersion beta:\n")
print(anova_result)

# Visualisation de la dispersion
plot(betadisper_result, main = "Dispersion beta entre les parcs")

# Visualisation des distances moyennes avec des couleurs par groupe
boxplot(betadisper_result, main = "Distances moyennes aux centroids par parc", 
        col = rainbow(length(unique(data$Parc))))

###########Des ANALYSES COMPLÉMENTAIRES 

# 1. Analyse des corrélations entre les horizons
horizon_presence <- data[, c("O", "A", "T1", "T2", "T3", "J")]
cor_matrix <- cor(horizon_presence)

# Affichage la matrice de corrélation
cat("Matrice de corrélation entre horizons:\n")
print(cor_matrix)

# Visualisation la matrice de corrélation
cor_melted <- melt(cor_matrix)
cor_plot <- ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  labs(title = "Corrélations entre horizons", x = "Horizon", y = "Horizon") +
  theme_minimal()
print(cor_plot)

# 2. Comparaison des épaisseurs entre les parcs
data_long <- melt(data[, c("Parc", "O_thickness", "A_thickness", "T1_thickness", 
                           "T2_thickness", "T3_thickness", "J_thickness")], 
                  id.vars = "Parc", variable.name = "Horizon", value.name = "Epaisseur")

# Boxplots par parc et par horizon
boxplot_plot <- ggplot(data_long, aes(x = Parc, y = Epaisseur, fill = Horizon)) +
  geom_boxplot() +
  labs(title = "Comparaison des épaisseurs par parc et horizon", 
       x = "Parc", y = "Épaisseur (cm)") +
  theme_minimal()
print(boxplot_plot)
