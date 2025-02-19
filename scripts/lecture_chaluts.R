
switch.datetime <- function(x) { paste(vect[x],collapse='') }

as.num <- function(x) { x <- as.numeric(as.character(x)) }

change.res <- function(x, new.res) { new.res * round(x/new.res) }

library(Rcpp)
library(R.matlab)
library(ggplot2)
library(matlab)
library(ggpubr)
library(ggthemes)

library(dplyr)

library(stringr)
library(scales)
library(reshape2)
library(data.table)
library(gridExtra)

library(devtools)

library(raster)
library(ncdf4)
library(rhdf5)
library(stars)
library(PupillometryR)

library(sf)
library("readxl")
library(ggrepel)

library(FactoMineR)
library(statmod)
library(vegan)

coco <- read.csv("C:/aurore/eddy/data_acoustic/correspondance.csv", sep = ";")
vect <- rep("NA", 58)
vect[as.numeric(coco$symbol)] <- coco$corres

list_camp_maracas = c('M7B', 'M7C', 'M7D')

setwd('D:/maracas')
rotate <- function(x) t(apply(x, 2, rev))

`%ni%` = Negate(`%in%`)

###################################################

df_bathy <- read.table(file = "D:/maracas/data_maracas/xyz_500m.txt",
                       sep = ",", dec = ".", head = FALSE)
names(df_bathy) <- c("lon", "lat", "depth")

df_bathy <- df_bathy %>% filter(lon %between% c(167.8, 168.6) &
                                  lat %between% c(-23.7, -22.7) &
                                  depth <= 0)

###################################################
#################### lecture donnees #########################
###################################################
# EI_M7B_90_160m

chalut_df_A <- read.csv2("data_maracas/chaluts/Fiche_all_maracas.csv")
chalut_df_A <- chalut_df_A %>%  dplyr::filter(trip_id != 'MARACAS7A') 


ggplot(chalut_df_A, aes (x =  Lon_SE , y =Lat_SS, color = trip_id  )) + geom_point()

df_bathy
ras_bathy <- raster(rotate(rotate(rotate(acast(df_bathy,    lon  ~ lat,
                                    value.var = 'depth',
                                    fun.aggregate = mean)))), 
                   xmn = min(df_bathy$lon),xmx = max(df_bathy$lon),
                   ymn = min(df_bathy$lat), ymx = max(df_bathy$lat), 
                   CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(ras_bathy)

df_xy <- chalut_df_A[,c('station_no', 'Lon_SE', 'Lat_SE')]
names(df_xy)[c(2,3)] <- c("x", "y")
coordinates(df_xy) <- ~ x + y
proj4string(df_xy) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

chalut_df_A$bathy <- raster::extract(ras_bathy, df_xy)

#################### lecture species composition

chalut_df_B <- read_xlsx("data_maracas/chaluts/MARACAS7_All micronekton_template_analysis_forages.xlsx")
chalut_df_B <- chalut_df_B %>%  dplyr::filter(`trip name` != 'Maracas7A') 

names(chalut_df_B) <- str_replace_all(names(chalut_df_B), " ", "_")
names(chalut_df_B)
chalut_df_B <- chalut_df_B %>% 
  dplyr::select(trip_name, station, set_no,scientific_name,  nb_individuals, weight_1_g)

chalut_df_B$espace_presence <- str_detect(chalut_df_B$scientific_name, "\\s")
chalut_df_B$genre <- str_extract(chalut_df_B$scientific_name, pattern = ".+(?=\\s)")
chalut_df_B$genre <- ifelse(chalut_df_B$espace_presence  == TRUE,
                            chalut_df_B$genre,
                            chalut_df_B$scientific_name)

chalut_df_B <- chalut_df_B %>%  dplyr::filter(genre != 'Rubbish (human' &
                                                genre != 'Unidentified' &
                                                genre != 'Unidentified gelatinous') 

# df_genre <- data.frame(genre =  unique(chalut_df_B$genre)[order(unique(chalut_df_B$genre))])
# 
# library('taxize')
# dede <- sapply(df_genre$genre, tax_name, get = c("phylum"), db = "ncbi")
# df_genre$phylum <- unlist(dede[3,])
# 
# df_genre[df_genre$genre == "Abyla", "phylum"] <- 'Cnidaria'
# df_genre[df_genre$genre == "Bolitaeninae", "phylum"] <- 'Mollusca'
# df_genre[df_genre$genre == "Ctenophora", "phylum"] <- 'Cnidaria'
# df_genre[df_genre$genre == "Eupogonesthes", "phylum"] <- 'Chordata'
# df_genre[df_genre$genre == "Hydromyles", "phylum"] <- 'Mollusca'
# df_genre[df_genre$genre == "Ommastrephinae", "phylum"] <- 'Mollusca'
# df_genre[df_genre$genre == "Salpinae", "phylum"] <- 'Chordata'
# df_genre[df_genre$genre == "Scintillosergia", "phylum"] <- 'Arthropoda'
# df_genre[df_genre$genre == "Shrimp", "phylum"] <- 'Arthropoda'
# df_genre[df_genre$genre == "Taoniinae", "phylum"] <- 'Mollusca'
# 
# table(df_genre$phylum)
# df_genre[is.na(df_genre$phylum),]
# save(file = "data_maracas/chaluts/df_genre.Rdata", df_genre)

load(file = "data_maracas/chaluts/df_genre.Rdata")

chalut_df_B <- merge(chalut_df_B, df_genre)
table(chalut_df_B$phylum)

chalut_df_B_bis <- chalut_df_B %>% 
  dplyr::filter(espace_presence == TRUE)
chalut_df_B_bis$scientific_name <- str_replace_all(chalut_df_B_bis$scientific_name,
                                                   " ", "_")
chalut_df_B_bis$species <- sapply(strsplit(chalut_df_B_bis$scientific_name, split='_', fixed=TRUE), function(x) (x[2]))

chalut_df_B_bis <- chalut_df_B_bis %>% 
  dplyr::filter(species != 'sp.') %>% 
  dplyr::select(-trip_name )


################### merge 
chalut_df_A <- chalut_df_A %>% 
  dplyr::select( trip_id, station_no, set_no , day_night_id, set_end_time_UTC , 
                 Lat_SE, Lon_SE, trawl_duration , observed_depth_avg , 
                 volume_filtered_observed_m3, bathy)
names(chalut_df_B_bis)[names(chalut_df_B_bis) == "station"] <-  "station_no"

df_all <- merge(chalut_df_B_bis, chalut_df_A)
df_all$nb_ind_norm <- df_all$nb_individuals/df_all$volume_filtered_observed_m3*1000
df_all <- df_all %>% 
  dplyr::filter(!is.na(nb_ind_norm))

n_distinct(df_all$scientific_name)


df_all_richesse <- df_all %>% 
  dplyr::group_by(trip_id, set_no,day_night_id, set_end_time_UTC,
                  Lon_SE, Lat_SE , observed_depth_avg , 
                  volume_filtered_observed_m3, bathy) %>% 
  dplyr::summarise(n_sp = n_distinct(scientific_name))
df_all_richesse$vertical_layer <- base::cut(df_all_richesse$observed_depth_avg,
                                            c(0, 200, 800 ))

pl_R = ggplot() +  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth), color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  scale_color_viridis_c(name = "Species number") +
  geom_point(data = df_all_richesse, aes(x = Lon_SE, y = Lat_SE, color = n_sp),
             size = 2) +
  facet_grid( vertical_layer ~ day_night_id) 
pl_R
ggsave(filename ="D:/maracas/R_figures/N_species_chaluts.jpg",
       width = 3, height = 2.3, 
       scale = 3, plot= pl_R)

pl_R = ggplot() +  stat_contour(data = df_bathy,
                                aes(x = lon, y = lat, z = depth), color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  scale_color_viridis_b(name = "Species number") +
  geom_point(data = df_all_richesse, aes(x = Lon_SE, y = Lat_SE)) +
  geom_label_repel(data = df_all_richesse, aes(x = Lon_SE, y = Lat_SE, label = set_no) )+
  facet_grid( vertical_layer ~ day_night_id) 
pl_R
ggsave(filename ="D:/maracas/R_figures/N_species_chaluts_position.jpg",
       width = 3, height = 2.5, 
       scale = 3, plot= pl_R)

pl_R = ggplot() +  stat_contour(data = df_bathy,
                                aes(x = lon, y = lat, z = depth), color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                           limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  scale_color_viridis_c(name = "Bathymetry (m)") +
  geom_point(data = df_all_richesse, aes(x = Lon_SE, y = Lat_SE, color = bathy),
             size = 3) +
  facet_grid( vertical_layer ~ day_night_id) 
pl_R
ggsave(filename ="D:/maracas/R_figures/N_species_chaluts_bathy.jpg",
       width = 3, height = 2.5, 
       scale = 3, plot= pl_R)




##### distance aux monts sous marins 

load("D:/maracas/data_maracas/env_stack.RData")

df_chal_sf <- st_as_sf(df_all_richesse, coords = c("Lon_SE", "Lat_SE"), 
                       crs = 4326)


dist_star <- st_as_stars(env_stack$DIS_SUMMIT)
dist_star <- st_transform(dist_star, crs = 4326) # project to same crs

i <- st_extract(dist_star, df_chal_sf)
df_all_richesse$dist_summ <- i$DIS_SUMMIT

ggplot() +  
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth), color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df_all_richesse, aes(x = Lon_SE, y = Lat_SE,
                                         color = dist_summ/1000  ),
             size = 2) + scale_color_viridis_c()



lm_trawl_depth_quasi <- glm(data = df_all_richesse,
                     n_sp ~ observed_depth_avg,
                     family = quasipoisson)
lm_trawl_depth_pois <- glm(data = df_all_richesse,
                            n_sp ~ observed_depth_avg,
                            family = poisson)
null.model <- glm(n_sp ~ 1, data = df_all_richesse, family = quasipoisson)
anova(null.model, lm_trawl_depth_quasi,lm_trawl_depth_pois, test = "Chisq")

df_newdata <- data.frame(observed_depth_avg = seq(min(df_all_richesse$observed_depth_avg),
                                                  max(df_all_richesse$observed_depth_avg), 50))
df_newdata$predict <- predict(lm_trawl_depth_quasi, newdata = df_newdata, "response")
df_newdata$predict_se <- predict(lm_trawl_depth_quasi, newdata = df_newdata, "response",  se.fit = TRUE)$se.fit


df_all_richesse$residuals <- qresid(lm_trawl_depth_quasi)

hist(df_all_richesse$residuals )
summary(lm_trawl_depth_quasi)
# plot(lm_trawl_depth)

df_all_richesse$dist_summ2 <- df_all_richesse$dist_summ/1000
dim(df_all_richesse)
lm_trawl_dist <- glm(data = df_all_richesse,
                     residuals ~ dist_summ2,
                     family = gaussian)
summary(lm_trawl_dist)
df_newdata2 <- data.frame(dist_summ2 = seq(min(df_all_richesse$dist_summ2),
                                                  max(df_all_richesse$dist_summ2), 
                                           length.out = 50))
df_newdata2$predict <- predict(lm_trawl_dist, newdata = df_newdata2, "response")
df_newdata2$predict_se <- predict(lm_trawl_dist, newdata = df_newdata2, "response",  se.fit = TRUE)$se.fit


a = ggplot(df_all_richesse, aes(x = observed_depth_avg, y = n_sp)) + geom_point() +
  xlab("Trawling depth (m)") + ylab("Number of species") + 
  theme_minimal() +
  geom_ribbon(data = df_newdata, aes(x = observed_depth_avg, y = predict, 
                                     ymin = predict - predict_se,
                                     ymax = predict + predict_se), fill = "grey90", color = NA, alpha = 0.5) +
  geom_line(data = df_newdata, aes(x = observed_depth_avg, y = predict), color = "darkblue") +
  
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90')) 
a

b = ggplot(df_all_richesse, aes(x = dist_summ2, y = residuals)) + geom_point() +
  xlab("Distance to summit (km)") + ylab("Quantile residuals") + 
  theme_minimal() +
  geom_ribbon(data = df_newdata2, aes(x = dist_summ2, y = predict, 
                                      ymin = predict - predict_se,
                                      ymax = predict + predict_se), fill = "grey90", color = NA, alpha = 0.5) +
  geom_line(data = df_newdata2, aes(x = dist_summ2, y = predict), color = "darkblue") +
  
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90')) 
b



d = ggarrange(a, b, nrow = 2)
d
ggsave(filename ="D:/maracas/R_figures/chaluts_relations_bathy_distance_summit.jpg",
       width = 2, height = 2, 
       scale = 3, plot= d)


##### classif 
df_all <- df_all[df_all$scientific_name != "Pyrosoma_atlanticum", ]
trawl_sp_matrice <- reshape2::acast(data = df_all,
                                    set_no  ~ scientific_name,
                                    value.var = 'nb_ind_norm',
                                    fun.aggregate = mean)
trawl_sp_matrice[is.na(trawl_sp_matrice)] <- 0
dim(trawl_sp_matrice)
# trawl_sp_matrice_bray <- vegan::vegdist(trawl_sp_matrice, method="bray")

########### try hellinger pca 
trawl_sp_matrice_hell <- decostand(trawl_sp_matrice, method = 'hellinger')
trawl_sp_matrice_dist <- dist(trawl_sp_matrice_hell)

# res_acp_taxo <- PCA(trawl_sp_matrice_hell, scale.unit = FALSE, 
#                     graph = FALSE)
# 
# df_pca <- as.data.frame( res_acp_taxo$ind$coord)
# df_pca$set_no <- rownames(df_pca)
# df_pca <- merge(df_all_richesse, df_pca)
# 
# ggplot(df_pca, aes(x = Dim.1, y = Dim.2, color = observed_depth_avg)) + 
#   geom_hline(yintercept = 0, color = 'grey70') + 
#   geom_vline(xintercept = 0, color = 'grey70') +
#   geom_point() + 
#   ylab("PC2 (12%)") +  xlab("PC1 (21%)") + theme_classic()   +
#   theme(panel.background = element_rect(fill = 'white', color = 'black'),
#         panel.grid.major = element_line(color = 'grey90')) 

############ clustering  

library(NbClust)
NbClust(data = as.matrix(trawl_sp_matrice_hell), 
        distance = 'euclidean', 
        index = c("all"),
        method = "complete")
# Only frey, mcclain, cindex, sihouette and dunn can be computed.

############ 
hclust_avg <- stats::hclust((trawl_sp_matrice_dist), method = 'complete')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 8)
df_clus <- as.data.frame(cut_avg)
df_clus$set_no <- rownames(df_clus)

df_clus <- merge(df_clus, df_all_richesse)
# ggplot(df_clus, aes(x = Axis.1, y = Axis.3, color = as.factor(clus) )) + 
#   geom_hline(yintercept = 0, color = 'grey70') + 
#   geom_vline(xintercept = 0, color = 'grey70') +
#   geom_point() + 
#   ylab("PC2 (12%)") +  xlab("PC1 (21%)") + theme_classic()   +
#   theme(panel.background = element_rect(fill = 'white', color = 'black'),
#         panel.grid.major = element_line(color = 'grey90')) 


pl_R <- ggplot() +  stat_contour(data = df_bathy,
                         aes(x = lon, y = lat, z = depth), color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df_clus, aes(x = Lon_SE, y = Lat_SE, color = as.factor(cut_avg   )   ),
             size = 2) +
  scale_color_manual(values = c("#000000", "#FFB300", "#56B4E9", "#00AE7E",
                                 "#9E98C9", "#0072B2", "#FF7504", "#FB59FB"),
                     name = 'Clusters') +
  facet_grid( vertical_layer ~ day_night_id) 

pl_R
ggsave(filename ="D:/maracas/R_figures/N_species_chaluts_cluster.jpg",
       width = 2.5, height = 2.5, 
       scale = 3, plot= pl_R)



df_clus <- as.data.frame(cut_avg)
df_clus$set_no <- rownames(df_clus)

df_all2 <- merge(df_all, df_clus)

a = ggplot(df_all2[df_all2$cut_avg %in% c(2, 7), ],
       aes(x = scientific_name, y = nb_ind_norm, color = as.factor(cut_avg) )) + 
  geom_point() + coord_flip() + 
  scale_color_manual(values = c("#FFB300", "#FF7504"),
                     name = 'Clusters') +
  facet_grid(phylum~cut_avg,  space = "free_y",scales = 'free')
a
ggsave(filename ="D:/maracas/R_figures/N_species_chaluts_cluster_composition_1.jpg",
       width = 2.5, height = 4, 
       scale = 3, plot= a)
unique(df_all2$scientific_name)[order(unique(df_all2$scientific_name))]

b = ggplot(df_all2[df_all2$cut_avg %ni% c(2, 7), ],
           aes(x = scientific_name, y = nb_ind_norm, color = as.factor(cut_avg) )) + 
  geom_point() + coord_flip() + 
  scale_color_manual(values = c("#000000","#56B4E9",  "#00AE7E",
                                "#9E98C9", "#0072B2", "#FB59FB"),
                     name = 'Clusters') +
  facet_grid(phylum~cut_avg,  space = "free_y",scales = 'free')
b
ggsave(filename ="D:/maracas/R_figures/N_species_chaluts_cluster_composition_2.jpg",
       width = 2.5, height = 3, 
       scale = 3, plot= b)
unique(df_all2$scientific_name)[order(unique(df_all2$scientific_name))]

df_smal_ass <- df_all2 %>% 
  dplyr::group_by(cut_avg,  scientific_name) %>% 
  dplyr::summarise(mean_abun = mean(nb_ind_norm )) %>% 
  arrange(desc(mean_abun)) %>% 
  dplyr::slice(1:10)
dim(df_smal_ass)
df_phyl  <- df_all2 %>% 
  dplyr::group_by(phylum, scientific_name) %>% 
  dplyr::summarise(mean_abun = mean(nb_ind_norm )) %>% 
  dplyr::select(-mean_abun)
df_smal_ass <- merge(df_smal_ass, df_phyl)

df_smal_ass$scientific_name <- str_replace(df_smal_ass$scientific_name, "_", " ")
c = ggplot(df_smal_ass, aes(x = scientific_name, y =mean_abun, 
                            color = as.factor(cut_avg) )) + 
  geom_segment(data = df_smal_ass, aes(x = scientific_name, y = 0.003,
                                       xend = scientific_name, yend = mean_abun, 
                                       color = as.factor(cut_avg) )) +
  geom_point() + coord_flip() + 
  scale_color_manual(values = c("#000000", "#FFB300", "#56B4E9", "#00AE7E",
                                "#9E98C9", "#0072B2", "#FF7504", "#FB59FB"),
                     name = 'Clusters') +
  xlab('Scientific name') + ylab('Mean abundance') +
  facet_grid(phylum ~ cut_avg,  space = "free_y",scales = 'free')
c
ggsave(filename ="D:/maracas/R_figures/N_species_chaluts_cluster_composition_3.jpg",
       width = 3.5, height = 3, 
       scale = 3, plot= c)

############ pcoa  

library(ape)
pcoa_bray <- pcoa(sqrt(trawl_sp_matrice_bray))

eig.val_taxo <- (pcoa_bray$values)
head(eig.val_taxo)

# df_pcoa
df_pcoa <- as.data.frame(pcoa_bray$vectors[, 1:5])
df_pcoa$set_no <- rownames(df_pcoa)

df_pcoa <- merge(df_all_richesse, df_pcoa)

ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, color = observed_depth_avg)) + 
  geom_hline(yintercept = 0, color = 'grey70') + 
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point() + 
  ylab("PC2 (12%)") +  xlab("PC1 (21%)") + theme_classic()   +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90')) 
