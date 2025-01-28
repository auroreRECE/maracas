
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


coco <- read.csv("C:/aurore/eddy/data_acoustic/correspondance.csv", sep = ";")
vect <- rep("NA", 58)
vect[as.numeric(coco$symbol)] <- coco$corres

list_camp_maracas = c('M7B', 'M7C', 'M7D')

###################################################
#################### lecture donnees #########################
###################################################
# EI_M7B_90_160m

rep_resolution = '160m'

for (i in 1:3){
  print(list_camp_maracas[i])
  print(rep_resolution)
  
  file <- paste('D:/maracas/data_maracas/EI_',list_camp_maracas[i],
                "_90_", rep_resolution, ".mat",sep='')
  Sv <- h5read(file,"Sv_surface")
  Sa <- h5read(file,"Sa_surface")
  dep <- h5read(file,"depth_surface")
  night <- h5read(file,"Night1Sunrise2Day3Sunset4")
  date <- as.data.frame(h5read(file,"DateTime"))
  date2 <-  data.frame(date = apply(date, 1, switch.datetime))
  time_incr <-  data.frame(t_inc = h5read(file,"Time")[1,])
  lat <- h5read(file,"Latitude")
  lon <- h5read(file,"Longitude")
  solar <- h5read(file,"SolarAltitude")
  depth_bottom  <- h5read(file,"depth_bottom")
  
  
  dim1 <- dim(Sv)[2]
  dim_prof <- dim(Sv)[1]

  time <- c(1:dim1)

  clean_mas <- h5read(file,"Mask_clean")
  clean_mas <- t(clean_mas[1:dim_prof,1:dim1,1])

  df_night <- data.frame(time = time, moment = night[1,])
  df_night$moment <- as.factor(df_night$moment)
  levels(df_night$moment) <- c("Night", "Sunrise", "Day", "Sunset")

  df_night$date <- date2$date
  df_night <- data.table(df_night)
  df_night$date <- as.POSIXct(strptime(df_night$date,
                                       format = "%Y/%m/%d-%H:%M:%S",
                                       tz = 'UTC'))
  df_night$lon <- lon[1,]
  df_night$lat <- lat[1,]
  df_night$solar <- solar[1,]
  df_night$depth_bottom <- depth_bottom[1,,1]
  # df_night$date <- as.POSIXct(format(df_night$date, tz = my_tz),
  # tz = my_tz)

  ## Sv
  Sv_38 <- t(Sv[1:dim_prof,1:dim1,1])
  colnames(Sv_38) <- - dep[1:dim_prof]
  rownames(Sv_38) <- time
  Sv_38 <- Sv_38 *  clean_mas
  df_Sv_38 <- reshape2::melt(Sv_38, na.rm = TRUE)
  names(df_Sv_38) <- c("time", "depth", "Sv")
  df_Sv_38$freq <- "38"

  Sv_70 <- t(Sv[1:dim_prof,1:dim1,2])
  colnames(Sv_70) <- - dep[1:dim_prof]
  rownames(Sv_70) <- time
  Sv_70 <- Sv_70 *  clean_mas
  df_Sv_70 <- reshape2::melt(Sv_70, na.rm = TRUE)
  names(df_Sv_70) <- c("time", "depth", "Sv")
  df_Sv_70$freq <- "70"

  Sv_120 <- t(Sv[1:dim_prof,1:dim1,3])
  colnames(Sv_120) <- - dep[1:dim_prof]
  rownames(Sv_120) <- time
  Sv_120 <- Sv_120 *  clean_mas
  df_Sv_120 <- reshape2::melt(Sv_120, na.rm = TRUE)
  names(df_Sv_120) <- c("time", "depth", "Sv")
  df_Sv_120$freq <- "120"

  Sv_200 <- t(Sv[1:dim_prof,1:dim1,4])
  colnames(Sv_200) <- - dep[1:dim_prof]
  rownames(Sv_200) <- time
  Sv_200 <- Sv_200 *  clean_mas
  df_Sv_200 <- reshape2::melt(Sv_200, na.rm = TRUE)
  names(df_Sv_200) <- c("time", "depth", "Sv")
  df_Sv_200$freq <- "200"

  df_Sv <- rbind(df_Sv_38, df_Sv_70, df_Sv_120, df_Sv_200)
  df_Sv$freq <- as.factor(df_Sv$freq)
  df_Sv <- data.table(df_Sv)
  df_Sv$sv <- 10^(df_Sv$Sv/10)

  ## Sa
  Sa_38 <- t(Sa[1:dim_prof,1:dim1,1])
  colnames(Sa_38) <- - dep[1:dim_prof]
  rownames(Sa_38) <- time
  Sa_38 <- Sa_38 *  clean_mas
  df_Sa_38 <- reshape2::melt(Sa_38, na.rm = TRUE)
  names(df_Sa_38) <- c("time", "depth", "Sa")
  df_Sa_38$freq <- "38"

  Sa_70 <- t(Sa[1:dim_prof,1:dim1,2])
  colnames(Sa_70) <- - dep[1:dim_prof]
  rownames(Sa_70) <- time
  Sa_70 <- Sa_70 *  clean_mas
  df_Sa_70 <- reshape2::melt(Sa_70, na.rm = TRUE)
  names(df_Sa_70) <- c("time", "depth", "Sa")
  df_Sa_70$freq <- "70"

  Sa_120 <- t(Sa[1:dim_prof,1:dim1,3])
  colnames(Sa_120) <- - dep[1:dim_prof]
  rownames(Sa_120) <- time
  Sa_120 <- Sa_120 *  clean_mas
  df_Sa_120 <- reshape2::melt(Sa_120, na.rm = TRUE)
  names(df_Sa_120) <- c("time", "depth", "Sa")
  df_Sa_120$freq <- "120"

  Sa_200 <- t(Sa[1:dim_prof,1:dim1,4])
  colnames(Sa_200) <- - dep[1:dim_prof]
  rownames(Sa_200) <- time
  Sa_200 <- Sa_200 *  clean_mas
  df_Sa_200 <- reshape2::melt(Sa_200, na.rm = TRUE)
  names(df_Sa_200) <- c("time", "depth", "Sa")
  df_Sa_200$freq <- "200"


  df_Sa <- rbind(df_Sa_38, df_Sa_70, df_Sa_120, df_Sa_200)
  df_Sa$freq <- as.factor(df_Sa$freq)
  df_Sa <- data.table(df_Sa)

  ## merge
  df_Sv <- merge(df_Sv, df_Sa, by = c("time", "depth", "freq"))
  df_Sv <- merge(df_Sv, df_night, by = c("time"))

  # df_Sv <- df_Sv %>% filter(depth <= -5)

  df_Sv$camp <- list_camp_maracas[i]
nana = paste0("df_tot_4freqMAR_90_", rep_resolution)


  if(list_camp_maracas[i] == "M7B"){
    assign(nana, df_Sv, .GlobalEnv)
  } else {li <- list(get(nana), df_Sv)
    temp <- rbindlist(li)
    assign(nana, temp, .GlobalEnv)}
print(nana)
}

# df_tot_4freqMAR_90_160m
# df_tot_4freqMAR_90_500m
# df_tot_4freqMAR_90_1000m

save(file = 'D:/maracas/data_maracas/df_tot_4freqMAR_90_160m.Rdata',
     df_tot_4freqMAR_90_160m)

save(file = 'D:/maracas/data_maracas/df_tot_4freqMAR_90_500m.Rdata',
     df_tot_4freqMAR_90_500m)

save(file = 'D:/maracas/data_maracas/df_tot_4freqMAR_90_1000m.Rdata',
     df_tot_4freqMAR_90_1000m)

mean(diff(df_tot_4freqMAR_90_160m$date))
mean(diff(df_tot_4freqMAR_90_500m$date))
mean(diff(df_tot_4freqMAR_90_1000m$date))


# dede = df_tot_4freqMAR_90_160m %>% 
#   dplyr::filter(freq == "38") 
# 
# dede = dede %>% 
#   dplyr::group_by(camp, date, lon, lat) %>% 
#   summarise(bottom_depth = mean(depth_bottom, na.rm = T))
# dede$date = as.character(dede$date)
# 
# write.csv2(dede, file = 'C:/aurore/MARACAS_bottom_depth.csv',
#            row.names = FALSE)

###################################################
#################### creation netcdf #########################
###################################################

load(file = "D:/maracas/data/df_tot_4freqMAR_90_1000m.Rdata") 
df_tot_4freqMAR$date_julien_UTC <- as.num(julian(df_tot_4freqMAR$date))


do.netcdf = function(rep_camp = 'M7B'){
  print(rep_camp)
  temp = df_tot_4freqMAR %>% dplyr::filter(camp == rep_camp)
  
  var_38 <- temp %>% dplyr::filter(freq == "38")
  mat_38 <-  acast(var_38, date_julien_UTC ~ depth, value.var = "Sv", mean)
  
  
  var_70 <- temp %>% dplyr::filter(freq == "70")
  mat_70 <-  acast(var_70, date_julien_UTC ~ depth, value.var = "Sv", mean)
  mat_compl_70 <- matrix(nrow = nrow(mat_38), ncol = ncol(mat_38)-ncol(mat_70),
                         dimnames = list(rownames(mat_38),
                                         colnames(mat_38)[1:(ncol(mat_38)-ncol(mat_70))]),
                         data = NA)
  mat_ent_70 <- cbind(mat_compl_70, mat_70)
  
  
  
  var_120 <- temp %>% dplyr::filter(freq == "120")
  mat_120 <-  acast(var_120, date_julien_UTC ~ depth, value.var = "Sv", mean)
  mat_compl_120 <- matrix(nrow = nrow(mat_38), ncol = ncol(mat_38)-ncol(mat_120),
                          dimnames = list(rownames(mat_38),
                                          colnames(mat_38)[1:(ncol(mat_38)-ncol(mat_120))]),
                          data = NA)
  mat_ent_120 <- cbind(mat_compl_120, mat_120)
  
  
  
  
  var_200 <- temp %>% dplyr::filter(freq == "200")
  mat_200 <-  acast(var_200, date_julien_UTC ~ depth, value.var = "Sv", mean)
  mat_compl_200 <- matrix(nrow = nrow(mat_38), ncol = ncol(mat_38)-ncol(mat_200),
                          dimnames = list(rownames(mat_38),
                                          colnames(mat_38)[1:(ncol(mat_38)-ncol(mat_200))]),
                          data = NA)
  mat_ent_200 <- cbind(mat_compl_200, mat_200)
  
  
  
  df_complement <- var_38 %>% 
    dplyr::group_by(date, date_julien_UTC) %>%
    summarize(lon = mean(lon),
              lat = mean(lat))
  
  
  date_char <-  as.character(df_complement$date)
  lon <-  df_complement$lon
  lat <-  df_complement$lat
  
  depth <- as.num(colnames(mat_38))
  date_julien <-  as.num(df_complement$date_julien_UTC)
  
  ntime <- length(date_julien)
  ndepth <- length(depth)
  
  freq <- c(1:4)
  
  depth     <- ncdim_def("depth", units = "meters", 
                         depth)
  time_jul  <- ncdim_def("time_dec" , 
                         units = "Decimal days since 1970-01-01 00:00:00", 
                         date_julien)
  freq      <- ncdim_def("freq" , units = "4 frequences, 1:38kHz, 2:70kHz, 
                         3:120kHz, 4:200kHz", 
                         freq)
  
  
  dimnchar <- ncdim_def("nchar", "", 1:50, create_dimvar=FALSE)
  var.sv  <- ncvar_def( "Sv"  , "", list(time_jul, depth, freq))
  
  
  var.lon     <- ncvar_def( "lon"    , "", list(time_jul))
  var.lat     <- ncvar_def( "lat"    , "", list(time_jul))
  var.date.char  <- ncvar_def( "date" , "date UTC character", list(dimnchar, time_jul), 
                               prec="char")
  
  name_netcdf <- paste0('C:/aurore/these/christophe/netcdf_', rep_camp, '.nc')
  nc_var <- nc_create(name_netcdf, list(var.sv,  
                                        var.lon,var.lat,  
                                        var.date.char))
  
  mamaSv <- array(c(mat_38, mat_ent_70, mat_ent_120, mat_ent_200),
                  dim = c(dim(mat_38)[1], dim(mat_38)[2], 4))
  
  ncvar_put(nc_var, var.sv    , mamaSv)
  
  
  ncvar_put(nc_var, var.lon     , lon   )
  ncvar_put(nc_var, var.lat     , lat   )
  ncvar_put(nc_var, var.date.char  , date_char)
  
  nc_close(nc_var)
  
  
}
sapply(list_camp_maracas, do.netcdf)


# var_ncdf <- nc_open("C:/aurore/these/christophe/netcdf_M7C.nc")
# test <- ncvar_get(var_ncdf, 'Sv')
# image(test[,,1])
# dim(test[,,1])

###################################################
#################### ADCP #########################
###################################################

df_bathy <- read.table(file = "D:/maracas/data_maracas/xyz_500m.txt",
                       sep = ",", dec = ".", head = FALSE)
names(df_bathy) <- c("lon", "lat", "depth")

df_bathy <- df_bathy %>% filter(lon %between% c(167.8, 168.6) &
                                  lat %between% c(-23.7, -22.7) &
                                  depth <= 0)

#################### 
load(file = "D:/maracas/data_maracas/df_tot_4freqMAR_90_1000m.Rdata") 


df38 <- df_tot_4freqMAR_90_1000m %>% 
  dplyr::filter(freq == '38')

df38small  <- df38 %>% 
  group_by(lon, lat, camp, moment) %>%
  summarise(Sa = mean(Sa, na.rm = TRUE),
            Sv = mean(Sv, na.rm = TRUE))



plN = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5), 
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5), 
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df38small[df38small$moment %in%
                                c('Night'), ],
             aes(x = lon, y = lat, col = Sv), size = 0.8) +
  facet_grid(moment ~ camp) + 
  scale_color_viridis_c()
  
plD = ggplot() + 
  stat_contour(data = df_bathy, aes(x = lon, y = lat, z = depth), 
                              color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5), 
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5), 
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df38small[df38small$moment %in%
                                c('Day'), ],
             aes(x = lon, y = lat, col = Sv), size = 0.8) +
  facet_grid(moment ~ camp) + 
  scale_color_viridis_c()

ii = ggarrange(plD, plN, nrow = 2)

ggsave(filename ="D:/maracas/R_figures/Sv_moyen_0_800m.tiff",
       width = 4, height = 3, units = "in",dpi = 150,
       scale = 2, plot= ii)

#################### 

df38$vertical_layer <- base::cut(-df38$depth, c(0, 200, 800 ))

df38_2 <- df38 %>% 
  dplyr::filter(moment %in% c("Day", "Night") &
                  lon %between% c(167.8, 168.6) & 
                  lat %between% c(-23.7, -22.7)  ) %>% 
  dplyr::group_by( lon, lat, moment, vertical_layer) %>% 
  dplyr::summarise(sv_mea = mean(sv, na.rm = TRUE),
                   sv_q50 = quantile(sv, probs = 0.5, na.rm = TRUE))
df38_2 <- data.frame(df38_2)
df38_2$Sv_mean <- 10 * log10(df38_2$sv_mea)
df38_2$Sv_q50 <- 10 * log10(df38_2$sv_q50)

df38_2$nasc_mean <- df38_2$sv_mea * 4 * pi * 1852 * 1852 * 20
df38_2$nasc_q50 <- df38_2$sv_q50 * 4 * pi * 1852 * 1852 * 20

plN_1 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df38_2[df38_2$moment == "Night" &
                             df38_2$vertical_layer == "(0,200]", ], 
             aes(x = lon, y = lat, color = nasc_q50), size = 0.9) +
  facet_grid( vertical_layer~ moment) + 
  scale_color_viridis_c(limits = c(0, 200), na.value = "#FDE725FF")
plN_2 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df38_2[df38_2$moment == "Day" &
                             df38_2$vertical_layer == "(0,200]", ], 
             aes(x = lon, y = lat, color = nasc_q50), size = 0.9) +
  facet_grid( vertical_layer~ moment) + 
  scale_color_viridis_c(limits = c(0, 60), na.value = "#FDE725FF")
plN_3 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df38_2[df38_2$moment == "Night" &
                             df38_2$vertical_layer == "(200,800]", ], 
             aes(x = lon, y = lat, color = nasc_q50), size = 0.9) +
  facet_grid( vertical_layer~ moment) + 
  scale_color_viridis_c(limits = c(0, 250), na.value = "#FDE725FF")
plN_4 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_point(data = df38_2[df38_2$moment == "Day" &
                             df38_2$vertical_layer == "(200,800]", ], 
             aes(x = lon, y = lat, color = nasc_q50), size = 0.9) +
  facet_grid( vertical_layer~ moment) + 
  scale_color_viridis_c(limits = c(0, 80), na.value = "#FDE725FF")



ii = ggarrange(plN_2,plN_1,  plN_4,plN_3,  ncol = 2, nrow = 2)

ggsave(filename ="D:/maracas/R_figures/quantile_nasc_couches_verticales.jpg",
       width = 4, height = 3, 
       scale = 3, plot= ii)


#####################################################################
df38$lon2 <- change.res(df38$lon, 0.01)
df38$lat2 <- change.res(df38$lat, 0.01)

df38_2 <- df38 %>% 
  dplyr::filter(moment %in% c("Day", "Night") &
                  lon %between% c(167.8, 168.6) & 
                  lat %between% c(-23.7, -22.7)  ) %>% 
  dplyr::group_by( lon2, lat2, moment, vertical_layer) %>% 
  dplyr::summarise(sv_mea = mean(sv, na.rm = TRUE),
                   sv_q50 = quantile(sv, probs = 0.5, na.rm = TRUE))
df38_2 <- data.frame(df38_2)
df38_2$nasc_mean <- df38_2$sv_mea * 4 * pi * 1852 * 1852 * 20
df38_2$nasc_q50 <- df38_2$sv_q50 * 4 * pi * 1852 * 1852 * 20

plN_1 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Night" &
                             df38_2$vertical_layer == "(0,200]", ], 
             aes(x = lon2, y = lat2, fill = nasc_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 200), na.value = "#FDE725FF")
plN_2 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Day" &
                             df38_2$vertical_layer == "(0,200]", ], 
             aes(x = lon2, y = lat2, fill = nasc_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 40), na.value = "#FDE725FF")
plN_3 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Night" &
                             df38_2$vertical_layer == "(200,800]", ], 
             aes(x = lon2, y = lat2, fill = nasc_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 100), na.value = "#FDE725FF")
plN_4 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Day" &
                             df38_2$vertical_layer == "(200,800]", ], 
             aes(x = lon2, y = lat2, fill = nasc_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 70), na.value = "#FDE725FF")



ii = ggarrange(plN_1, plN_2, plN_3, plN_4, ncol = 2, nrow = 2)

ggsave(filename ="D:/maracas/R_figures/quantile_nasc_couches_verticales_large.jpg",
       width = 4, height = 3, 
       scale = 3, plot= ii)


#####################################################################
df38$lon2 <- change.res(df38$lon, 0.01)
df38$lat2 <- change.res(df38$lat, 0.01)

df38_2 <- df38 %>% 
  dplyr::filter(moment %in% c("Day", "Night") &
                  lon %between% c(167.8, 168.6) & 
                  lat %between% c(-23.7, -22.7)  ) %>% 
  dplyr::group_by( lon2, lat2, moment, vertical_layer) %>% 
  dplyr::summarise(sv_mea = mean(sv, na.rm = TRUE),
                   sv_q50 = quantile(sv, probs = 0.5, na.rm = TRUE))
df38_2 <- data.frame(df38_2)

plN_1 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Night" &
                            df38_2$vertical_layer == "(0,200]", ], 
            aes(x = lon2, y = lat2, fill = sv_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 2.5e-7))
plN_2 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Day" &
                            df38_2$vertical_layer == "(0,200]", ], 
            aes(x = lon2, y = lat2, fill = sv_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 5e-8))
plN_3 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Night" &
                            df38_2$vertical_layer == "(200,800]", ], 
            aes(x = lon2, y = lat2, fill = sv_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 1.5e-7))
plN_4 = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2[df38_2$moment == "Day" &
                            df38_2$vertical_layer == "(200,800]", ], 
            aes(x = lon2, y = lat2, fill = sv_q50)) +
  facet_grid( vertical_layer ~ moment) + 
  scale_fill_viridis_c(limits = c(0, 6e-8))



ii = ggarrange(plN_1, plN_2, plN_3, plN_4, ncol = 2, nrow = 2)

ggsave(filename ="D:/maracas/R_figures/couches_moyennes_large_MEAN.jpg",
       width = 4, height = 3, 
       scale = 3, plot= ii)


###################################################
#################### ADCP #########################
###################################################

df_for_couche <- read.csv2(file = 'D:/maracas/data_maracas/couche_code.csv')

df38$depth= abs(df38$depth)
df_for_couche$depth= abs(df_for_couche$depth)

df38 = merge(df38, df_for_couche)
df38$camp = as.factor(df38$camp)
df38$large_couche = as.factor(df38$large_couche)
levels(df38$large_couche ) = c("Epepelagic",  "Upper meso" , "Lower meso")

aa = ggplot(data = df38[df38$moment %in% c('Day', 'Night'), ], 
       aes(x = moment, y = Sv, fill = camp)) + 
  geom_boxplot() +
  facet_grid(~ large_couche) +
  theme_minimal() + coord_flip() +
  scale_fill_viridis_d(option = 'cividis')


ggsave(filename ="D:/maracas/R_figures/boxplotSV.jpg",
       width = 4, height = 3, units = "in",dpi = 150,
       scale = 2, plot= aa)


###################################################
############## figures differences ################
###################################################
df38$lon2 <- change.res(df38$lon, 0.01)
df38$lat2 <- change.res(df38$lat, 0.01)

df38_2 <- df38 %>% 
  dplyr::filter(moment %in% c("Day", "Night") &
                  lon %between% c(167.8, 168.6) & 
                  lat %between% c(-23.7, -22.7)  ) %>% 
  dplyr::group_by( camp, lon2, lat2, moment, vertical_layer) %>% 
  dplyr::summarise(sv_mea = mean(sv, na.rm = TRUE),
                   sv_q50 = quantile(sv, probs = 0.5, na.rm = TRUE))
df38_2 <- data.frame(df38_2)
df38_2$nasc_mean <- df38_2$sv_mea * 4 * pi * 1852 * 1852 * 20
df38_2$nasc_q50 <- df38_2$sv_q50 * 4 * pi * 1852 * 1852 * 20


plpl = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df38_2,  aes(x = lon2, y = lat2, fill = nasc_mean)) +
  facet_grid(  camp ~  vertical_layer +moment) + 
  scale_fill_viridis_c(limits = c(0, 400))

plpl

df38_2_night_surf <- df38_2 %>% 
  dplyr::filter(moment == 'Night' & vertical_layer == "(0,200]") %>% 
  dplyr::mutate(nasc_surf_night = nasc_mean   ) %>% 
  dplyr::select( camp, lon2,  lat2,nasc_surf_night )

df38_2_day_prof <- df38_2 %>% 
  dplyr::filter(moment == 'Day' & vertical_layer == "(200,800]") %>% 
  dplyr::mutate(nasc_prof_day = nasc_mean   ) %>% 
  dplyr::select( camp, lon2,  lat2,nasc_prof_day )

df_test <- merge(df38_2_night_surf, df38_2_day_prof)
df_test$diff <- df_test$nasc_surf_night  - df_test$nasc_prof_day


plpl = ggplot() + 
  stat_contour(data = df_bathy,
               aes(x = lon, y = lat, z = depth),
               color="grey60", size=0.25) +
  xlab('') + ylab('') + 
  scale_x_continuous(breaks = seq(167.5, 168.6, 0.5),
                     limits = c(167.8, 168.6)) +
  scale_y_continuous(breaks = seq(-23.5, -22.7, 0.5),
                     limits = c(-23.7, -22.7)) +
  coord_equal() + theme_minimal() +
  geom_tile(data = df_test,  aes(x = lon2, y = lat2, fill = diff)) +
  facet_wrap( ~  camp) + 
  scale_fill_distiller(type = "div", limits = c(-250, 250))
plpl

###################################################
#################### solene code #########################
###################################################
load("D:/maracas/data_maracas/env_stack.RData")

# iso_sf <- st_read(dsn = "D:/maracas/data_maracas/bathy/contour.shp")
# iso_sf <- st_crop(iso_sf, raster::extent(167.7, 168.7, -23.8, -22.7))

st_area <- st_read(dsn = "D:/maracas/data_maracas/bathy/maracas7_study_area.shp")
st_area$ELEV <- NULL

# remove parts outside study area
st_area_utm <- st_transform(st_area, crs = 32758)
env_stack$SLOPE <- raster::mask(env_stack$SLOPE, st_area_utm)
env_stack$SUMMIT <- raster::mask(env_stack$DIS_SUMMIT , st_area_utm)

# convert the raster to a star object
summ_star <- st_as_stars(env_stack$DEPTH)




# convert the raster to a star object
slope_star <- st_as_stars(env_stack$SLOPE)

slope_sf <- st_contour(slope_star, contour_lines = FALSE,
                       breaks = c(0.15, max(slope_star$SLOPE, na.rm = T)),
                       na.rm = F)[2, ] # le deuxieme element c'est l'interieur des polygons de ce contour, le premier element c'était l'exterieur.
slope_sf$id <- "slope" ; slope_sf$SLOPE <- NULL # rename the polygons

slope_sf <- st_cast(slope_sf, to = "POLYGON")
slope_sf <- slope_sf[order(as.numeric(st_area(slope_sf)), decreasing = T)[1:3], ]
slope_sf <- st_transform(slope_sf, crs = 4326) # project to same crs
slope_sf$id <- "slope" ; slope_sf$Min <- NULL ; slope_sf$Max <- NULL 

plot(slope_sf, col = "#5CC5EF")

summit_sf <- st_read("D:/maracas/data_maracas/contour_summits_EK60analysis.shp")
summit_sf <- st_cast(summit_sf, to = "POLYGON")
summit_sf$id <- "summit" ; summit_sf$ELEV <- NULL # rename the polygons

# slope_sf <- st_read("D:/maracas/data_maracas/contour_slopes_EK60analysis.shp")
# slope_sf <- st_cast(slope_sf, to = "POLYGON")
# slope_sf <- st_make_valid(slope_sf)

st_area$id <- "open_water"

hab_sf <- rbind(st_area, slope_sf, summit_sf)
plot(hab_sf, pal = canva_pal(palette = "Pool party"))

df38$depth_bin <- round(df38$depth/10, 0)*10
df38$date_char <- as.character(df38$date)
df38_bin <- data.frame(df38) %>%
  group_by(lon, lat, date_char) %>%
  dplyr::summarise(Sv = 10 * log10(mean(sv, na.rm = TRUE))) %>% 
  dplyr::select(-Sv)

# convert to spatial object and remove parts outside the study area
df38_sf <- st_as_sf(df38_bin, coords = c("lon", "lat"), crs = 4326)
i <- st_intersects(df38_sf, st_area)
df38_sf <- df38_sf[apply(i, 1, any), ]
df38_sf$habitat <- "open_water"

sf_use_s2(FALSE)
i <- st_intersects(df38_sf, slope_sf)
df38_sf[apply(i, 1, any), ]$habitat <- "slope"

sf_use_s2(TRUE)
i <- st_intersects(df38_sf, summit_sf)
df38_sf[apply(i, 1, any), ]$habitat <- "summit"


df_try <- do.call(rbind, st_geometry(df38_sf)) %>% 
  as_tibble() %>% setNames(c("lon","lat"))
df_try <- data.frame(df_try)
df_try$habitat <- df38_sf$habitat
df_try$date_char <- df38_sf$date_char


df38_bis <- merge(df38, df_try, by = c("lon", "lat", "date_char"))
df38_bis <- df38_bis %>% dplyr::filter(moment %in% c('Day', 'Night'))
df38_bis$nasc <- df38_bis$sv * 4 * pi * 1852 * 1852 * 20
df38_bis <- df38_bis %>% dplyr::filter(nasc <= 400)

g3_n <- ggplot(df38_bis[df38_bis$moment == 'Night', ],
               aes(x = habitat, y = nasc, fill = habitat)) +
  geom_point(aes(color = habitat),
             position = position_jitter(w = .1),
             alpha = 0.1, size = 0.01) +
  geom_boxplot(outlier.alpha = 0, col = "grey50", width = 0.2) +
  ylim(0, 150) + coord_flip() +
  geom_flat_violin(position = position_nudge(x = .1),
                   trim = TRUE, 
                   alpha = 0.7, 
                   scale = "width",
                   color = NA) +
  facet_grid( vertical_layer ~ moment, scales = 'free') + 
  scale_fill_canva(palette = "Pool party", name = "Habitat") + 
  scale_color_canva(palette = "Pool party", name = "Habitat") + 
  theme_minimal()  + theme(axis.text.y = element_blank(),
                           axis.title.y  = element_blank(),
                           panel.background = element_rect(color = 'black',
                                                           fill ='white'))

g3_d <- ggplot(df38_bis[df38_bis$moment == 'Day', ],
               aes(y = nasc, fill = habitat, x = habitat)) +
  geom_point(aes(color = habitat),
             position = position_jitter(w = .1),
             alpha = 0.1, size = 0.01) +
  geom_boxplot(outlier.alpha = 0, col = "grey50", width = 0.2) +
  ylim(0, 50) + coord_flip() +
  geom_flat_violin(position = position_nudge(x = .1),
                   trim = TRUE, 
                   alpha = 0.7, 
                   scale = "width",
                   color = NA) +
  facet_grid(vertical_layer ~ moment, scales = 'free') + 
  scale_fill_canva(palette = "Pool party", name = "Habitat") + 
  scale_color_canva(palette = "Pool party", name = "Habitat") + 
  theme_minimal()  + theme(axis.text.y  = element_blank(),
                           axis.title.y  = element_blank(),
                           panel.background = element_rect(color = 'black',
                                                           fill ='white'))

a = ggarrange(g3_d,g3_n, common.legend = TRUE)
a
ggsave(filename ="D:/maracas/R_figures/boxplot_seamounts.jpg",
       width = 3, height = 1.8, 
       scale = 3, plot= a)



m <- lm(nasc ~ camp  + moment:habitat:vertical_layer, 
        df38_bis)
summary(m)



df38_bis$lon2 <- change.res(df38_bis$lon, 0.01)
df38_bis$lat2 <- change.res(df38_bis$lat, 0.01)

df38_bis2 <- df38_bis %>% 
  dplyr::group_by( lon2, lat2, moment, depth_bin, camp, habitat,vertical_layer ) %>% 
  dplyr::summarise(sv_mea = mean(sv, na.rm = TRUE),
                   sv_q50 = quantile(sv, probs = 0.5, na.rm = TRUE))
df38_bis2$nasc <- df38_bis2$sv_mea * 4 * pi * 1852 * 1852 * 20


m <- glm(nasc ~ camp  + moment:vertical_layer*habitat, 
         df38_bis2,
        family = 'Gamma')
summary(m)



#####################################################################


# g1 <- ggplot(df38_bis, aes(x = 1, y = Sv, fill = habitat)) +
#   geom_boxplot(outlier.alpha = 0, col = "lightgrey") +
#   facet_grid( moment ~ vertical_layer, scales = 'free') +ylim(-100, -60) +
#   scale_fill_canva(palette = "Pool party", name = "Habitat") + 
#   theme_minimal() 




df38_bis <- merge(df38, df_try, by = c("lon", "lat", "date_char"))
df38_bis <- df38_bis %>% dplyr::filter(moment %in% c('Day', 'Night'))

g3_n <- ggplot(df38_bis[df38_bis$moment == 'Night', ],
               aes(x = habitat, y = Sv, fill = habitat)) +
  geom_point(aes(color = habitat),
             position = position_jitter(w = .1),
             alpha = 0.1, size = 0.01) +
  geom_boxplot(outlier.alpha = 0, col = "grey50", width = 0.2) +
  ylim(-100, -50) +
  coord_flip() +
  geom_flat_violin(position = position_nudge(x = .1),
                   trim = TRUE, 
                   alpha = 0.7, 
                   scale = "width",
                   color = NA) +
  facet_grid( vertical_layer ~ moment, scales = 'free') + 
  scale_fill_canva(palette = "Pool party", name = "Habitat") + 
  scale_color_canva(palette = "Pool party", name = "Habitat") + 
  theme_minimal()  + theme(axis.text.y = element_blank(),
                           axis.title.y  = element_blank(),
                           panel.background = element_rect(color = 'black',
                                                           fill ='white'))

g3_d <- ggplot(df38_bis[df38_bis$moment == 'Day', ],
               aes(y = Sv, fill = habitat, x = habitat)) +
  geom_point(aes(color = habitat),
             position = position_jitter(w = .1),
             alpha = 0.1, size = 0.01) +
  geom_boxplot(outlier.alpha = 0, col = "grey50", width = 0.2) +
  ylim(-100, -50) +
  coord_flip() +
  geom_flat_violin(position = position_nudge(x = .1),
                   trim = TRUE, 
                   alpha = 0.7, 
                   scale = "width",
                   color = NA) +
  facet_grid(vertical_layer ~ moment, scales = 'free') + 
  scale_fill_canva(palette = "Pool party", name = "Habitat") + 
  scale_color_canva(palette = "Pool party", name = "Habitat") + 
  theme_minimal()  + theme(axis.text.y  = element_blank(),
                           axis.title.y  = element_blank(),
                           panel.background = element_rect(color = 'black',
                                                           fill ='white'))

a = ggarrange(g3_d,g3_n, common.legend = TRUE)

ggsave(filename ="D:/maracas/R_figures/boxplot_seamounts_Sv.jpg",
       width = 3, height = 1.8, 
       scale = 3, plot= a)


