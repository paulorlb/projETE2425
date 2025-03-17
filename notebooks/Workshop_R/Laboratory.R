##########################################################
##### A gentle introduction to spatial data analysis #####
#####       Laboratory on spatial analysis           #####
#####           Universidade de Aveiro               #####
#####               April, 9th 2024                  #####
##########################################################

##### Libraries
# Eurostat
library(eurostat)
# Data manipulation and graphical libraries
library(tidyverse)
library(lubridate)
library(ggpubr)
library(ggplot2)
# Spatial analysis libraries
library(sf)
library(spdep)


##### Working directory
setwd("H:/Il mio Drive/Teaching/Erasmus_UAA_2024")

##### Load spatial data
load("Panel_Eurostat_UAA2024.RData")
str(NUTS2021)

##### Convert data.frame to sf object
NUTS2021_sf <- st_as_sf(NUTS2021)

##### Extract geometries for continental Portugal
View(euro_sf)
PT_sf <- euro_sf %>%
  filter(CNTR_CODE %in% c("PT"),
         Longitude > -15)


##### Spatial contiguity matrix
#   The matrix usually denoted by W is a N by N positive and symmetric matrix which denotes fore each observation (row)
#   those locations (columns) that belong to its neighborhood set as nonzero elements (Anselin and Bera (1998),
#   Arbia (2014))
# https://www.paulamoraga.com/tutorial-neighborhoodmatrices/ 
# https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html
# https://rpubs.com/erikaaldisa/spatialweights
# https://spatialanalysis.github.io/lab_tutorials/Contiguity_Spatial_Weights.html#r-packages-used

# Convert sf to sp data
PT_sp <- sf::as_Spatial(PT_sf)

# Extract centroids' coordinates
coords <- st_coordinates(st_centroid(PT_sf))
colnames(coords) <- c("Longitude","Latitude")

# Creating Queen's Contiguity based neighbors
contnb_q <- spdep::poly2nb(PT_sp, queen = TRUE, row.names = PT_sp$NUTS_NAME)
summary(contnb_q)
str(contnb_q)
contmat_q <- nb2mat(contnb_q, zero.policy = TRUE, style = "W")
colnames(contmat_q) <- PT_sp$NUTS_NAME
contmat_q

# Creating Rook's Contiguity based neighbors
contnb_r <- spdep::poly2nb(PT_sp, queen = FALSE, row.names = PT_sp$NUTS_NAME)
summary(contnb_r)
str(contnb_r)
contmat_r <- nb2mat(contnb_r, zero.policy = TRUE)

# List of neighbors
contnb_q_listw <- spdep::nb2listw(neighbours = contnb_q, style = "W", zero.policy = TRUE)
contnb_r_listw <- spdep::nb2listw(neighbours = contnb_r, style = "W", zero.policy = TRUE)

# Mapping
png(filename = "ContMat_PT.png",width = 1800, height = 1200, res = 150)
par(mfrow = c(1,2))
plot(x = st_geometry(PT_sf), border = 'black', col="white", lwd=2, main = "Contiguity criterion (Queen) for Portugal")
plot(contnb_q, coords = coords, pch = 19, cex = 0.8, add = TRUE, col = 'red')
plot(x = st_geometry(PT_sf), border = 'black', col="white", lwd=2, main = "Contiguity criterion (Rook) for Portugal")
plot(contnb_r, coords = coords, pch = 19, cex = 0.8, add = TRUE, col = 'darkgreen')
par(mfrow = c(1,1))
dev.off()


##### Fixed distance weight matrix
# Setting the longlat=TRUE option determines the distance measurement in kilometres
# (according to the Great Circle distnace), while longlat = FALSE or omitting results in expressing
# the distance in degrees.
# Centroid within 60 kilometers
wm_d60 <- dnearneigh(coords, d1 = 0, d2 = 60, longlat = TRUE)
# Centroid within 100 kilometers
wm_d100 <- dnearneigh(coords, d1 = 0, d2 = 100, longlat = TRUE)
# Centroid within 200 kilometers
wm_d200 <- dnearneigh(coords, d1 = 0, d2 = 200, longlat = TRUE)
# Maps
png(filename = "FixDist_PT.png",width = 1800, height = 1200, res = 150)
par(mfrow = c(1,3))
plot(st_geometry(PT_sf), border = "darkgrey")
plot(wm_d60, coords, add = TRUE, pch = 19, cex = 1)
title(main = "Neighbors within 60 km", cex.main = 1.1)
plot(st_geometry(PT_sf), border = "darkgrey")
plot(wm_d100, coords, col = "red", add = TRUE, pch = 19, cex = 1)
title(main = "Neighbors within 100 km", cex.main = 1.1)
plot(st_geometry(PT_sf), border = "darkgrey")
plot(wm_d200, coords, col = "blue", add = TRUE, pch = 19, cex = 1)
title(main = "Neighbors within 200 km", cex.main = 1.1)
par(mfrow = c(1,1))
dev.off()


##### Inverse-distance weighting
distnb <- nbdists(nb = contnb_q, coords = coords)
# IDW distance
IDW_dist <- lapply(distnb, function(x) 1/x)
IDW_dist_listw <- nb2listw(contnb_q, glist = IDW_dist)
IDW_dist_mat <- listw2mat(IDW_dist_listw)
# Exponential distance
exp_dist <- lapply(distnb, function(x) exp(-1.5*x))
exp_dist_listw <- nb2listw(contnb_q, glist = exp_dist)
exp_dist_mat <- listw2mat(exp_dist_listw)

# Heatmap of weights
png(filename = "IDW_heatmap.png",width = 1800, height = 1200, res = 150)
lattice::levelplot(t(IDW_dist_mat))
dev.off()

# Mapping weights
exp_dist_mat
rownames(exp_dist_mat)
PT_sf <- PT_sf %>%
  mutate("Aveiro - Exponential" = exp_dist_mat[22,],
         "Aveiro - IDW" = IDW_dist_mat[22,],
         "Coimbra - Exponential" = exp_dist_mat[23,],
         "Coimbra - IDW" = IDW_dist_mat[23,])

p_dist <- PT_sf %>%
  dplyr::select("Aveiro - Exponential","Aveiro - IDW","Coimbra - Exponential","Coimbra - IDW") %>%
  pivot_longer(names_to = "Distance",values_to = "weight",
               cols = c("Aveiro - Exponential","Aveiro - IDW","Coimbra - Exponential","Coimbra - IDW")) %>%
  st_as_sf() %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = weight)) + 
  facet_wrap(~ Distance, ncol = 4) + 
  scale_fill_binned(type = "viridis",
                    breaks = seq(from=0,to=0.4,by=0.05), 
                    name = "Weight" ) +
  theme_bw() + 
  theme(plot.title = element_text(size = 20)) + 
  labs(title = "Distance-based weights for Aveiro and Coimbra",
       subtitle =  latex2exp::TeX("Exponential weighting VS IDW weighting"),
       x = "", y = "")
p_dist
ggpubr::ggexport(plotlist = list(p_dist), filename = "Distance_Aveiro_Coimbra.png", width = 1800, height = 1200, res = 150)


##### k-Nearest Neighbors distance weight matrix
# 1-NN (the nearest neighbor)
NN1_sf <- spdep::knearneigh(coords, k=1)
NN1_nb <- spdep::knn2nb(NN1_sf)
NN1_cnt <- spdep::nb2mat(NN1_nb, style = "B")
# 2-NN (the nearest neighbor)
NN2_sf <- spdep::knearneigh(coords, k=2)
NN2_nb <- spdep::knn2nb(NN2_sf)
NN2_cnt <- spdep::nb2mat(NN2_nb, style = "B")
# 5-NN (the nearest neighbor)
NN5_sf <- spdep::knearneigh(coords, k=5)
NN5_nb <- spdep::knn2nb(NN5_sf)
NN5_cnt <- spdep::nb2mat(NN5_nb, style = "B")
# Maps
png(filename = "kNN_PT.png",width = 1800, height = 1200, res = 150)
par(mfrow = c(1,3))
plot(st_geometry(PT_sf), border = "darkgrey")
plot(NN1_nb, coords, col = "black", add = TRUE, pch = 19, cex = 1)
title(main = "1-NN Neighbors", cex.main = 1.1)
plot(st_geometry(PT_sf), border = "darkgrey")
plot(NN2_nb, coords, col = "red", add = TRUE, pch = 19, cex = 1)
title(main = "2-NN Neighbors", cex.main = 1.1)
plot(st_geometry(PT_sf), border = "darkgrey")
plot(NN5_nb, coords, col = "blue", add = TRUE, pch = 19, cex = 1)
title(main = "5-NN Neighbors", cex.main = 1.1)
par(mfrow = c(1,1))
dev.off()



######################################################
########## Spatial autocorrelation analysis ##########
######################################################

##### Extract data for Portugal
NUTS2021_PT <- NUTS2021_sf %>%
  filter(NUTS0_Code %in% c("PT"))

##### Compute per capita GDP 
NUTS2021_PT <- NUTS2021_PT %>%
  mutate(GDPpercapita = GDP_Mln_PPS/(Pop_T_Tot/1000000))

##### Dynamic plot of per capita GDP
p2 <- NUTS2021_PT %>%
  select(Year,NUTS2_Code,GDPpercapita,geometry) %>%
  st_as_sf() %>%
  ggplot() + 
  geom_sf(mapping = aes(fill = GDPpercapita)) + 
  coord_sf(xlim = c(-10, -6), ylim = c(36.5, 43), expand = FALSE) +
  facet_wrap(~Year, ncol = 8) + 
  scale_fill_binned(type = "viridis",
                    breaks = round(unname(quantile(NUTS2021_PT$GDPpercapita, probs = seq(0.1,1,by=0.07))),2), 
                    name = "Per capita €" ) +
  theme_bw() + 
  theme(plot.title = element_text(size = 20)) + 
  labs(title = "Example with provincial data for Portugal",
       subtitle =  latex2exp::TeX("Per capita GDP (€ PPS) in Portugal (NUTS-3 classification)"),
       x = "", y = "")
p2 <- annotate_figure(p = p2, bottom = text_grob("Data source: Eurostat", face = "italic", size = 10))
p2
ggpubr::ggexport(plotlist = list(p2), filename = "GDP_Portugal.png", width = 1800, height = 1200, res = 150)


##### Extract a single-year data (continental Portugal)
DataPT_2021 <- NUTS2021_PT %>%
  filter(Year == 2021, Longitude > -15)

##### Compute the Global Moran's statistic (with Queen's contiguity neighbors)
moran.test(DataPT_2021$GDPpercapita, contnb_q_listw)
moran.plot(DataPT_2021$GDPpercapita, contnb_q_listw)
png(filename = "MoranScatter_PT.png",width = 1800, height = 1200, res = 150)
moran.plot(DataPT_2021$GDPpercapita, contnb_q_listw,pch = 19,
           xlab = "GDP per capita (€)",
           ylab = "Lag-1 GDP per capita (€)",
           main = "Moran's scatterplot for the Portuguese provinces (2021)")
dev.off()

##### Compute the Moran's scatterplot manually
x <- DataPT_2021$GDPpercapita                               # creating a variable for analysis
zx <- scale(x)                                              # variable standardisation
wzx <- lag.listw(contnb_q_listw, zx)                        # spatial lag of x
mean(zx)                                                    # mean control = 0
sd(zx)                                                      # standard deviation control = 1
morlm <- lm(wzx~zx)                                         # linear regression model
slope <- morlm$coefficients[2]                              # directional coefficient
intercept <- morlm$coefficients[1]                          # constant term
par(pty="s")                                                # square chart window
plot(zx, wzx, xlab="zx",ylab="spatial lag zx", pch="*")     # Scatterplot
abline(intercept, slope)                                    # regression line
abline(h=0, lty=2)                                          # horizontal line at y = 0
abline(v=0, lty=2)                                          # vertical line at x = 0


##### Map of belonging to the quarters of the Moran scatterplot
# Assigning regions to quadrants
cond1 <- ifelse(zx>=0 & wzx>=0, 1,0)  # I quarter
cond2 <- ifelse(zx>=0 & wzx<0, 2,0)   # II quarter
cond3 <- ifelse(zx<0 & wzx<0, 3,0)    # III quarter
cond4 <- ifelse(zx<0 & wzx>=0, 4,0)   # IV quarter
cond.all <- as.data.frame(cond1+cond2+cond3+cond4)
# Map
brks<-c(1,2,3,4)
cols<-c("grey25", "grey60", "grey85", "grey45")

png(filename = "MoranMap_PT.png",width = 1800, height = 1200, res = 150)
plot(st_geometry(PT_sf), col=cols[findInterval(cond.all$V1, brks)])
legend("bottomleft",
       legend=c("I quarter - HH - high surrounded by high",
                "II quarter - LH - low surrounded by high",
                "III quarter - LL - low surrounded by low",
                "IV quarter - HL - high surrounded by low"),
       fill=cols, bty="n", cex=0.8)
title(main="Moran's map for the Portoguese provinces (2021)")
dev.off()


##### Compute the Global Moran's statistic (with 5-NN neighbors)
moran.test(DataPT_2021$GDPpercapita, spdep::mat2listw(NN5_cnt))
moran.plot(DataPT_2021$GDPpercapita, spdep::mat2listw(NN5_cnt,row.names = PT_sp$NUTS_NAME))
