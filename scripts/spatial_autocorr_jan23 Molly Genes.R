# Yellow-Naped Amazon Spatial Autocorrelation
# Molly Dupin 
# Updated January 2023

# set and check working directory
setwd("~/Desktop/R/MS_thesis/Files_for_final_analysis/Vocal_data")
getwd()

# install.packages("sf")
# install.packages("dplyr")
# install.packages("rgdal")
# install.packages("spData")
# install.packages("rgeos")
# install.packages("vegan")
# install.packages("sp")
# install.packages("raster")

# Run prior to loading rgdal to ignore Proj4 warnings
# options("rgdal_show_exportToProj4_warnings"="none")

# Load packages
X <- c("sp", "rgdal", "rgeos", "dplyr", "raster", "dplyr", "spData", "vegan")
invisible(lapply(X, library, character.only = TRUE))

# Read in and convert the dataframe to Spatial Points object.
sels <- read.csv("ReducedVocalDataYNA.csv")
glimpse(sels)

sels$Latitude <- as.factor(sels$Latitude)
sels$Longitude <- as.factor(sels$Longitude)

sels$Latitude <- as.numeric(sels$Latitude)
sels$Longitude <- as.numeric(sels$Longitude)

mat <- data.frame(lon = sels$Longitude, lat = sels$Latitude)
str(mat)

sp_pts <- SpatialPoints(mat, proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs +init:epsg31972"))

# Reproject using EPSG codes
# EPSG codes are specific to regions across the world. Codes can be accessed at the following site: http://epsg.io or by searching for other EPSG repositories online.
epsg <- rgdal::make_EPSG()
epsg[grep("^31972$", epsg$code), ]

sp_pts <- sp::spTransform(sp_pts, CRSobj = CRS(epsg$prj4[grep("^31972$", epsg$code)]))
bbox(sp_pts)
proj4string(sp_pts)

# Calculate the pairwise distances
geo_dists <- rgeos::gDistance(sp_pts, byid = TRUE)
str(geo_dists)

saveRDS(geo_dists, "pairwise_distances.RDS")
geo_dists <- readRDS("pairwise_distances.RDS")
View("geo_dists")

# Projection is in meters, so we convert to km
geo_dists <- round(geo_dists/1000)

# When we plot our data, we want the points to be connected in one continuous line. In order to achieve that, we identified "breaks" in distances using the following code. Break intervals are measured in km, and once a break limit was chosen we checked the number of existing 0's by creating a table of the classes. If a 0 shows in one of the classes, that class has no data and the break number should be altered. Although not essential, minimizing zeros maintains continuity in the plotted line for the graphic.
geo_vect <- geo_dists[lower.tri(geo_dists)]
geo_vect <- geo_vect[geo_vect != 1 & !is.na(geo_vect)]
range(geo_vect)

geo_classes <- cut(geo_vect, breaks = seq(0, round(max(geo_vect)), 100))
table(geo_classes)

classes <- seq(0, max(geo_dists), 100)
length(classes)

# Run a mantel correlation on the data
xc <- readRDS("SPCC2016_2018_2019.RDS")
correl_SPCC_100km <- vegan::mantel.correlog(D.eco = xc, D.geo = geo_dists, break.pts = classes, cutoff = FALSE, r.type = "pearson", nperm = 999, mult = "holm", progressive = TRUE)
str(correl_SPCC_100km)

# Save the correlation as an .RDS file so you don't have to run it multiple times in the future
saveRDS(correl_SPCC_100km, "correl_SPCC_100km.RDS")

# Read in .RDS file if you've already saved one
correl <- readRDS("correl_SPCC_100km.RDS")
str(correl)

# Create a dataframe from the mantel correlation matrix, remove any existing NA's in the p-values. Create a vector of significance values for the plot and match with  appropriate p-values in the dataframe.
df1 <- as.data.frame(correl$mantel.res)
df1 <- df1[!is.na(df1$`Pr(corrected)`), ]

sig <- ifelse(df1$`Pr(corrected)` <= 0.05, "sig", "non-sig")

p_val <- df1$`Pr(corrected)`

distance <- df1$class.index

corr <- df1$Mantel.cor

sig_cat <- sig
sig_cat[sig_cat == "sig" & corr < 0] <- "Lower"
sig_cat[sig_cat == "sig" & corr > 0] <- "Higher"
sig_cat[sig_cat == "non-sig"] <- "Neutral"
unique(sig_cat)

corr_df <- data.frame(distance = distance, corr = corr, p_val = p_val, sig = sig, sig_cat = sig_cat)
corr_df$sig <- factor(corr_df$sig, levels = c("sig", "non-sig"))
corr_df$sig_cat <- factor(corr_df$sig_cat, levels = c("Lower", "Neutral", "Higher"))

# Adjust the x- and y-axis scales using ggplot: c(0, 15)
# install.packages("scales")
# install.packages("ggplot2")
scales_x <- list(Central_America = scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, 5), labels = seq(0, 15, 5)))

scales_y <- list(SPCC = scale_y_continuous(limits = c(-0.5, 0.5)))

str(corr_df)

# Plot the spatial autocorrelogram using ggplot
# jpeg("Spatial_auto100km.tiff", units="in", width=10, height=10, res=300) this never works
ggplot(corr_df, aes(x = distance, y = corr)) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.5) +
  geom_point(aes(color = sig_cat, fill = sig_cat), size = 4, shape = 21) +
  scale_color_manual(values = c(alpha("darkblue", 0.6), gray.colors(12)[6], "red")) +
  scale_fill_manual(values = c(alpha("darkblue", 0.6), gray.colors(12)[6], "red")) +
  geom_line(colour = "black", size = 0.25) +
  theme_bw() +
  theme(panel.background = element_rect(color = "white", fill = "white"),
        panel.grid.major.x = element_line(size = 0.1, color = "black"), 
        axis.line = element_line(color = "black", size = 0.35), 
        axis.text.x = element_text(size = 15, angle = 0, face = "bold"),
        axis.text.y = element_text(size = 15, angle = 0, face = "bold"),
        axis.title = element_text(size = 20), 
        axis.ticks = element_line(size = 0.15), 
        legend.text = element_text(hjust = 0, size = 18), 
        legend.title = element_text(hjust = 0, size = 18), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        legend.position = "top", 
        strip.text = element_text(size = 10)) +
  xlab("Pairwise Geographic Distance (km)") + 
  ylab("Mantel Spatial Correlation of Acoustic Similarity") + 
  guides(fill = guide_legend(title = "Spatial Correlation Values"), color = guide_legend(title = "Spatial Correlation Values"))

dev.off()
