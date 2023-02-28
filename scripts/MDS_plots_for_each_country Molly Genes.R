# Molly Dupin
# New Mexico State University
# Amazona auropalliata acoustic analysis
# MDS plots

setwd("~/Desktop/R/MS_thesis/Files_for_final_analysis/Vocal_data")
getwd()

X <- c("warbleR", "Rraven", "parallel", "devtools", "dplyr", "ggfortify", "ggplot2", "pbapply")
invisible(lapply(X, library, character.only = TRUE))

sels <- imp_raven(path = "~/Desktop/R/MS_thesis/Files_for_final_analysis", name.from.file = TRUE, ext.case = "lower", rm_dup = TRUE, freq.cols = FALSE) 
str(sels)
class(sels)
View(sels)

write.csv(sels, file = "Sels.csv", row.names = FALSE)
#############################################################
# Start here when running. Read in .csv
Sys.getlocale()
Sys.setlocale(locale="C")

sels <- read.csv("ReducedVocalDataYNA.csv", header = TRUE)
View(sels)
# View(sels)
# run SPCC on all calls
# first check number of cores
parallel::detectCores()
# 4

# set cores to be total - 2
cores <- parallel::detectCores() - 2
cores # 2

# run the spectrographic cross-correlation on calls
# grep out specific group (region, country, etc)
xc_mat <- xcorr(X = sels[grep("Mexico", sels$Country), ], wl = 512, ovlp = 85, bp = c(0.5, 2.5), parallel = cores, pb = TRUE, cor.mat = TRUE, na.rm = TRUE)
str(xc_mat)
View(xc_mat)

# let's check these calls out using multidimensional scaling (MDS) to visualize SPCC acoustic space
?stats::cmdscale

xc_dist <- stats::as.dist(1-xc_mat, diag = TRUE, upper = FALSE)
class(xc_dist)

# use only 2 dimensions b/c just visualizing
mds_res <- stats::cmdscale(d = xc_dist, k = 2, add = TRUE)

# create a data frame for ggplotting, using the new coordinates determined by MDS
# first extract site names from xc_mat dimension names
nms <- dimnames(xc_mat)[[1]]
site <- sapply(1:length(nms), function(x){
  strsplit(nms[x], split = "_")[[1]][5]
})

temp <- filter(sels, grepl("Mexico", sels$Country))
View(temp)
Call_Type <- temp$Call_Type
View(Call_Type)

df <- data.frame(X = mds_res$points[, 1], Y = mds_res$points[, 2], Call_Type = Call_Type)
str(df)

hulls <- plyr::ddply(df, "Call_Type", function(x){x[chull(x$X, x$Y), ]})

theme_AcousticSpace <- function(){
  theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        panel.grid.major = element_line(size = 3, colour = "white"),
        panel.grid.minor = element_line(size = 0.75, colour = "white"),
        axis.line = element_line(size = 0.45, colour = "black"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))
}

# make a color palette object with different colors for # sites
cbPalette <- c("navajowhite4", "yellow", "darkorange", "orangered4")

cbPalette2 <- c("cornflowerblue", "coral4", "chartreuse1", "darkorchid4", "coral", "green4", "cadetblue1", "darkslateblue", "black", "darkolivegreen", "darkmagenta", "azure4", "firebrick2", "blue", "deeppink", "grey50", "hotpink4", "khaki1", "springgreen1", "cyan", "plum1", "steelblue", "wheat4")


jpeg("Mexico_SPCC_MDS.jpeg", units = NULL, width = 10, height = 10, res = 300)
ggplot(data = df, 
       aes(x = X, y = Y, color = Call_Type, fill = Call_Type)) +
  geom_point(size = 3) + #shape = Call_Type
  geom_polygon(data = hulls, aes(x = X, y = Y, fill = Call_Type, color = Call_Type), alpha = 0.2, size = 0.2) +
  # ggtitle("YNA vocal data clustered by visually classified call types from Costa Rica, 2016") +
  scale_color_manual(values = cbPalette) + 
  scale_fill_manual(values = cbPalette) +
  xlab("MDS Dimension 1") + 
  ylab("MDS Dimension 2") +
  theme_AcousticSpace() +
  theme(legend.position = "top") 
dev.off()

# kernel

# load necessary packages
X <- c("warbleR", "ggplot2", "pbapply", "parallel", "data.table", "tidyverse", "mclust", "aricode", "ks", "overlap", "scales", "viridis", "egg", "ggplotify", "grid", "effsize", "dplyr")
invisible(lapply(X, library, character.only = TRUE))

# create path object
path <- "~/Desktop/R/MS_thesis/Files_for_final_analysis/Vocal_data"
# selection table for analysis must be in .rds format
sels <- read.csv("ReducedVocalDataYNA.csv")
# create extended selection table
# ccs_est <- selection_table(sels, max.dur = TRUE, path = path, extended = TRUE,  # confirm.extended = FALSE, parallel = 2)
# glimpse(ccs_est)
# saveRDS(sels, "VocalDataYNA_EST.RDS")
ccs_est <- readRDS(file.path(path, "VocalDataYNA_EST.RDS"))
glimpse(ccs_est)
# Duplicated selection
# wh <- which(duplicated(paste(ccs_est$sound.files, ccs_est$selec, sep = "_")))
# as.character(ccs_est$sound.files[wh])
# # Remove duplicated sel-- "2016_07_02_PM_PENA_TW_B6_P2.wav"
# ccs_est <- ccs_est[-wh,]
# ccs_est <- ccs_est %>% 
#   dplyr::select(-c(X))
# names(ccs_est)
# re-save as RDS/duplicate manually removed from original csv and EST remade

sels <- read.csv("ReducedVocalDataYNA.csv", header = TRUE)
View(sels)

check_res <- check_sels(ccs_est)
# make sure all selections are good- do not proceed otherwise

# check available cores
parallel::detectCores()
# set to use 2 less cores than available
cores <- parallel::detectCores() - 2

# Read in the xcorr matrix- no need to do multiple for each region because we will incorporate that into the code further down
xc_mat <- readRDS(file.path(path, "CombinedRegionsXcorr.RDS"))
str(xc_mat) # check structure

# Perform MDS on SPCC matrices 
# Convert SPCC similarity matrix to a distance matrix and dist object
xc_mat_dist <- stats::as.dist(1 - xc_mat, diag = TRUE, upper = TRUE)
str(xc_mat_dist)

# Perform MDS to get 2 dimensions 
# k may need to be changed depending on visible patterns in low dimensional acoustic space
# circular/half circle patterns may be due to non-optimized dimensions. Increase to 10/15. And then just take the first 2 of those to make the points look spread out
# the higher k is set, the lower maxit should be set
iso <- invisible(MASS::isoMDS(xc_mat_dist, k = 4, maxit = 1000, trace = FALSE))
str(iso)
saveRDS(iso, "YNA_CombinedMDS_iso.RDS")
iso <- readRDS("YNA_CombinedMDS_iso.RDS") # read in iso matrices here if you've already made them

# now combine the call metadata with the new mds coordinates for the acoustic space plot
mds_df <- ccs_est %>%
  dplyr::select(sound.files, Year, Country, Region, site, Call_Type, site_year) %>% 
  # incorporate the mds into this new data frame
  dplyr::mutate(
    X = iso$points[,1], Y = iso$points[,2])
glimpse(mds_df)

# This can be saved as an RDS or csv when finished
write.csv(mds_df, file.path(path, "MDS_dataframe_YNA.csv"), row.names = FALSE)
mds_df <- read.csv("MDS_dataframe_YNA.csv")

# read in the .csv from above if you don't wanna do mds again
regions <- mds_df %>%
  pull(Region) %>%
  unique()

country <- mds_df %>%
  pull(Country) %>%
  unique()

# Below are colors grace used, adjust to needs, colors by call TYPE
cols <- c("springgreen3", "gold4", "slateblue1", "darkorchid4")

# "red", "blue", "gold", "purple4", "orange", "green4", "lightpink", "deepskyblue", "magenta", "lightgoldenrod", "grey", "brown", "black", "red4", "orchid4", "tan1", "hotpink4", "indianred", "orchid", "lightseagreen", "lightseablue", "darkolivegreen1", "rosybrown4", "darkolivegreen4", "tan4", "honeydew4"

x <- (2) # change number based on level you want to print 1,2,3,4,5 for all
y <- 1
mds_df_tmp <- mds_df %>%
  filter(Country == country[x]) %>%
  droplevels() 
View(mds_df_tmp)
ggplot(mds_df_tmp, aes(x = X, y = Y, color = Call_Type)) +
  stat_density2d(bins=20) +
  scale_color_manual(values = cols) +
  xlab("MDS Dim 1") +
  ylab("MDS Dim 2") +
  theme(panel.background = element_rect(fill = "white", color = "white"), 
        plot.background = element_rect("white", color = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        axis.line = element_line(size = 0.8, colour = "black"),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        axis.text = element_text(size = 15),
        legend.position = "right") 
dev.off()  
jpeg("Combined_Kernel.jpeg", units = NULL, width = 10, height = 11, res = 300)








