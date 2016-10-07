############################################
###### Open Source Software Workshop #######
################## DEMO ####################
############################################

############################################
########### BEHIND THE SCENES ##############
############################################

#PACKAGES
library(spdep)
library(strptimer)
library(reshape2)
library(dplyr)
library(plyr)
library(ggmap)
library(ggplot2)
library(readr)
library(mapproj)
library(sp)
library(spatstat)
library(ape)

#NECESSARY FUNCTIONS
#Source: http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
##############################################################
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
#
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
#
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

############################################
################ THE DEMO ##################
############################################

###################
#~~~~ ArcGIS ~~~~~#
###################

#Item 1:
#       Create Space Time Cube 

#Item 2:
#       Emerging Hot Spot Analysis 

#Item 3:
#       Enrich Layer 

###################
#~~~~~~~ R ~~~~~~~#
###################

#PACKAGES
library(arcgisbinding)
arc.check_product()

#BUILD THE BRIDGE
enrich_df <- arc.open('C:/Users/marj8502/Documents/ArcGIS/Packages/Pro_Internal_Use_Only/commondata/rbridge_demo1/SF_Crime_Enrich.shp')
enrich_select_df <- arc.select(enrich_df, fields = c('FID', 'Shape', 'SUM_VALUE', 'N36_BUS', 'N35_BUS', 'N34_BUS', 'N39_BUS', 'N16_BUS', 'N24_BUS', 'N11_BUS', 'N13_BUS', 'N37_BUS', 'N15_BUS', 'N14_BUS', 'N19_BUS','N20_BUS', 'FedLandPt', 'NPSLandPt', 'NLCDfrstPt', 'TOTPOP10', 'ACRES_FACT'))
enrich_spdf <- arc.data2sp(enrich_select_df)
col_names <- c("OBJECTID", "Crime_Counts", "Accommodations", "Food_Service", "Entertainment", "Auto", "Clothing", "Bank", "Electronics", "Grocery", "Restaurant", "Gas", "Health", "Misc", "Non_Store", "FedLandPct", "NPSLandPct", "NLCDFrstPct", "Population", "Acres")
colnames(enrich_spdf@data) <- col_names

#RATES 
#Item 1: 
#       Compute the global empirical Bayes estimate for rates "shrunk" to the 
#       overall mean. 

n <- enrich_spdf@data$Crime_Counts
x <- enrich_spdf@data$Population

library(spdep)
EB <- EBest(n, x)
p <- EB$raw
b <- attr(EB, "parameters")$b
a <- attr(EB, "parameters")$a
v <- a + (b/x)
v[v < 0] <- b/x
z <- (p - b)/sqrt(v)

enrich_spdf@data$EB_Rate <- z

#R ARCGIS BRIDGE
#Convert from Spatial Points Data Frame to Shapefile 
arcdf_train2 <- arc.sp2data(enrich_spdf)
#Write to an Esri Shapefile
arc.write('C:/Users/marj8502/Documents/ArcGIS/Packages/Pro_Internal_Use_Only/commondata/rbridge_demo1/RATES_2014.shp', arcdf_train2)


###################
#~~~~ ArcGIS ~~~~~#
###################

#Item 1:
#       Optimized Hot Spot Analysis


###################
#~~~~~~~ R ~~~~~~~#
###################

#BUILD THE BRIDGE
rate_df <- arc.open('C:/Users/marj8502/Documents/ArcGIS/Packages/Pro_Internal_Use_Only/commondata/rbridge_demo1/CRIME_RATES.shp')
rate_select_df <- arc.select(rate_df, fields = c('FID', 'Shape', "Crime_Coun", "Accommodat", "Food_Servi", "Entertainm", "Auto", "Clothing", "Bank", "Electronic", "Grocery", "Restaurant", "Gas", "Health", "Misc", "Non_Store", "FedLandPct", "NPSLandPct", "NLCDFrstPc", "Population", "Acres", "EB_Rate"))
rate_spdf <- arc.data2sp(rate_select_df)
col_names <- c('FID', "Crime_Counts", "Accommodation", "Food_Service", "Entertainment", "Auto", "Clothing", "Bank", "Electronic", "Grocery", "Restaurant", "Gas", "Health", "Misc", "Non_Store", "FedLandPct", "NPSLandPct", "NLCDFrstPct", "Population", "Acres", "EB_Rate")
colnames(rate_spdf@data) <- col_names

#Correlation Matrix
corr_sub <- rate_spdf@data[c("EB_Rate", "Accommodation", "Food_Service", "Entertainment", "Auto", "Clothing", "Bank", "Electronic", "Grocery", "Restaurant", "Gas", "Health", "FedLandPct", "NPSLandPct", "NLCDFrstPct", "Acres")]
cormax <- round(cor(corr_sub), 2)

#Pretty it up
upper_tri <- get_upper_tri(cormax)
melted_cormax <- melt(upper_tri, na.rm = TRUE)

#Reorder 
cormax <- reorder_cormat(cormax)
upper_tri <- get_upper_tri(cormax)
melted_cormax <- melt(upper_tri, na.rm = TRUE)
ggheatmap <- ggplot(melted_cormax, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
print(ggheatmap)

#Add correlation coefficients
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

#SPATIAL WEIGHTS MATRIX
rate_list <- poly2nb(rate_spdf, snap = 300, queen = TRUE)
W <- nb2listw(rate_list, style = "W", zero.policy = TRUE)

#SIMULTANEOUS AUTOREGRESSIVE MODEL
#Maximum likelihood estimation 
sf_sar <- spautolm(EB_Rate ~ NLCDFrstPct + Electronic, data = rate_spdf@data, W, weights = Population, zero.policy = TRUE)
summary(sf_sar)

