# libraries
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(imageRy)
library(purrr)
library(plotly)

# set working directory
setwd("C:/Suz/Università/Data science/Internship/Ridgeline/files")
path <- "C:/Suz/Università/Data science/Internship/Ridgeline/files"

# function im.ridgelinecroptime
im.ridgelinecrop <- function(im_list, scale, n=5, direction="horizontal", time_axis="y") {
  
  if(!is.list(im_list)) stop("im_list must be a list of raster images")
  if(!is.numeric(scale)) stop("scale must be numeric")
  if(!is.numeric(n)) stop("n must be numeric")
  direction <- match.arg(direction, choices = c("horizontal", "vertical"))
  time_axis <- match.arg(time_axis, choices = c("y", "x"))
  
  # Lists to save dataframes and plots
  all_df_list <- list()
  all_plots_list <- list()
  all_partitions_list <- list()

  # To manage one-image lists
  if (is.null(names(im_list))) {
    names(im_list) <- "SingleTime Image"
  }
  
  for (t in seq_along(im_list)) {
    im <- im_list[[t]]
    time_label <- names(im_list)[t]  
    
    e <- ext(im)
    partitions <- vector("list", n)
    plots <- list()
    partitions_image <- list()
    
    if (direction == "horizontal") {
      height <- (ymax(e) - ymin(e)) / n
      for (i in 1:n) {
        y_top <- ymax(e) - (i - 1) * height
        y_bottom <- ymax(e) - i * height
        sub_e <- ext(xmin(e), xmax(e), y_bottom, y_top)
        partitions[[i]] <- crop(im, sub_e)
        plots[[paste0("partition", i)]] <- local({
          p <- partitions[[i]]
          ii <- i
          function() {
            plot(p, main = paste("Partition", ii))
          }
        })
        partitions_image[[paste0("partition", i)]] <- local({
          p <- partitions[[i]]
        })
      }
    } else if (direction == "vertical") {
      width <- (xmax(e) - xmin(e)) / n
      for (i in 1:n) {
        x_left  <- xmin(e) + (i - 1) * width
        x_right <- xmin(e) + i * width
        sub_e <- ext(x_left, x_right, ymin(e), ymax(e))
        partitions[[i]] <- crop(im, sub_e)
        
        plots[[paste0("partition", i)]] <- local({
          p <- partitions[[i]]
          ii <- i
          function() {
            plot(p, main = paste("Partition", ii))
          }
        })
        partitions_image[[paste0("partition", i)]] <- local({
          p <- partitions[[i]]
        })
      }
    }
    
    df_list <- list()
    for (i in seq_along(partitions)) {
      df <- as.data.frame(partitions[[i]]) %>%
        pivot_longer(
          cols= everything(),
          names_to = "layer",
          values_to = "values"
        ) %>%
        mutate(partition = i,
               time = time_label)  
      df_list[[i]] <- df
    }
    all_df_list[[t]] <- bind_rows(df_list)
    all_plots_list[[time_label]] <- plots
    all_partitions_list[[time_label]] <- partitions_image

  }
  
  df_long <- bind_rows(all_df_list)
  
  df_long$time <- factor(
    df_long$time,
    levels = names(im_list)
  )
  
  if (time_axis=="x"){
    pl <- ggplot(df_long, aes(x = values, y = layer, fill = layer)) +
      geom_density_ridges(scale = scale, rel_min_height = 0.01, alpha=0.5) +
      scale_fill_viridis_d(option = "viridis") +
      facet_grid(partition ~ time, space = "free_y") + 
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 9),
        panel.spacing.y = unit(0, "cm"),
        panel.spacing.x = unit(1,"cm")
      )
    
  } else if (time_axis=="y"){
    pl <- ggplot(df_long, aes(x = values, y = layer, fill = layer)) +
      geom_density_ridges(scale = scale, rel_min_height = 0.01, alpha=0.5) +
      scale_fill_viridis_d(option = "viridis") +
      facet_grid(time ~ partition, space = "free_y") + 
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 9),
        panel.spacing.y = unit(0, "cm"),
        panel.spacing.x = unit(1,"cm")
      )
  }
  
  
  return(list(
    pl = pl,
    plots = all_plots_list,
    partitions_image = all_partitions_list
  ))
}


# image download: B2, B3, B4, B8, B11
library(tools)
files <- list.files(path, pattern = "\\.tif$", full.names=TRUE)
for (f in files){
  nome <- file_path_sans_ext(basename(f))
  assign(nome,rast(f))
}

###############################
## Spectral Indices Analysis ##
###############################


# Observing RGB images of the 12 months
#mesi <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
for (m in mesi) {
  
  RGB <- get(paste0(m, "RGB"))   

  BLUE  <- RGB[[1]]   
  GREEN <- RGB[[2]]   
  RED   <- RGB[[3]] 
  
  plotRGB(RGB, 1, 2, 3, stretch="hist")
}

# selection of the seasonal images and computation of the spectral indices
seasons <- c("jan","mar", "jun","oct")
for (m in seasons) {
  
  RGB <- get(paste0(m, "RGB"))   
  B8  <- get(paste0(m, "B8"))    
  B11 <- get(paste0(m, "B11"))   
  
  BLUE  <- RGB[[1]]   
  GREEN <- RGB[[2]]   
  RED   <- RGB[[3]] 
  
  ndvi <- (B8 - RED) / (B8 + RED)
  ndwi <- (GREEN - B8) / (GREEN + B8)
  ndbi <- (B11 - B8) / (B11 + B8)
  
  assign(paste0("ndvi", m), ndvi)
  assign(paste0("ndwi", m), ndwi)
  assign(paste0("ndbi", m), ndbi)
  
  plotndvi <- plot(ndvi, main = paste("NDVI", toupper(m)))
  plotndwi <- plot(ndwi, main = paste("NDWI", toupper(m)))
  plotndbi <- plot(ndbi, main = paste("NDBI", toupper(m)))
  plotrgb <- plotRGB(RGB, 1, 2, 3, stretch="hist")
  
  assign(paste0("plotndvi", m), plotndvi)
  assign(paste0("plotndwi", m), plotndwi)
  assign(paste0("plotndbi", m), plotndbi)
  assign(paste0("plotrgb", m), plotrgb)
  
}

names(ndvijan) <- "NDVI"
names(ndvimar) <- "NDVI"
names(ndvijun) <- "NDVI"
names(ndvioct) <- "NDVI"

names(ndwijan) <- "NDWI"
names(ndwimar) <- "NDWI"
names(ndwijun) <- "NDWI"
names(ndwioct) <- "NDWI"

names(ndbijan) <- "NDBI"
names(ndbimar) <- "NDBI"
names(ndbijun) <- "NDBI"
names(ndbioct) <- "NDBI"

# NDVI = (NIR-RED)/(NIR+RED)
January <-ndvijan
March <-ndvimar
June <-ndvijun
October <-ndvioct
raster_names_ndvi <- c("January","March","June","October")
raster_list_ndvi <- mget(raster_names_ndvi)

ndvitime <- im.ridgelinecrop(raster_list_ndvi,4,1,direction = "vertical")
ndvitime$pl +
  labs(title = "NDVI Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )

ndvitimecrop <- im.ridgelinecrop(raster_list_ndvi,3,4,direction = "vertical")
ndvitimecrop$pl+
  labs(title = "NDVI Partitioned Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )

# NDWI = (GREEN-NIR)/(GREEN+NIR)
January <-ndwijan
March <-ndwimar
June <-ndwijun
October <-ndwioct
raster_names_ndwi <- c("January","March","June","October")
raster_list_ndwi <- mget(raster_names_ndwi)
ndwitime <- im.ridgelinecrop(raster_list_ndwi,3,1,direction = "vertical")
ndwitime$pl+
  labs(title = "NDWI Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )

ndwitimecrop <- im.ridgelinecrop(raster_list_ndwi,3,4,direction = "vertical")
ndwitimecrop$pl+
  labs(title = "NDWI Partitioned Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )

# NDBI = (SWIR-RED)/(SWIR+RED)
January <-ndbijan
March <-ndbimar
June <-ndbijun
October <-ndbioct
raster_names_ndbi <- c("January","March","June","October")
raster_list_ndbi <- mget(raster_names_ndbi)
ndbitime <- im.ridgelinecrop(raster_list_ndbi,3,1,direction = "vertical")
ndbitime$pl+
  labs(title = "NDBI Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )

ndbitimecrop <- im.ridgelinecrop(raster_list_ndbi,3,4,direction = "vertical")
ndbitimecrop$pl+
  labs(title = "NDBI Partitioned Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )


##################
## RGB analysis ##
##################

# RGB analysis just to show the possibility
# It won't be the focus of the thesis

#raster_names <- c("junRGB","decRGB")
#raster_list <- mget(raster_names)

#monthstimecrop <- im.ridgelinecroptime(raster_list,1,4,direction = "vertical")
#monthstimecrop$pl

##########################
### Spectral as layers ###
##########################

# to have a comprehensive view of the spectral indices
# it is possible to set the indices as layers and
# to use the function to plot the cropped images together
# with the time and the spectral indices

spectral_mesi <- list()
seasons <- c("jan","mar","jun","oct")

for (m in seasons) {
  
  ndvi <- get(paste0("ndvi", m))
  ndwi <- get(paste0("ndwi", m))
  ndbi <- get(paste0("ndbi", m))
  
  month <- c(ndvi, ndwi, ndbi)
  names(month) <- c("ndvi", "ndwi", "ndbi")
  
  spectral_mesi[[m]] <- month
}


January <- spectral_mesi$jan
March <- spectral_mesi$mar
June <- spectral_mesi$jun
October <- spectral_mesi$oct

raster_names_spectral <- c("January","March","June","October")
spectral_months <- mget(raster_names_spectral)

spectraltime <- im.ridgelinecrop(spectral_months,3,1,"vertical")
spectraltime$pl+
  labs(title = "Spectral Indices Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )

spectral <- im.ridgelinecrop(spectral_months,3,4,"vertical")
spectral$pl+
  labs(title = "Spectral Indices Partitioned Ridgeline Plots") +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
    )
  )

# plots of the RGB cropped images
#for (s in seasons){
#  image <- get(paste0(s,"RGB"))
#  plotRGB(image, 1,2,3, stretch="hist")
#  e <- ext(image)
#  x_lines <- seq(e$xmin, e$xmax, length.out = 5)[-c(1,5)]
#  abline(v = x_lines, col = "white", lwd = 5)
#}

#plotRGB(junRGB, 1,2,3, stretch="hist")

####################
#### Clustering ####
####################


#ndvi
jan_partition1 <- ndvitimecrop$partitions_image$January$partition1
jan_partition2 <- ndvitimecrop$partitions_image$January$partition2
jan_partition3 <- ndvitimecrop$partitions_image$January$partition3
jan_partition4 <- ndvitimecrop$partitions_image$January$partition4
mar_partition1 <- ndvitimecrop$partitions_image$March$partition1
mar_partition2 <- ndvitimecrop$partitions_image$March$partition2
mar_partition3 <- ndvitimecrop$partitions_image$March$partition3
mar_partition4 <- ndvitimecrop$partitions_image$March$partition4
jun_partition1 <- ndvitimecrop$partitions_image$June$partition1
jun_partition2 <- ndvitimecrop$partitions_image$June$partition2
jun_partition3 <- ndvitimecrop$partitions_image$June$partition3
jun_partition4 <- ndvitimecrop$partitions_image$June$partition4
oct_partition1 <- ndvitimecrop$partitions_image$October$partition1
oct_partition2 <- ndvitimecrop$partitions_image$October$partition2
oct_partition3 <- ndvitimecrop$partitions_image$October$partition3
oct_partition4 <- ndvitimecrop$partitions_image$October$partition4

############################
#Silhouette
############################
library(cluster)
# confrontre il numero di cluster con silhouette: hierarchical e k-means
# sample

par(ann=TRUE, mar=c(5,4,4,2)+0.1, mgp=c(3,1,0))

image_values <- na.omit(terra::as.matrix(ndvijan))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical January Full Image")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans January Full Image")


image_values <- na.omit(terra::as.matrix(jan_partition1))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical January Partition 1")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans January Partition 1")


image_values <- na.omit(terra::as.matrix(jan_partition2))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical January Partition 2")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans January Partition 2")

image_values <- na.omit(terra::as.matrix(jan_partition3))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical January Partition 3")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans January Partition 3")

image_values <- na.omit(terra::as.matrix(jan_partition4))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical January Partition 4")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans January Partition 4")

# june

# confrontre il numero di cluster con silhouette: hierarchical e k-means
# sample

image_values <- na.omit(terra::as.matrix(ndvijun))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical June Full Image")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans June Full Image")

image_values <- na.omit(terra::as.matrix(jun_partition1))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical June Partition 1")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans June Partition 1")

image_values <- na.omit(terra::as.matrix(jun_partition2))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical June Partition 2")


# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans June Partition 2")

image_values <- na.omit(terra::as.matrix(jun_partition3))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical June Partition 3")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans June Partition 3")

image_values <- na.omit(terra::as.matrix(jun_partition4))
set.seed(1234)
sample_data <- image_values[sample(nrow(image_values),5000), c("NDVI")]

# uso hierarchical per scegliere i clusters
d <- dist(sample_data, method = "euclidean")
hc <- hclust(d, method = "ward.D2")

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:6){
  tclusk[[k]] <- cutree(hc,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=d)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=d))$avg.width
}
plot(2:6,tasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette Hierarchical June Partition 4")

# kmeans with euclidean
pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:6){
  pclusk[[k]] <- kmeans(sample_data,k,nstart = 100)
  psil[[k]] <- silhouette(pclusk[[k]]$cluster,d)
  pasw[k] <- summary(psil[[k]])$avg.width
}

plot(2:6,pasw[2:6],type="l",xlab="Number of clusters",ylab="ASW",main = "Silhouette kmeans June Partition 4")


# January 4 clusters
# NDVI clustering for band 1
library(cluster)

# partition 1
image_values <- na.omit(terra::as.matrix(jan_partition1))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jan_partition1[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_1 <- kmeans_result$cluster

# partition 2
image_values <- na.omit(terra::as.matrix(jan_partition2))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jan_partition2[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_2 <- kmeans_result$cluster

# partition 3
image_values <- na.omit(terra::as.matrix(jan_partition3))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jan_partition3[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_3 <- kmeans_result$cluster


# partition 4
image_values <- na.omit(terra::as.matrix(jan_partition4))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jan_partition4[[1]]
values(classified_image) <- kmeans_result$cluster
plot(classified_image, axes = FALSE)
cluster_vec_4 <- kmeans_result$cluster

# full image
image_values <- na.omit(terra::as.matrix(ndvijan))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- ndvijan[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
full_classified <- plot(classified_image, axes = FALSE)
cluster_vec_full <- kmeans_result$cluster

# plot of the diferent indices
cluster_perc1_jan4 <- prop.table(table(cluster_vec_1))*100
cluster_perc2_jan4 <- prop.table(table(cluster_vec_2))*100
cluster_perc3_jan4 <- prop.table(table(cluster_vec_3))*100
cluster_perc4_jan4 <- prop.table(table(cluster_vec_4))*100
cluster_percfull_jan4 <- prop.table(table(cluster_vec_full))*100

# divido il plot del cluster dell full image e guardo la percentuale di pixel nei clusters
e <- ext(classified_image)
width <- (xmax(e) - xmin(e)) / 4
part <- vector("list", 4)

for (i in 1:4) {
    x_left  <- xmin(e) + (i - 1) * width
    x_right <- xmin(e) + i * width
    sub_e <- ext(x_left, x_right, ymin(e), ymax(e))
    part[[i]] <- crop(classified_image, sub_e)
}

part1 <- part[[1]]
part2 <- part[[2]]
part3 <- part[[3]]
part4 <- part[[4]]

plot(part1)

freq_table1 <- as.data.frame(freq(part1))
total_pixels <- sum(freq_table1$count)
freq_table1$percentage_jan4 <- round(100 * freq_table1$count / total_pixels,2)

freq_table2 <- as.data.frame(freq(part2))
total_pixels <- sum(freq_table2$count)
freq_table2$percentage_jan4 <- round(100 * freq_table2$count / total_pixels,2)

freq_table3 <- as.data.frame(freq(part3))
total_pixels <- sum(freq_table3$count)
freq_table3$percentage_jan4 <- round(100 * freq_table3$count / total_pixels,2)

freq_table4 <- as.data.frame(freq(part4))
total_pixels <- sum(freq_table4$count)
freq_table4$percentage_jan4 <- round(100 * freq_table4$count / total_pixels,2)

# confronti
cluster_percfull_jan4

freq_table1$percentage_jan4
cluster_perc1_jan4
freq_table2$percentage_jan4
cluster_perc2_jan4
freq_table3$percentage_jan4
cluster_perc3_jan4
freq_table4$percentage_jan4
cluster_perc4_jan4



# June 4 clusters
# NDVI clustering for band 1

# partition 1
image_values <- na.omit(terra::as.matrix(jun_partition1))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jun_partition1[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_1 <- kmeans_result$cluster

# partition 2
image_values <- na.omit(terra::as.matrix(jun_partition2))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jun_partition2[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_2 <- kmeans_result$cluster

# partition 3
image_values <- na.omit(terra::as.matrix(jun_partition3))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jun_partition3[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_3 <- kmeans_result$cluster

# partition 4
image_values <- na.omit(terra::as.matrix(jun_partition4))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- jun_partition4[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_4 <- kmeans_result$cluster

# full image
image_values <- na.omit(terra::as.matrix(ndvijun))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 4)
classified_image <- ndvijun[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
full_classified <- plot(classified_image, axes = FALSE)
cluster_vec_full <- kmeans_result$cluster


# plot of the diferent indices
cluster_perc1_jun4 <- prop.table(table(cluster_vec_1))*100
cluster_perc2_jun4 <- prop.table(table(cluster_vec_2))*100
cluster_perc3_jun4 <- prop.table(table(cluster_vec_3))*100
cluster_perc4_jun4 <- prop.table(table(cluster_vec_4))*100
cluster_percfull_jun4 <- prop.table(table(cluster_vec_full))*100

# divido il plot del cluster dell full image e guardo la percentuale di pixel nei clusters
e <- ext(classified_image)
width <- (xmax(e) - xmin(e)) / 4
part <- vector("list", 4)

for (i in 1:4) {
  x_left  <- xmin(e) + (i - 1) * width
  x_right <- xmin(e) + i * width
  sub_e <- ext(x_left, x_right, ymin(e), ymax(e))
  part[[i]] <- crop(classified_image, sub_e)
}
part1 <- part[[1]]
part2 <- part[[2]]
part3 <- part[[3]]
part4 <- part[[4]]

freq_table1 <- as.data.frame(freq(part1))
total_pixels <- sum(freq_table1$count)
freq_table1$percentage_jun4 <- round(100 * freq_table1$count / total_pixels,2)

freq_table2 <- as.data.frame(freq(part2))
total_pixels <- sum(freq_table2$count)
freq_table2$percentage_jun4 <- round(100 * freq_table2$count / total_pixels,2)

freq_table3 <- as.data.frame(freq(part3))
total_pixels <- sum(freq_table3$count)
freq_table3$percentage_jun4 <- round(100 * freq_table3$count / total_pixels,2)

freq_table4 <- as.data.frame(freq(part4))
total_pixels <- sum(freq_table4$count)
freq_table4$percentage_jun4 <- round(100 * freq_table4$count / total_pixels,2)

# confronti
cluster_percfull_jun4

freq_table1$percentage_jun4
cluster_perc1_jun4
freq_table2$percentage_jun4
cluster_perc2_jun4
freq_table3$percentage_jun4
cluster_perc3_jun4
freq_table4$percentage_jun4
cluster_perc4_jun4



# January 3 clusters
# NDVI clustering for band 1
library(cluster)

# partition 1
image_values <- na.omit(terra::as.matrix(jan_partition1))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jan_partition1[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_1 <- kmeans_result$cluster

# partition 2
image_values <- na.omit(terra::as.matrix(jan_partition2))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jan_partition2[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_2 <- kmeans_result$cluster

# partition 3
image_values <- na.omit(terra::as.matrix(jan_partition3))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jan_partition3[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_3 <- kmeans_result$cluster

# partition 4
image_values <- na.omit(terra::as.matrix(jan_partition4))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jan_partition4[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_4 <- kmeans_result$cluster

# full image
image_values <- na.omit(terra::as.matrix(ndvijan))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- ndvijan[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
full_classified <- plot(classified_image, axes = FALSE)
cluster_vec_full <- kmeans_result$cluster


# plot of the diferent indices
cluster_perc1_jan3 <- prop.table(table(cluster_vec_1))*100
cluster_perc2_jan3 <- prop.table(table(cluster_vec_2))*100
cluster_perc3_jan3 <- prop.table(table(cluster_vec_3))*100
cluster_perc4_jan3 <- prop.table(table(cluster_vec_4))*100
cluster_percfull_jan3 <- prop.table(table(cluster_vec_full))*100

# divido il plot del cluster dell full image e guardo la percentuale di pixel nei clusters
e <- ext(classified_image)
width <- (xmax(e) - xmin(e)) / 4
part <- vector("list", 4)

for (i in 1:4) {
  x_left  <- xmin(e) + (i - 1) * width
  x_right <- xmin(e) + i * width
  sub_e <- ext(x_left, x_right, ymin(e), ymax(e))
  part[[i]] <- crop(classified_image, sub_e)
}
part1 <- part[[1]]
plot(part4)
part2 <- part[[2]]
part3 <- part[[3]]
part4 <- part[[4]]

freq_table1 <- as.data.frame(freq(part1))
total_pixels <- sum(freq_table1$count)
freq_table1$percentage_jan3 <- round(100 * freq_table1$count / total_pixels,2)

freq_table2 <- as.data.frame(freq(part2))
total_pixels <- sum(freq_table2$count)
freq_table2$percentage_jan3 <- round(100 * freq_table2$count / total_pixels,2)

freq_table3 <- as.data.frame(freq(part3))
total_pixels <- sum(freq_table3$count)
freq_table3$percentage_jan3 <- round(100 * freq_table3$count / total_pixels,2)

freq_table4 <- as.data.frame(freq(part4))
total_pixels <- sum(freq_table4$count)
freq_table4$percentage_jan3 <- round(100 * freq_table4$count / total_pixels,2)

# confronti
cluster_percfull_jan3

freq_table1$percentage_jan3
cluster_perc1_jan3
freq_table2$percentage_jan3
cluster_perc2_jan3
freq_table3$percentage_jan3
cluster_perc3_jan3
freq_table4$percentage_jan3
cluster_perc4_jan3



# June 3 clusters
# NDVI clustering for band 1
library(cluster)


# partition 1
image_values <- na.omit(terra::as.matrix(jun_partition1))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jun_partition1[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_1 <- kmeans_result$cluster

# partition 2
image_values <- na.omit(terra::as.matrix(jun_partition2))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jun_partition2[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_2 <- kmeans_result$cluster

# partition 3
image_values <- na.omit(terra::as.matrix(jun_partition3))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jun_partition3[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_3 <- kmeans_result$cluster

# partition 4
image_values <- na.omit(terra::as.matrix(jun_partition4))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- jun_partition4[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
cluster_vec_4 <- kmeans_result$cluster

# full image
image_values <- na.omit(terra::as.matrix(ndvijun))
set.seed(1234)
kmeans_result <- kmeans(image_values, centers = 3)
classified_image <- ndvijun[[1]]
values(classified_image) <- kmeans_result$cluster

plot(classified_image, axes = FALSE)
full_classified <- plot(classified_image, axes = FALSE)
cluster_vec_full <- kmeans_result$cluster

# plot of the diferent indices
cluster_perc1_jun3 <- prop.table(table(cluster_vec_1))*100
cluster_perc2_jun3 <- prop.table(table(cluster_vec_2))*100
cluster_perc3_jun3 <- prop.table(table(cluster_vec_3))*100
cluster_perc4_jun3 <- prop.table(table(cluster_vec_4))*100
cluster_percfull_jun3 <- prop.table(table(cluster_vec_full))*100

# divido il plot del cluster dell full image e guardo la percentuale di pixel nei clusters
e <- ext(classified_image)
width <- (xmax(e) - xmin(e)) / 4
part <- vector("list", 4)

for (i in 1:4) {
  x_left  <- xmin(e) + (i - 1) * width
  x_right <- xmin(e) + i * width
  sub_e <- ext(x_left, x_right, ymin(e), ymax(e))
  part[[i]] <- crop(classified_image, sub_e)
}
part1 <- part[[1]]
part2 <- part[[2]]
part3 <- part[[3]]
part4 <- part[[4]]

freq_table1 <- as.data.frame(freq(part1))
total_pixels <- sum(freq_table1$count)
freq_table1$percentage_jun3 <- round(100 * freq_table1$count / total_pixels,2)

freq_table2 <- as.data.frame(freq(part2))
total_pixels <- sum(freq_table2$count)
freq_table2$percentage_jun3 <- round(100 * freq_table2$count / total_pixels,2)

freq_table3 <- as.data.frame(freq(part3))
total_pixels <- sum(freq_table3$count)
freq_table3$percentage_jun3 <- round(100 * freq_table3$count / total_pixels,2)

freq_table4 <- as.data.frame(freq(part4))
total_pixels <- sum(freq_table4$count)
freq_table4$percentage_jun3 <- round(100 * freq_table4$count / total_pixels,2)

# confronti
cluster_percfull_jun3

freq_table1$percentage_jun3
cluster_perc1_jun3
freq_table2$percentage_jun3
cluster_perc2_jun3
freq_table3$percentage_jun3
cluster_perc3_jun3
freq_table4$percentage_jun3
cluster_perc4_jun3
