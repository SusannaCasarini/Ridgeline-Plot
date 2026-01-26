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
im.ridgelinecroptime <- function(im_list, scale, n=5, direction="horizontal") {
  
  if(!is.list(im_list)) stop("im_list must be a list of raster images")
  if(!is.numeric(scale)) stop("scale must be numeric")
  if(!is.numeric(n)) stop("n must be numeric")
  direction <- match.arg(direction, choices = c("horizontal", "vertical"))
  
  # Liste per salvare tutti i dataframe e i plots
  all_df_list <- list()
  all_plots_list <- list()
  
  # To manage one-image lists
  if (is.null(names(im_list))) {
    names(im_list) <- "SingleTime Image"
  }
  
  
  # Loop su tutte le immagini
  for (t in seq_along(im_list)) {
    im <- im_list[[t]]
    time_label <- names(im_list)[t]  
    
    e <- ext(im)
    bands <- vector("list", n)
    plots <- list()
    
    if (direction == "horizontal") {
      height <- (ymax(e) - ymin(e)) / n
      for (i in 1:n) {
        y_top <- ymax(e) - (i - 1) * height
        y_bottom <- ymax(e) - i * height
        sub_e <- ext(xmin(e), xmax(e), y_bottom, y_top)
        bands[[i]] <- crop(im, sub_e)
        plots[[paste0("band", i)]] <- local({
          b <- bands[[i]]
          ii <- i
          function() {
            plot(b, main = paste("Band", ii))
          }
        })
      }
    } else if (direction == "vertical") {
      width <- (xmax(e) - xmin(e)) / n
      for (i in 1:n) {
        x_left  <- xmin(e) + (i - 1) * width
        x_right <- xmin(e) + i * width
        sub_e <- ext(x_left, x_right, ymin(e), ymax(e))
        bands[[i]] <- crop(im, sub_e)
        
        plots[[paste0("band", i)]] <- local({
          b <- bands[[i]]
          ii <- i
          function() {
            plot(b, main = paste("Band", ii))
          }
        })
      }
    }
    
    # Creazione dataframe per ciascuna immagine
    df_list <- list()
    for (i in seq_along(bands)) {
      df <- as.data.frame(bands[[i]]) %>%
        pivot_longer(
          cols= everything(),
          names_to = "layer",
          values_to = "values"
        ) %>%
        mutate(band = i,
               time = time_label)  # aggiungo la variabile time
      df_list[[i]] <- df
    }
    all_df_list[[t]] <- bind_rows(df_list)
    all_plots_list[[time_label]] <- plots
    
  }
  
  # Unisco tutti i dati
  df_long <- bind_rows(all_df_list)
  
  # Plot ridgeline con faceting per immagine e banda
  
  if (n!=1){
    pl <- ggplot(df_long, aes(x = values, y = layer, fill = layer)) +
      geom_density_ridges(scale = scale, rel_min_height = 0.01, alpha=0.5) +
      scale_fill_viridis_d(option = "viridis") +
      facet_grid(band ~ time) + 
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 9),
        panel.spacing = unit(0.1, "lines"),
        labs(
          title = "Bands ridgeline plots for time and space",
          x = "Time",
          y = "Layer",
          fill = "Layer2",
        )
        
      )
    
  } else if (n==1){
    pl <- ggplot(df_long, aes(x = values, y = layer, fill = layer)) +
      geom_density_ridges(scale = scale, rel_min_height = 0.01, alpha=0.5) +
      scale_fill_viridis_d(option = "viridis") +
      facet_grid(time~.) + 
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8),
        strip.text = element_text(size = 9),
        panel.spacing = unit(0.1, "lines"),
        labs(
          title = "Bands ridgeline plots for time",
          x = "Vaues",
          y = "Time",
          fill = "Layer2",
        )
        
      )
  }
  
  
  return(list(
    pl = pl,
    plots = all_plots_list
  ))
}

# image download: B2, B3, B4, B8, B11

library(tools)
files <- list.files(path, pattern = "\\.tif$", full.names=TRUE)
for (f in files){
  nome <- file_path_sans_ext(basename(f))
  assign(nome,rast(f))
}

###################
## NDVI analysis ##
###################

# NDVI, NDWI and NDBI computation
# NDVI = (NIR-RED)/(NIR+RED)
# NDWI = (GREEN-NIR)/(GREEN+NIR)
# NDBI = (SWIR-RED)/(SWIR+RED)

mesi <- c("jan", "feb", "mar", "apr", "may", "jun",
          "jul", "aug", "sep", "oct", "nov", "dec")

for (m in mesi) {
  
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
  
  par(mfrow = c(2, 2))
  plot(ndvi, main = paste("NDVI", toupper(m)))
  plot(ndwi, main = paste("NDWI", toupper(m)))
  plot(ndbi, main = paste("NDBI", toupper(m)))
  im.plotRGB(RGB, 1, 2, 3)
  par(mfrow = c(1, 1))
}

# NDVI analysis in time
raster_names <- c("ndvimar","ndvijun","ndvisep","ndvidec")

raster_list <- mget(raster_names)
ndvitime <- im.ridgelinecroptime(raster_list,4,1,direction = "vertical")
ndvitime$pl

# NDVI analysis in time and space cutting the image
ndwitimecrop <- im.ridgelinecroptime(raster_list,0.5,4,direction = "vertical")
ndwitimecrop$pl

# NDWI analysis in time
raster_names <- c("ndwimar","ndwijun","ndwisep","ndwidec")
raster_list <- mget(raster_names)
ndwitime <- im.ridgelinecroptime(raster_list,3,1,direction = "vertical")
ndwitime$pl

# NDWI analysis in time and space cutting the image
ndwitimecrop <- im.ridgelinecroptime(raster_list,3,4,direction = "vertical")
ndwitimecrop$pl

# NDBI analysis in time
raster_names <- c("ndbimar","ndbijun","ndbisep","ndbidec")
raster_list <- mget(raster_names)
ndbitime <- im.ridgelinecroptime(raster_list,3,1,direction = "vertical")
ndbitime$pl

# NDBI analysis in time and space cutting the image
ndbitimecrop <- im.ridgelinecroptime(raster_list,3,4,direction = "vertical")
ndbitimecrop$pl

##################
## RGB analysis ##
##################

# RGB analysis in time
raster_names <- c("marRGB","junRGB","sepRGB","decRGB")
raster_list <- mget(raster_names)
monthstime <- im.ridgelinecroptime(raster_list,1,1,direction = "vertical")
monthstime$pl

# RGB analysis in time and space cutting the image
ndbitimecrop <- im.ridgelinecroptime(raster_list,1,4,direction = "vertical")
ndbitimecrop$pl

######################
## Cluster Analysis ##
######################

# Cluster analysis using the previous informations
# Ridgeline plots from cluster analysis

im.classify(junRGB,4)

