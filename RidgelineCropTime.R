library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(imageRy)
library(purrr)
library(plotly)


# function

im.ridgelinecroptime <- function(im_list, scale, n=5, direction="horizontal") {
  
  if(!is.list(im_list)) stop("im_list must be a list of raster images")
  if(!is.numeric(scale)) stop("scale must be numeric")
  if(!is.numeric(n)) stop("n must be numeric")
  direction <- match.arg(direction, choices = c("horizontal", "vertical"))
  
  all_df_list <- list()
  all_plots_list <- list()

  # To manage one-image lists
  if (is.null(names(im_list))) {
    names(im_list) <- "SingleTime Image"
  }
  
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
    
    # dataframe for each image
    df_list <- list()
    for (i in seq_along(bands)) {
      df <- as.data.frame(bands[[i]]) %>%
        pivot_longer(
          cols= everything(),
          names_to = "layer",
          values_to = "values"
        ) %>%
        mutate(band = i,
               time = time_label)
      df_list[[i]] <- df
    }
    all_df_list[[t]] <- bind_rows(df_list)
    all_plots_list[[time_label]] <- plots
    
  }
  
  # Unisco tutti i dati
  df_long <- bind_rows(all_df_list)
  
  # Plot ridgeline con faceting per immagine e banda
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
        y = "Layer"
      )
      
    )
    
  
  return(list(
    pl = pl,
    plots = all_plots_list
  ))
}

# application 

# importing the data
en1 <- im.import("EN_01.png")
en2 <- im.import("EN_02.png")
en3 <- im.import("EN_03.png")
en4 <- im.import("EN_04.png")
en5 <- im.import("EN_05.png")
en6 <- im.import("EN_06.png")
en7 <- im.import("EN_07.png")
en8 <- im.import("EN_08.png")
en9 <- im.import("EN_09.png")
en10 <- im.import("EN_10.png")
en11 <- im.import("EN_11.png")
en12 <- im.import("EN_12.png")

# the input im_list has to be a list of raster objects with the same name for the bands
names(en1) <- c("B2", "B3", "B4")
names(en2) <- c("B2", "B3", "B4")
names(en3) <- c("B2", "B3", "B4")
names(en4) <- c("B2", "B3", "B4")
names(en5) <- c("B2", "B3", "B4")
names(en6) <- c("B2", "B3", "B4")
names(en7) <- c("B2", "B3", "B4")
names(en8) <- c("B2", "B3", "B4")
names(en9) <- c("B2", "B3", "B4")
names(en10) <- c("B2", "B3", "B4")
names(en11) <- c("B2", "B3", "B4")
names(en12) <- c("B2", "B3", "B4")

raster_names <- c("en1","en2","en3","en4","en5","en6","en7","en8","en9","en10","en11","en12")
raster_list <- mget(raster_names)
class(raster_list)

# apply the function
enf <- im.ridgelinecroptime(raster_list,2,n=4,direction = "horizontal")

# recall the space-time plot
enf$pl

# recall the single cropped images
enf$plots$en1$band1()
