library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(imageRy)

mato <- im.import("matogrosso_ast_2006209_lrg.jpg")

# function
# n: number of bands into which you want to divide the image
# default for n is 5

im.ridgelinecrop <- function(im, scale, palette = c(
  "viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo"),
  n) {
  
  palette <- palette[1]
  n <- 5
  
  if(!is(im, "SpatRaster")) stop("im must be a SpatRaster")
  if(!is.numeric(scale)) stop("scale must be numeric")
  if(!is.character(palette)) stop("palette must be a character")
  if(!palette %in% c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")) stop("palette must be one of the color options in the viridis package (viridis, magma, plasma, inferno, cividis, mako, rocket, turbo)")
  if(!is.numeric(n)) stop("n must be numeric")

  n <- n
  
  # image division
  e <- ext(im)
  height <- (ymax(e) - ymin(e)) / n
  bands <- vector("list", n)
  
  for (i in 1:n) {
    y_top <- ymax(e) - (i - 1) * height
    y_bottom <- ymax(e) - i * height
    sub_e <- ext(xmin(e), xmax(e), y_bottom, y_top)
    bands[[i]] <- crop(im, sub_e)
  }
  
# dataframe creation
# three columns: layer, values, band
  
  df_list <- list()
  
  for (i in seq_along(bands)) {
    df <- as.data.frame(bands[[i]]) %>%
      pivot_longer(
        cols= everything(),
        names_to = "layer",
        values_to = "values"
      ) %>%
      mutate(band = i)   
    df_list[[i]] <- df  
  }
  
  df_long <- bind_rows(df_list)
  
  
  pl <- ggplot(df_long, aes(x = values, y = layer, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = scale, rel_min_height = 0.01) +
    scale_fill_viridis_c(option = palette) +
    facet_wrap(~ band, nrow = 1) +
    theme_minimal()
  
  return(pl)
  
# ALTERNATIVE OPTION: bands in the y axis
  
#  df_list <- list()
  
#  for (i in seq_along(bands)) {
#    df <- as.data.frame(bands[[i]], wide = FALSE) %>%
#      pivot_longer(
#        cols = -c(wide),
#        names_to = "layer",
#        values_to = "values"
#      )%>%
#      mutate(band = factor(i))  
  
#    df_list[[i]] <- df  
#  }  
  
#  df_long <- bind_rows(df_list)
  
#  pl <- ggplot(df_long, aes(x = values, y = band, fill = after_stat(x))) + # slice al posto di layer (opz 2)
#    geom_density_ridges_gradient(scale = scale, rel_min_height = 0.01) +
#    scale_fill_viridis_c(option = palette) +
#    facet_wrap(~ layer, nrow = 1) +
#    theme_minimal()
  
# return(pl)
  
}  

im.ridgelinecrop(mato,1,"viridis",n=5)

im.ridgelinecrop(mato,1,"viridis")
