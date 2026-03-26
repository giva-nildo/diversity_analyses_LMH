### =================================
### Obtaining and saving the mini map
### =================================
# Givanildo Rodrigues da Silva
# 03/25/2026

library(dplyr)
library(ggplot2)
library(sf)
install.packages("rnaturalearthdata")
library(rnaturalearthdata)
library(countrycode)
library(tibble)

country_df <- read_xlsx("Supplementary Tables.xlsx", sheet = "Supplementary Table S1", range = "A2:M452")
country_df <- country_df[,8]


# 1) Unique countries in your data
country_present <- country_df %>%
  distinct(Origin) %>%
  mutate(
    iso3 = countrycode(Origin, origin = "country.name", destination = "iso3c")
  )

# optional: check unmatched names
country_present %>% filter(is.na(iso3))

# 2) World map
world <- ne_countries(scale = "medium", returnclass = "sf")

# 3) Join and create binary flag
world_map <- world %>%
  left_join(country_present, by = c("iso_a3" = "iso3")) %>%
  mutate(in_data = !is.na(Origin))

#-----------------------------
# 3) Mini map similar to your example
#-----------------------------
mini_map <- ggplot() +
  geom_sf(data = world_map, fill = "grey92", color = "white", linewidth = 0.15) +
  geom_sf(
    data = subset(world_map, in_data),
    fill = "grey65",
    color = "white",
    linewidth = 0.2
  ) +
  coord_sf(
    expand = T
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

mini_map

ggsave(
  filename = "map.pdf",
  plot = mini_map,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)
