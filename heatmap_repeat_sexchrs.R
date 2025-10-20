library(tidyverse)

# Preparar datos
enrichment_data <- read_tsv("enrichment.tsv")
species_order <-rev(read_tsv("birds.renamed.lst")

  region_order <- c("Start","End","Middle","Total")
heatmap_data <- enrichment_data %>%
  select(species, starts_with("enrichment_")) %>%
  pivot_longer(cols = starts_with("enrichment_"),
               names_to = "region",
               values_to = "enrichment") %>%
  mutate(
    chromosome = case_when(
      str_detect(region, "_Z$") ~ "Z",
      str_detect(region, "_W$") ~ "W",
      TRUE ~ "Other"
    ),
    region_clean = str_replace(region, "enrichment_", "") %>%
      str_replace("_Z$|_W$", "") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    # Transformación logarítmica - agregar 1 para evitar log(0)
    enrichment_log = log10(enrichment + 1)
  ) %>%
  group_by(species) %>%
  mutate(
    region_clean= factor(region_clean,levels=region_order)
  ) %>%
  mutate(
    species=factor(species,levels=species_order)
  )
  

# Heatmap con escala log
ggplot(heatmap_data, aes(x = region_clean, y = species, fill = enrichment_log)) +
  geom_tile(color = "white", linewidth = 0.5) +
  facet_grid(. ~ chromosome, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "white",
    mid = "red", 
    high = "darkred",
    midpoint = 1.5,# log10(1) = 0
    name = "log10(Enrichment + 1)",
    na.value = "grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Region",
    y = "Species",
    title = "Enrichment in Z and W Chromosomes (Log Scale)"
  )

ggsave("enrichment_heatmap_log.svg", width = 10, height = 12)
