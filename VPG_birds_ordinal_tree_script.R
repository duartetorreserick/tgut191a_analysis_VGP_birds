# --- set working directory safely ---
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  this_path <- rstudioapi::getSourceEditorContext()$path
  if (!is.null(this_path) && nzchar(this_path)) {
    setwd(dirname(this_path))
  }
} else if (!interactive()) {
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(script_path) > 0) {
    setwd(dirname(normalizePath(script_path)))
  }
}

enrichment_data <- read_tsv("/Users/jarvislab/Documents/VGP/Zebrafinch/annotation_density/enrichment_final.tsv")

# Package installation and loading
install_load <- function(packages){
  for (p in packages) {
    if (p %in% rownames(installed.packages())) {
      library(p, character.only=TRUE)
    } else {
      install.packages(p, repos = "http://cran.us.r-project.org")
      library(p,character.only = TRUE)
    }
  }
}

required_packages <- c("readxl", "data.table", "tidyverse", "ape", "ggplot2", 
                       "TreeTools", "devtools", "stringr", "ggnewscale", 
                       "BiocManager", "RColorBrewer", "googlesheets4", "ggforce", 
                       "gtable", "ggpubr", "grid", "remotes")

install_load(required_packages)

# Install Bioconductor packages
BiocManager_packages <- c("treeio", "ggtree", "ggtreeExtra", "tidytree")
for (p in BiocManager_packages) {
  if(!require(p, quietly=TRUE, character.only=TRUE)){
    BiocManager::install(p)
    library(p, character.only=TRUE)
  }
}

# Install GitHub packages
remotes::install_github("ropensci/bold", quiet=TRUE)
remotes::install_github("ropensci/taxize", quiet=TRUE)
remotes::install_github("ropensci/rotl", quiet=TRUE)
devtools::install_github("phylotastic/datelife", quiet=TRUE)
devtools::install_github("phylotastic/datelifeplot", quiet=TRUE)

library(rotl)

# Load VGP ordinal list
gs4_auth(token = NULL, scopes = "https://www.googleapis.com/auth/spreadsheets.readonly", 
         email = "duartetorreserick@gmail.com")
ordinal_list <- read_sheet(ss = "17aOjpVgclwdDcDccx7Jpuy4IL6PUcVYntB-iwfj0aZ4", sheet = 1)

# Filter and process ordinal list
ordinal_list <- ordinal_list %>%
  dplyr::filter(Lineage %in% c("Birds", "Aves")) %>%
  as.data.table()

setnames(ordinal_list, "Orders Scientific Name (inferred >50 MYA divergence times)", "order_50MYA")
ordinal_list <- ordinal_list %>% 
  mutate(order_50MYA = gsub(")", "", order_50MYA)) %>%
  separate(order_50MYA, into = c("order", "suborder1"), sep = " \\(| > |\\|", extra = "drop")

# Name corrections
ordinal_list$`Scientific Name updated` <- ordinal_list$`Scientific Name`
ordinal_list$`Scientific Name updated`[ordinal_list$`Scientific Name` == "Ammospiza nelsoni"] <- "Ammodramus nelsoni"

# Filter for complete records
ordinal_list <- ordinal_list %>%
  filter(!is.na(`UCSC Browser main haplotype`) & 
           `UCSC Browser main haplotype` != "" & 
           str_trim(`UCSC Browser main haplotype`) != "")

# Select relevant columns
selected_columns <- ordinal_list %>%
  select(Order = order, 
         `Scientific Name updated`,
         `Assembly ID`,
         `UCSC Browser main haplotype`)

# Join with repeat hits data
repeat_hits <- read_tsv("/Users/jarvislab/Documents/VGP/Zebrafinch/tree/repeat_hits_genome.tsv",
                        col_names = c("Assembly ID", "count"))

selected_columns <- selected_columns %>%
  left_join(repeat_hits, by = "Assembly ID") %>%
  mutate(count = replace_na(as.integer(count), 0))

# Join with N50 data
n50 <- read_tsv("/Users/jarvislab/Documents/VGP/Zebrafinch/treeoflife/ordinal_list_with_n50.tsv",
                col_names = c("order", "Scientific Name", "Assembly ID", "Accession", "N50"))

selected_columns <- selected_columns %>%
  left_join(n50, by = "Assembly ID")

# Apply name corrections
name_corrections <- data.frame(
  original_name = c("Ammospiza maritima", "Oenanthe melanoleuca", "Coloeus monedula",
                    "Strigops habroptilus", "Amazona ochrocephala", "Guaruba guaruba",
                    "Dryobates pubescens", "Alca Torda", "Aegotheles albertisi",
                    "Apteryx mantelli", "Rhea pennata"),
  correct_name = c("Ammodramus maritimus", "Muscicapa mesoleuca", "Corvus monedula",
                   "Strigops habroptila", "Amazona ochrocephala", "Guaruba guarouba",
                   "Picoides pubescens", "Alca torda", "Aegotheles albertisi",
                   "Apteryx australis mantelli", "Pterocnemia pennata"),
  stringsAsFactors = FALSE
)

selected_columns_corrected <- selected_columns
for(i in 1:nrow(name_corrections)) {
  old_name <- name_corrections$original_name[i]
  new_name <- name_corrections$correct_name[i]
  selected_columns_corrected$`Scientific Name updated`[selected_columns_corrected$`Scientific Name updated` == old_name] <- new_name
}

# Create phylogenetic tree
VGP_ordinal_resolved_names <- rotl::tnrs_match_names(names = selected_columns_corrected$`Scientific Name updated`)
in_tree <- rotl::is_in_tree(VGP_ordinal_resolved_names$ott_id)
VGP_ordinal_subtree <- rotl::tol_induced_subtree(VGP_ordinal_resolved_names$ott_id[in_tree])

# Clean tip labels
VGP_ordinal_subtree$tip.label <- gsub("_ott\\d+", "", VGP_ordinal_subtree$tip.label)
VGP_ordinal_subtree$tip.label <- gsub("_", " ", VGP_ordinal_subtree$tip.label)

# Color generation function
make_distinct_hcl <- function(n, chroma = 70, lum1 = 60, lum2 = 80, hue_start = 10) {
  h <- seq(hue_start, hue_start + 360, length.out = n + 1)[1:n] %% 360
  l <- rep(c(lum1, lum2), length.out = n)
  grDevices::hcl(h = h, c = chroma, l = l)
}

# Get tip labels and prepare data
temp_tree <- ggtree(VGP_ordinal_subtree, ladderize = TRUE)
tip_labels <- temp_tree$data %>% filter(isTip) %>% arrange(y) %>% pull(label)

# Prepare barplot data for repeats
barplot_data <- tibble(label = tip_labels) %>%
  left_join(selected_columns_corrected, by = c("label" = "Scientific Name updated")) %>%
  mutate(count = replace_na(count, 0),
         count_log = log10(count + 1),
         Order = replace_na(Order, "Unknown"))

# Prepare N50 data
barplot_n50_data <- tibble(label = tip_labels) %>%
  left_join(selected_columns_corrected, by = c("label" = "Scientific Name updated")) %>%
  mutate(N50 = replace_na(N50, 0),
         N50_log = log10(N50 + 1),
         Order = replace_na(Order, "Unknown"))

# Create order data for tree annotation
ord_data <- barplot_data %>%
  select(label, taxon_order = Order)

# Generate colors
tree_order_sequence <- barplot_data %>%
  arrange(match(label, tip_labels)) %>%
  select(Order) %>%
  distinct() %>%
  filter(Order != "Unknown") %>%
  pull(Order)

order_colors_tree_ordered <- make_distinct_hcl(length(tree_order_sequence))
names(order_colors_tree_ordered) <- tree_order_sequence

present_orders <- unique(barplot_data$Order)
present_orders <- present_orders[!is.na(present_orders) & present_orders != "Unknown"]
present_orders <- rev(present_orders)

# Create tree with both barplots
p_tree_page_size <- ggtree(VGP_ordinal_subtree, ladderize = TRUE, branch.length = "branch.length") %<+% ord_data +
  theme_tree2() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  geom_tiplab(size = 2, hjust = -0.01, fontface = "italic",
              aes(color = taxon_order), show.legend = FALSE) +
  scale_color_manual(values = order_colors_tree_ordered) +
  geom_treescale(x = 0, y = 1, width = 0.005, fontsize = 2) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5, 120, 5, 5)) +
  # First barplot: Repeat count
  geom_fruit(data = barplot_data,
             geom = geom_col,
             mapping = aes(y = label, x = count_log, fill = Order),
             orientation = "y",
             offset = 1.2,
             pwidth = 0.8,
             axis.params = list(axis = "x",
                                breaks = c(0, 2, 4, 6,8),
                                labels = c("1", "100", "10K", "1M","1B"),
                                text.size = 2,
                                title = "Repeat Count")) +
  scale_fill_manual(values = order_colors_tree_ordered, name = "Order", breaks = present_orders, guide = "none") 
  # new_scale_fill() +
  # # Second barplot: N50 values
  # geom_fruit(data = barplot_n50_data,
  #            geom = geom_col,
  #            mapping = aes(y = label, x = N50_log, fill = Order),
  #            orientation = "y",
  #            offset = 2,
  #            pwidth = 0.8,
  #            axis.params = list(axis = "x",
  #                               breaks = c(0, 2, 4, 6, 8,10),
  #                               labels = c("0", "20M", "40M", "60M", "80M","100M"),
  #                               text.size = 2,
  #                               title = "Contig N50 (bp)")) +
  # scale_fill_manual(values = order_colors_tree_ordered, name = "Order", breaks = present_orders, guide = "none")

print(p_tree_page_size)

# Save the plot
ggplot2::ggsave(filename = "tree_with_n50_barplot.svg",
                plot = p_tree_page_size,
                device = svglite::svglite,
                width = 8.5, height = 11,
                units = "in", bg = "white")


# Define presence thresholds
THRESH_PRESENT <- 1       # ≥1 hit = detected
THRESH_ROBUST  <- 10000     # ≥100 hits = robustly present (tune as needed)

df <- selected_columns_corrected %>%
  select(assembly = `Assembly ID`, Order, count) %>%
  mutate(count = ifelse(is.na(count), 0L, as.integer(count)),
         present = count >= THRESH_PRESENT,
         robust  = count >= THRESH_ROBUST)

YY_total            <- n_distinct(df$assembly)
XX_present          <- df %>% filter(present) %>% summarise(n = n_distinct(assembly)) %>% pull(n)
XX_robust           <- df %>% filter(robust)  %>% summarise(n = n_distinct(assembly)) %>% pull(n)
XX_robust_passerine <- df %>% filter(robust, Order == "Passeriformes") %>% summarise(n = n_distinct(assembly)) %>% pull(n)

pct_present <- sprintf("%.1f%%", 100 * XX_present / YY_total)
pct_robust  <- sprintf("%.1f%%", 100 * XX_robust  / YY_total)
pct_pass    <- if (XX_robust > 0) sprintf("%.1f%%", 100 * XX_robust_passerine / XX_robust) else "0.0%"
list(YY_total=YY_total, XX_present=XX_present, XX_robust=XX_robust, XX_robust_passerine=XX_robust_passerine,
     pct_present=pct_present, pct_robust=pct_robust, pct_pass=pct_pass)

