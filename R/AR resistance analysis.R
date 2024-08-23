source(here('R', 'librerie.R'))


# Joining Exons -----------------------------------------------------------

# Load necessary library
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# 
# library(Biostrings)

# Read FASTA files
exon1 <- readDNAStringSet(here("dati", "EX1.fas"))
exon2 <- readDNAStringSet(here("dati", "EX2.fas"))
exon3 <- readDNAStringSet(here("dati", "EX3.fas"))


# Join sequences line by line
combined_sequences <- DNAStringSet(paste0(as.character(exon1), as.character(exon2), as.character(exon3)))

# Keep headers from exon 1
names(combined_sequences) <- names(exon1)

# Salva il file FASTA combinato
writeXStringSet(combined_sequences, filepath = here("complete sequences.fas"))
#Changing headers name

fasta_file <- "./Complete sequences.fas"

fasta_seqs <- readDNAStringSet(fasta_file)

#function to convert in tibble

fasta_df <-  tibble(
  name = names(fasta_seqs),
  sequence = as.character(fasta_seqs)
)

#modify header names

fasta_df <- fasta_df %>% 
  mutate(name = str_replace(name, "EX1.*", "RECONSTRUCTED SEQUENCE"))

#Convert back to DNAStringSet

renamed_seqs <- DNAStringSet(fasta_df$sequence)
names(renamed_seqs) <- fasta_df$name

renamed_file <- "Complete sequences renamed.fas"
writeXStringSet(renamed_seqs, renamed_file)




# Mapping -----------------------------------------------------------------

library(tmap)    # for static and interactive maps


ER <- st_read("dati/limits_R_8_municipalities.geojson")

ER3 <- ER %>% 
  mutate(catture = if_else(name == "Cesena", 16, if_else(name == "Modena", 9, if_else(name %in% c("Bomporto", "Gatteo", "Crevalcore"), 3, if_else(name %in% c("Piacenza", "Carpaneto Piacentino", "Bagnara di Romagna", "Cervia", "Castel San Pietro Terme", "Forlì", "Granarolo dell'Emilia", "Lugo", "Ozzano dell'Emilia", "Russi"), 1,0))))) #%>% view()



#identifying provinces by color

province_colors <- scale_fill_manual(values = grey.colors(n = length(unique(ER$prov_name))))

# Calculate center of each province

# Convert dataframe in sf object
ER3_sf <- st_as_sf(ER3)

# Group by province and summarize municipalities' geometry by province
province_geom <- ER3_sf %>%
  group_by(prov_acr) %>%
  summarize(geometry = st_union(geometry))

# Calculate province centroids
province_centroids <- province_geom %>%
  st_centroid()


# Calculate frequencies of values in "catture"

ER3_filtered <-as.data.frame(ER3) %>% 
  filter(name %in% c("Modena", "Gatteo", "Cesena", "Piacenza", 
                     "Carpaneto Piacentino", "Bomporto", "Crevalcore", "Bagnara di Romagna", "Cervia", 
                     "Castel San Pietro Terme", "Granarolo dell'Emilia", "Forlì", "Lugo", "Ozzano dell'Emilia", "Russi"))

# Calcolare le frequencies
freq_table <- ER3_filtered %>%
  count(catture) %>%
  dplyr::rename(frequency = n)

# unite frequencies with inital df
ER3_filtered <- ER3_filtered %>%
  left_join(freq_table, by = "catture")

# create new combined colums for frequency
ER3_filtered <- ER3_filtered %>%
  mutate(catture_with_freq = paste(catture, " (n=", frequency, ")", sep = ""))

# Order levels for "catture" and duplicates removal
unique_levels <- ER3_filtered %>%
  distinct(catture, .keep_all = TRUE) %>%
  arrange(catture) %>%
  pull(catture_with_freq)

# Convert catture_with_freq in a factor with ordered levels
ER3_filtered$catture_with_freq <- factor(ER3_filtered$catture_with_freq, 
                                         levels = unique_levels)


# Creating the map
map <- tm_shape(ER3) +
  tm_polygons("prov_name", fill.scale = tm_scale_categorical(values = "grays", values.range = c(0.1, 0.75)), fill.legend = tm_legend_hide()) +
  tm_borders() +
  tm_shape(st_as_sf(ER3_filtered)) +
  tm_polygons("catture_with_freq", fill.scale = tm_scale_categorical(values = "reds", values.range = c(0.2, 1)), 
              fill.legend = tm_legend(title = "N° of samples (municipalities)")) +
  tm_shape(province_centroids) + # Aggiungi le etichette per le province
  tm_text("prov_acr", size = 1, col = "black", fontface = "bold") +
  tm_title_out("Sampling Map", position = tm_pos_out("center", "top")) +
  tm_layout(legend.position = c("right", "top"))

tmap_save(map, "Mappa_catture_rev-1.png")
