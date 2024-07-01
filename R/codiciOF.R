source(here('R', 'librerie.R'))


# Joining Exons -----------------------------------------------------------


#Extracting headers from fasta

# File path
fasta_file1 <- here("dati", "EX1.fas")

# Read lines from the FASTA file
fasta_lines <- readLines(fasta_file1)

# Extract headers
headers <- fasta_lines[grepl("^>", fasta_lines)]

# Print the extracted headers
headers %>% writeLines(here("headers.txt"))


# extracting sequences from fasta 

# File path

fasta_file2 <- here("dati", "EX2.fas")
fasta_file3 <- here("dati", "EX3.fas")
# Read lines from the FASTA file
fasta_lines_1 <- readLines(fasta_file1)
fasta_lines_2 <- readLines(fasta_file2)
fasta_lines_3 <- readLines(fasta_file3)


# Initialize a list to store concatenated sequences
list1 <- list()
list2 <- list()
list3 <- list()



#EXON 1
# Initialize variables to keep track of the current sequence and header
current_sequence1 <- ""
current_header1 <- ""

# Iterate over the lines of the FASTA file
for (line in fasta_lines_1) {
  if (startsWith(line, ">")) {
    if (current_header1 != "") {
      # Store the previous sequence
      list1 <- c(list1, current_sequence1)
    }
    current_header1 <- line
    current_sequence1 <- ""
  } else {
    current_sequence1 <- paste0(current_sequence1, line)
  }
}

# Store the last sequence
list1 <- c(list1, current_sequence1)


# EXON2
# Initialize variables to keep track of the current sequence and header
current_sequence2 <- ""
current_header2 <- ""

for (line in fasta_lines_2) {
  if (startsWith(line, ">")) {
    if (current_header2 != "") {
      # Store the previous sequence
      list2 <- c(list2, current_sequence2)
    }
    current_header2 <- line
    current_sequence2 <- ""
  } else {
    current_sequence2 <- paste0(current_sequence2, line)
  }
}

# Store the last sequence
list2 <- c(list2, current_sequence2)

# EXON3
# Initialize variables to keep track of the current sequence and header
current_sequence3 <- ""
current_header3 <- ""

for (line in fasta_lines_3) {
  if (startsWith(line, ">")) {
    if (current_header3 != "") {
      # Store the previous sequence
      list3 <- c(list3, current_sequence3)
    }
    current_header3 <- line
    current_sequence3 <- ""
  } else {
    current_sequence3 <- paste0(current_sequence3, line)
  }
}

# Store the last sequence
list3 <- c(list3, current_sequence3)



#combining the 3 exon sequences

# Initialize a list to store concatenated sequences
concatenated_sequences <- list()

# Combine sequences line by line
for (i in 1:length(list1)) {
  concatenated_seq <- paste(list1[i], list2[i], list3[i], sep = "")
  concatenated_sequences <- c(concatenated_sequences, concatenated_seq)
}

#adding headers
VKorc_sequences <- list()
for (i in 1:length(headers)) {
  VKorc_seq <- paste(headers[i], concatenated_sequences[i], sep = "\n")
  VKorc_sequences <- c(VKorc_sequences, VKorc_seq)
}

#exporting complete sequences
VKorc_sequences %>% paste(collapse = "\n") %>% writeLines("Complete sequences.fas")

#Changing headers name

fasta_file4 <- "./Complete sequences.fas"

fasta_seqs <- readDNAStringSet(fasta_file4)

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

ER2 <- ER %>% 
  mutate(catture = if_else(name %in% c("Modena", "Gatteo", "Cesena", "Piacenza", 
                                       "Carpaneto Piacentino", "Bomporto","Crevalcore", "Bagnara di Romagna", "Cervia", 
                                       "Castel San Pietro Terme", "Granarolo dell'Emilia", "Ozzano dell'Emilia", "Russi" ), "si", "no"))
ER3 <- ER %>% 
  mutate(catture = if_else(name == "Cesena", 16, if_else(name == "Modena", 9, if_else(name %in% c("Bomporto", "Gatteo", "Crevalcore"), 3, if_else(name %in% c("Piacenza", "Carpaneto Piacentino", "Bagnara di Romagna", "Cervia", "Castel San Pietro Terme", "Forlì", "Granarolo dell'Emilia", "Lugo", "Ozzano dell'Emilia", "Russi"), 1,0))))) #%>% view()


map1 <- tm_shape(ER) +
  tm_polygons("prov_name", fill.scale = tm_scale_continuous(values= "greys", midpoint = 28000))

map2 <- tm_shape(ER2) +
  tm_polygons("catture", fill.scale = tm_scale_categorical(values = "grays", values.range = c(0.1,0.7)))

#prova 3

province_colors <- scale_fill_manual(values = grey.colors(n = length(unique(ER$prov_name))))

map3 <- tm_shape(ER)+
  tm_polygons("prov_name", fill.scale = tm_scale_categorical(values = "grays", values.range = c(0.1,0.75)), fill.legend = tm_legend(title = "Provinces")) +
  tm_borders() +
  tm_shape(ER %>% filter(name %in% c("Modena", "Gatteo", "Cesena", "Piacenza", 
                                     "Carpaneto Piacentino", "Bomporto","Crevalcore", "Bagnara di Romagna", "Cervia", 
                                     "Castel San Pietro Terme", "Granarolo dell'Emilia", "Ozzano dell'Emilia", "Russi" ))) +
  tm_polygons(fill = "red", fill.legend = tm_legend(title = "Comuni catture", show = T, position = "bottom")) +
  tm_title_out("Sampling Map", position = tm_pos_out("center", "top"))


tm_shape(ER3)+
  tm_polygons("prov_name", fill.scale = tm_scale_categorical(values = "grays", values.range = c(0.1,0.75)), fill.legend = tm_legend(title = "Provinces")) +
  tm_borders() +
  tm_shape(ER3 %>% filter(name %in% c("Modena", "Gatteo", "Cesena", "Piacenza", 
                                     "Carpaneto Piacentino", "Bomporto","Crevalcore", "Bagnara di Romagna", "Cervia", 
                                     "Castel San Pietro Terme", "Granarolo dell'Emilia", "Ozzano dell'Emilia", "Russi" ))) +
  tm_polygons("catture", fill.scale = tm_scale_continuous(values = "reds", values.range = c(0.2,1)), fill.legend = tm_legend(title = "Samples", format = list(1,3,9,13))) +
  tm_title_out("Sampling Map", position = tm_pos_out("center", "top"))






#mappa con legenda modificata

# Calcola il punto medio di ciascuna provincia

# Convert dataframe in sf object
ER3_sf <- st_as_sf(ER3)

# Group by province and summarize municipalities' geometry by province
province_geom <- ER3_sf %>%
  group_by(prov_acr) %>%
  summarize(geometry = st_union(geometry))

# Calculate province centroids
province_centroids <- province_geom %>%
  st_centroid()

map4 <- tm_shape(ER3) +
  tm_polygons("prov_name", fill.scale = tm_scale_categorical(values = "grays", values.range = c(0.1, 0.75)), fill.legend = tm_legend_hide()) +
  tm_borders() +
  tm_shape(ER3 %>% filter(name %in% c("Modena", "Gatteo", "Cesena", "Piacenza", 
                                      "Carpaneto Piacentino", "Bomporto", "Crevalcore", "Bagnara di Romagna", "Cervia", 
                                      "Castel San Pietro Terme", "Granarolo dell'Emilia", "Forlì", "Lugo", "Ozzano dell'Emilia", "Russi" ))) +
  tm_polygons("catture", fill.scale = tm_scale_categorical(values = "reds", values.range = c(0.2,1)), 
              fill.legend = tm_legend(title = "N° of samples")) +
  tm_shape(province_centroids) + # Aggiungi le etichette per le province
  tm_text("prov_acr", size = 1, col = "black", fontface = "bold") +
  tm_title_out("Sampling Map", position = tm_pos_out("center", "top"))

tmap_save(map4, "Mappa_catture.png")






