source(here('R', 'librerie.R'))

library(Biostrings)


#COI sequences

# COI sequences -----------------------------------------------------------


fasta_file <- readDNAStringSet("dati/COI riunite 35Rr24Rn8Mm.fas")

# Rimuovere " FR" da tutti gli header
#new_headers <- gsub(" FR", "", names(fasta_file))

new_headers <- paste0("Seq", seq_along(names(fasta_file)), "_", names(fasta_file), "/2023")
short_headers <- paste0("Seq", seq_along(names(fasta_file)))
# Assegnare i nuovi header al DNAStringSet
names(fasta_file) <- new_headers

# Salvare il file FASTA modificato
writeXStringSet(fasta_file, "dati/COI ordinate nomi completi.fas")

#nomi con solo progressivi
names(fasta_file) <- short_headers
writeXStringSet(fasta_file, "dati/COI ordinate nomi corti.fas")

write(names(fasta_file), "COI headers.txt")



# VKORC1 sequences --------------------------------------------------------


fasta_file <- readDNAStringSet("dati/VKORC1 riunite 35Rr24Rn8Mm.fas")

# Rimuovere " RECONSTRUCTED SEQUENCE" da tutti gli header
#new_headers <- gsub(" RECONSTRUCTED SEQUENCE", "", names(fasta_file))
new_headers <- paste0("Seq", seq_along(names(fasta_file)), "_", names(fasta_file), "/2023")
short_headers <- paste0("Seq", seq_along(names(fasta_file)))
#new_headers <- paste0("Seq", seq_along(new_headers), "_", new_headers, "/2023")

# Assegnare i nuovi header al DNAStringSet
names(fasta_file) <- new_headers

# Salvare il file FASTA modificato
writeXStringSet(fasta_file, "dati/VKORC1 ordinate nomi completi.fas")

#nomi con solo progressivi
names(fasta_file) <- short_headers
writeXStringSet(fasta_file, "dati/VKORC1 ordinate nomi corti.fas")

write(names(fasta_file), "VKORC1 headers.txt")


# Rimuovere nomi e tenere solo sequenziali --------------------------------

fasta_file <- readDNAStringSet("dati/VKORC1 sequences nomi corretti.fas")
fasta_file2 <- readDNAStringSet("dati/COI sequences nomi corretti 2.fas")

short_headers <- paste0("Seq", seq_along(names(fasta_file)))
short_headers2 <- paste0("Seq", seq_along(names(fasta_file2)))
short_headers2



# Identificare mancante ---------------------------------------------------

VKORC_file <- readDNAStringSet("dati/VKORC1 sequences.fas")

new_headers_vkorc <- gsub(" RECONSTRUCTED SEQUENCE", "", names(VKORC_file))
names(VKORC_file) <- new_headers_vkorc

COI_file <- readDNAStringSet("dati/COI sequences.fas")

new_headers_COI <- gsub(" FR", "", names(COI_file))
names(COI_file) <- new_headers_COI

COI_headers <- names(COI_file)
VKORC_headers <- names(VKORC_file)

# Identificare la sequenza mancante in COI_file
missing_in_COI <- setdiff(VKORC_headers, COI_headers)

# Visualizzare la sequenza mancante
if(length(missing_in_COI) > 0) {
  cat("La sequenza mancante in COI_file ha l'header:", missing_in_COI, "\n")
} else {
  cat("Nessuna sequenza mancante in COI_file.\n")
}

# Se desideri allineare le sequenze che coincidono
common_headers <- intersect(COI_headers, VKORC_headers)

# Creare una nuova lista di sequenze allineate
aligned_COI <- COI_file[common_headers]
aligned_VKORC <- VKORC_file[common_headers]

# Salvare le sequenze allineate (opzionale)
writeXStringSet(aligned_COI, "dati/COI_sequences_ordinate.fasta")
writeXStringSet(aligned_VKORC, "dati/VKORC1_sequences_ordinate.fasta")



# Verifica ordine headers -------------------------------------------------

# Leggere i file FASTA
COI_sequences <- readDNAStringSet("dati/COI ordinate nomi completi.fas")
VKORC_sequences <- readDNAStringSet("dati/VKORC1 ordinate nomi completi.fas")

# Estrarre gli header (nomi) delle sequenze
COI_headers <- names(COI_sequences)
VKORC_headers <- names(VKORC_sequences)

# Confrontare l'ordine degli header
if (identical(COI_headers, VKORC_headers)) {
  cat("Gli header sono nello stesso ordine in entrambi i file.\n")
} else {
  cat("Gli header NON sono nello stesso ordine nei due file.\n")
  
  # Opzionale: mostrare quali header sono diversi
  different_indices <- which(COI_headers != VKORC_headers)
  cat("Differenze trovate agli indici:\n", different_indices, "\n")
  cat("COI_file headers:", COI_headers[different_indices], "\n")
  cat("VKORC_file headers:", VKORC_headers[different_indices], "\n")
}



# Prova joining esoni con Biostrings --------------------------------------

# Carica il pacchetto necessario
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

library(Biostrings)

# Leggi i file FASTA
exon1 <- readDNAStringSet(here("dati", "EX1.fas"))
exon2 <- readDNAStringSet(here("dati", "EX2.fas"))
exon3 <- readDNAStringSet(here("dati", "EX3.fas"))

# Assicurati che gli header siano nello stesso ordine
if (!identical(names(exon1), names(exon2)) || !identical(names(exon2), names(exon3))) {
  stop("Gli header dei file FASTA non corrispondono o non sono nello stesso ordine.")
}

# Unisci le sequenze corrispondenti riga per riga
combined_sequences <- DNAStringSet(paste0(as.character(exon1), as.character(exon2), as.character(exon3)))

# Mantieni gli header originali
names(combined_sequences) <- names(exon1)

# Salva il file FASTA combinato
writeXStringSet(combined_sequences, filepath = here("complete sequences2.fas"))

