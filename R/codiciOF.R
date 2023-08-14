source(here('R', 'librerie.R'))

#Voglio unire le sequenze di 3 diverse sequenze di 3 esoni per ricostruire la CDS del gene VKORC1


# File paths
file_path_1 <- here("dati","Ex1.txt")
file_path_2 <- here("dati","Ex2.txt")
file_path_3 <- here("dati","Ex3.txt")


# Read sequences from the FASTA files
sequences_1 <- readLines(file_path_1)
sequences_2 <- readLines(file_path_2)
sequences_3 <- readLines(file_path_3)

# Initialize a list to store concatenated sequences
concatenated_sequences <- list()

# Initialize variables to keep track of the current sequence
current_sequence <- ""

for (i in 1:length(sequences_1)) {
  if (!startsWith(sequences_1[i], ">")) {
    current_sequence <- paste0(current_sequence, sequences_1[i])
    current_sequence <- paste0(current_sequence, sequences_2[i])
    current_sequence <- paste0(current_sequence, sequences_3[i])
    
    concatenated_sequences <- c(concatenated_sequences, current_sequence)
    
    current_sequence <- ""  # Reset for the next sequence
  }
}

# Save the concatenated sequences to a new FASTA file
seqs <- as.character(concatenated_sequences)
writeLines(seqs, "concatenated_sequences.fasta")



# Extracting headers from fasta -------------------------------------------


# File path
fasta_file1 <- here("dati", "ALLINEAMENTO RATTI PROGETTO EX1.fas")

# Read lines from the FASTA file
fasta_lines <- readLines(fasta_file1)

# Extract headers
headers <- fasta_lines[grepl("^>", fasta_lines)]

# Print the extracted headers
cat("Extracted headers:\n")
cat(headers, sep = "\n") %>% as.character() %>% writeLines(here("headers.txt"))


# extracting sequences from fasta -----------------------------------------

# File path
#fasta_file <- "path_to_fasta_file.fasta"
fasta_file2 <- here("dati", "ALLINEAMENTO RATTI PROGETTO EX2.fas")
fasta_file3 <- here("dati", "ALLINEAMENTO RATTI PROGETTO EX3.fas")
# Read lines from the FASTA file
sequences1 <- readLines(fasta_file1)
sequences2 <- readLines(fasta_file2)
sequences3 <- readLines(fasta_file3)

# Initialize a variable to store the concatenated sequence
sequence <- ""

# Iterate over the lines of the FASTA file
for (line in fasta_lines) {
  if (!startsWith(line, ">")) {
    sequence <- paste0(sequence, line, "\n")
  }
}

#Initialize a list to store concatenated sequences
concatenated_sequences <- list()

# Initialize variables to keep track of the current sequence
current_sequence <- ""

for (i in 1:length(sequences_1)) {
  if (!startsWith(sequences_1[i], ">")) {
    current_sequence <- paste0(current_sequence, sequences_1[i])
    current_sequence <- paste0(current_sequence, sequences_2[i])
    current_sequence <- paste0(current_sequence, sequences_3[i])
    
    concatenated_sequences <- c(concatenated_sequences, current_sequence)
    
    current_sequence <- ""  # Reset for the next sequence
  }
}

# Save the concatenated sequences to a new FASTA file
seqs <- as.character(concatenated_sequences)
writeLines(seqs, "concatenated_sequences.fasta")


