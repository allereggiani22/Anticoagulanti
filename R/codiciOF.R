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
#cat("Extracted headers:\n")
#cat(headers, sep = "\n") 
headers %>% writeLines(here("headers.txt"))


# extracting sequences from fasta -----------------------------------------

# File path
#fasta_file <- "path_to_fasta_file.fasta"
fasta_file2 <- here("dati", "ALLINEAMENTO RATTI PROGETTO EX2.fas")
fasta_file3 <- here("dati", "ALLINEAMENTO RATTI PROGETTO EX3.fas")
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

#così funziona, rimane il problema del raddoppio del codone di giunzione



























# Initialize a variable to store the concatenated sequence
sequence1 <- ""
sequence2 <- ""
sequence3 <- ""

# Iterate over the lines of the FASTA file
for (line in fasta_lines_1) {
  if (!startsWith(line, ">")) {
    sequence1 <- paste0(sequence1, line, "\n")
  }
}

for (line in fasta_lines_2) {
  if (!startsWith(line, ">")) {
    sequence2 <- paste0(sequence2, line, "\n")
  }
}

for (line in fasta_lines_3) {
  if (!startsWith(line, ">")) {
    sequence3 <- paste0(sequence3, line, "\n")
  }
}

#now I have to split the txt sequences into vectors

list1 <- sequence1 %>% strsplit("\n")
list2 <- sequence2 %>% strsplit("\n")
list3 <- sequence3 %>% strsplit("\n")

#Initialize a list to store concatenated sequences
concatenated_sequences <- list()

# Initialize variables to keep track of the current sequence
current_sequence <- ""

for (i in 1:length(sequences1)) {
  if (!startsWith(sequences1[i], ">")) {
    current_sequence <- paste0(current_sequence, sequences1[i],"\n")
    current_sequence <- paste0(current_sequence, sequences2[i])
    current_sequence <- paste0(current_sequence, sequences3[i])
    
    concatenated_sequences <- c(concatenated_sequences, current_sequence)
    
    current_sequence <- ""  # Reset for the next sequence
  }
}

# Save the concatenated sequences to a new FASTA file
seqs <- as.character(concatenated_sequences)
writeLines(seqs, "concatenated_sequences.fasta")


