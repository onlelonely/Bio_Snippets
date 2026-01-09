# ---------------------------------------------
# Title: 5'UTR序列分析
# Description: From: Source/1. Atlas/🧬 Genomics & Molecular Bio/Gene Regulation/5'UTR序列分析.md (3 blocks)
# ---------------------------------------------

# --- Part 1 ---
# Load necessary libraries
library(CAGEr)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(dplyr)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)


# Load UCSC transcript db
txdb <- TxDb.Hsapiens.UCSC.hg38.refGene


# Define a function to extract TSS ±100 bp for all transcripts
get_tss_with_annotation_for_transcripts <- function(txdb, upstream = 100, downstream = 100) {
  # Get transcript regions and define their TSS (promoters)
  standardChrs <- paste0("chr", c(1:22, "X", "Y"))
  transcriptRegions <- transcripts(txdb)
  transcriptRegions <- transcriptRegions[seqnames(transcriptRegions) %in% standardChrs]
  tssRegions <- promoters(transcriptRegions, upstream = upstream, downstream = downstream)
  
  # Convert TSS regions to a data frame with transcript IDs
  tssRegionsDF <- as.data.frame(tssRegions)
  return(tssRegionsDF)
}

# Call the function to extract TSS ±100 bp for all transcripts
annotated_tss_data <- get_tss_with_annotation_for_transcripts(txdb)

# Save the data
write.csv(annotated_tss_data, "TSS_annotated_transcripts_plusminus_100bp.csv", row.names = FALSE)

# Load Biomart library
library(biomaRt)

# Get human genes from Ensembl
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Extract UCSC transcript IDs from annotated_tss_data
refseq_ids <- annotated_tss_data$tx_name

# Query biomaRt for mRNA IDs - split into two queries
gene_info_mrna1 <- getBM(
  attributes = c(
    "refseq_mrna",
    "hgnc_symbol"
  ),
  filters = "refseq_mrna",
  values = refseq_ids,
  mart = ensembl
)

gene_info_mrna2 <- getBM(
  attributes = c(
    "refseq_mrna",
    "entrezgene_id"
  ),
  filters = "refseq_mrna",
  values = refseq_ids,
  mart = ensembl
)

# Query biomaRt for ncRNA IDs - split into two queries
gene_info_nc1 <- getBM(
  attributes = c(
    "refseq_ncrna",
    "hgnc_symbol"
  ),
  filters = "refseq_ncrna",
  values = refseq_ids,
  mart = ensembl
)

gene_info_nc2 <- getBM(
  attributes = c(
    "refseq_ncrna",
    "entrezgene_id"
  ),
  filters = "refseq_ncrna",
  values = refseq_ids,
  mart = ensembl
)

# Merge mRNA results
gene_info_mrna <- merge(
  gene_info_mrna1,
  gene_info_mrna2,
  by = "refseq_mrna"
)

# Merge ncRNA results
gene_info_nc <- merge(
  gene_info_nc1,
  gene_info_nc2,
  by = "refseq_ncrna"
)

# Combine the results
gene_info_combined <- rbind(
  gene_info_mrna %>%
    rename(refseq_id = refseq_mrna) %>%
    filter(refseq_id != ""),
  gene_info_nc %>%
    rename(refseq_id = refseq_ncrna) %>%
    filter(refseq_id != "")
) %>% distinct()

# Merge the gene information with annotated_tss_data
annotated_tss_data_updated <- merge(
  annotated_tss_data,
  gene_info_combined,
  by.x = "tx_name",
  by.y = "refseq_id",
  all.x = TRUE
)

# Save the updated data
write.csv(annotated_tss_data_updated, "TSS_annotated_transcripts_plusminus_100bp_with_genes.csv", row.names = FALSE)

# remove entrezgene_id = NA
annotated_tss_data_updated <- annotated_tss_data_updated[!is.na(annotated_tss_data_updated$entrezgene_id), ]

# Save the updated data
write.csv(annotated_tss_data_updated, "TSS_annotated_transcripts_plusminus_100bp_with_genes_no_na.csv", row.names = FALSE)

# Load necessary libraries
library(Biostrings)
library(GenomicRanges)

# Function to extract sequences and save as FASTA
extract_sequences_to_fasta <- function(annotated_tss_data, fasta_output, genome = BSgenome.Hsapiens.UCSC.hg38) {
  # Convert the annotated data to a GRanges object
  gr <- GRanges(
    seqnames = annotated_tss_data$seqnames,
    ranges = IRanges(start = annotated_tss_data$start, end = annotated_tss_data$end),
    strand = annotated_tss_data$strand
  )
  
  # Extract sequences from the genome
  sequences <- getSeq(genome, gr)
  
# Create FASTA headers with better NA handling
  fasta_headers <- sapply(1:nrow(annotated_tss_data_updated), function(i) {
    symbol <- if (!is.na(annotated_tss_data_updated$hgnc_symbol[i])) 
      annotated_tss_data_updated$hgnc_symbol[i] else "NA"
    tx_name <- if (!is.na(annotated_tss_data_updated$tx_name[i])) 
      annotated_tss_data_updated$tx_name[i] else "NA"
    entrez <- if (!is.na(annotated_tss_data_updated$entrezgene_id[i])) 
      annotated_tss_data_updated$entrezgene_id[i] else "NA"
    paste(symbol, tx_name, entrez, sep = "|")
  })
  
  # Combine headers and sequences into a named DNAStringSet
  fasta <- DNAStringSet(sequences)
  names(fasta) <- fasta_headers
  
  # Save the sequences in FASTA format
  writeXStringSet(fasta, fasta_output)
  
  message("FASTA file saved as: ", fasta_output)
}

# Define output FASTA file
fasta_file <- "TSS_annotated_wholeTS_sequences.fasta"

# Call the function to extract sequences and save as FASTA
extract_sequences_to_fasta(annotated_tss_data_updated, fasta_file)

# --- Part 2 ---
# Load necessary libraries
library(CAGEr)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)


# Load UCSC transcript db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# Define a function to extract TSS ±100 bp for all transcripts
get_tss_with_annotation_for_transcripts <- function(txdb, upstream = 100, downstream = 100) {
  # Get transcript regions and define their TSS (promoters)
  standardChrs <- paste0("chr", c(1:22, "X", "Y"))
  transcriptRegions <- transcripts(txdb)
  transcriptRegions <- transcriptRegions[seqnames(transcriptRegions) %in% standardChrs]
  tssRegions <- promoters(transcriptRegions, upstream = upstream, downstream = downstream)
  
  # Convert TSS regions to a data frame with transcript IDs
  tssRegionsDF <- as.data.frame(tssRegions)
  return(tssRegionsDF)
}

# Call the function to extract TSS ±100 bp for all transcripts
annotated_tss_data <- get_tss_with_annotation_for_transcripts(txdb)

# Save the data
write.csv(annotated_tss_data, "TSS_annotated_transcripts_plusminus_100bp.csv", row.names = FALSE)

# Load biomaRt library
library(biomaRt)

# Get human genes from Ensembl
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Extract UCSC transcript IDs from annotated_tss_data
ucsc_ids <- annotated_tss_data$tx_name

# Query biomaRt for gene information
gene_info <- getBM(
  attributes = c(
    "ensembl_transcript_id_version", 
    "hgnc_symbol",
    "entrezgene_id"
  ),
  filters = "ensembl_transcript_id_version",
  values = ucsc_ids,
  mart = ensembl
)

# Merge the gene information with annotated_tss_data
annotated_tss_data_updated <- merge(
  annotated_tss_data,
  gene_info,
  by.x = "tx_name",
  by.y = "ensembl_transcript_id_version",
  all.x = TRUE
)

# Save the updated data
write.csv(annotated_tss_data_updated, "TSS_annotated_transcripts_plusminus_100bp_with_genes.csv", row.names = FALSE)

# remove entrezgene_id = NA
annotated_tss_data_updated <- annotated_tss_data_updated[!is.na(annotated_tss_data_updated$entrezgene_id), ]
# Save the updated data
write.csv(annotated_tss_data_updated, "TSS_annotated_transcripts_plusminus_100bp_with_genes_no_na.csv", row.names = FALSE)

# Load necessary libraries
library(Biostrings)
library(GenomicRanges)

# Function to extract sequences and save as FASTA
extract_sequences_to_fasta <- function(annotated_tss_data, fasta_output, genome = BSgenome.Hsapiens.UCSC.hg38) {
  # Convert the annotated data to a GRanges object
  gr <- GRanges(
    seqnames = annotated_tss_data$seqnames,
    ranges = IRanges(start = annotated_tss_data$start, end = annotated_tss_data$end),
    strand = annotated_tss_data$strand
  )
  
  # Extract sequences from the genome
  sequences <- getSeq(genome, gr)
  
# Create FASTA headers with better NA handling
  fasta_headers <- sapply(1:nrow(annotated_tss_data_updated), function(i) {
    symbol <- if (!is.na(annotated_tss_data_updated$hgnc_symbol[i])) 
      annotated_tss_data_updated$hgnc_symbol[i] else "NA"
    tx_name <- if (!is.na(annotated_tss_data_updated$tx_name[i])) 
      annotated_tss_data_updated$tx_name[i] else "NA"
    entrez <- if (!is.na(annotated_tss_data_updated$entrezgene_id[i])) 
      annotated_tss_data_updated$entrezgene_id[i] else "NA"
    paste(symbol, tx_name, entrez, sep = "|")
  })
  
  # Combine headers and sequences into a named DNAStringSet
  fasta <- DNAStringSet(sequences)
  names(fasta) <- fasta_headers
  
  # Save the sequences in FASTA format
  writeXStringSet(fasta, fasta_output)
  
  message("FASTA file saved as: ", fasta_output)
}

# Define output FASTA file
fasta_file <- "TSS_annotated_wholeTS_sequences.fasta"

# Call the function to extract sequences and save as FASTA
extract_sequences_to_fasta(annotated_tss_data_updated, fasta_file)

# --- Part 3 ---
# Load necessary libraries
library(biomaRt)
library(dplyr)

# Read the CSV file
transcripts <- read.csv("transcripts.csv")

# Connect to Ensembl using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query Ensembl for Canonical Transcripts
canonical_transcripts <- getBM(
  attributes = c("hgnc_symbol", "refseq_mrna", "transcript_is_canonical", "transcript_mane_select"),
  filters = "hgnc_symbol",
  values = unique(transcripts$GeneSymbol),
  mart = ensembl
)

# Filter only canonical transcripts
canonical_transcripts <- canonical_transcripts %>%
  filter(transcript_is_canonical == 1)

# Merge the canonical transcript info with your input
filtered_transcripts <- transcripts %>%
  inner_join(canonical_transcripts, by = c("GeneSymbol" = "hgnc_symbol", "Transcript" = "refseq_mrna"))

  # Save transcript_mane_select only in a csv, ignoring version and removing duplicates
  filtered_transcripts %>%
    mutate(transcript_mane_select = sub("\\..*", "", transcript_mane_select)) %>%
    distinct(transcript_mane_select) %>%
    filter(transcript_mane_select != "") %>%
    select(transcript_mane_select) %>%
    write.csv("transcript_mane_select.csv", row.names = FALSE)