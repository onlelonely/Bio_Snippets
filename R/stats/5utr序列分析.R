# ---------------------------------------------
# Title: 5'UTR序列分析
# Description: From: Source/1. Atlas/🧬 Genomics & Molecular Bio/Gene Regulation/5'UTR序列分析.md (3 blocks)
# ---------------------------------------------

# --- Part 1 ---
# Load necessary library
library(Biostrings)
library(parallel)

# Read the FASTA file
fasta_sequences <- readDNAStringSet("Results/TSS_annotated_TS_sequences.fasta")

# Get sequence names
seq_names <- names(fasta_sequences)

# Function to process one sequence
process_sequence <- function(seq_idx) {
    # Initialize vector for this transcript
    g_content <- numeric(40)
    positions <- seq(1, 200, by=5)
    
    # Get current sequence
    curr_seq <- fasta_sequences[seq_idx]
    
    # Process each 5bp window
    for (i in 1:40) {
        start_pos <- positions[i]
        window_seq <- subseq(curr_seq, start=start_pos, width=5)
        
        # Calculate G content for this window
        g_counts <- vcountPattern("G", window_seq)
        g_content[i] <- g_counts/5  # Divide by window size (5)
    }
    
    # Create data frame for this transcript
    result_df <- data.frame(
        transcript = seq_names[seq_idx],
        setNames(as.list(g_content), 
                paste0("g_content_", 
                      seq(1, 200, by=5)[1:40], "_",
                      seq(5, 200, by=5)[1:40]))
    )
    
    return(result_df)
}
# Set up parallel processing with 4 cores
num_cores <- 4
cl <- makeCluster(num_cores)

# Export required objects to the cluster
clusterExport(cl, c("fasta_sequences", "seq_names"))
clusterEvalQ(cl, library(Biostrings))

# Run parallel processing
all_results <- parLapply(cl, 1:length(fasta_sequences), process_sequence)

# Stop the cluster
stopCluster(cl)

# Combine all results
final_df <- do.call(rbind, all_results)

# Write to CSV
write.csv(final_df, "g_content_5bp_windows_by_transcript.csv", row.names=FALSE)

# --- Part 2 ---
# Define the processing function
process_sequence <- function(seq_idx) {
    # Initialize vectors to store base content for each position
    positions <- seq(1, 200, by=25)  # Calculate from sequence start
    
    # Get current sequence
    curr_seq <- fasta_sequences[seq_idx]
    
    # Create vectors to store content for each window
    window_contents <- lapply(positions, function(start_pos) {
        window_seq <- subseq(curr_seq, start=start_pos, width=25)
        
        # Calculate content for each base
        a_counts <- vcountPattern("A", window_seq)
        t_counts <- vcountPattern("T", window_seq)
        g_counts <- vcountPattern("G", window_seq)
        c_counts <- vcountPattern("C", window_seq)
        
        # Return contents as proportions
        c(a_counts/25, t_counts/25, g_counts/25, c_counts/25)
    })
    
    # Convert to numeric vectors
    contents <- do.call(rbind, window_contents)
    
    # Create data frame for this transcript with one row
    result_df <- data.frame(
        transcript = seq_names[seq_idx],
        setNames(as.list(c(
            contents[,1],  # A content
            contents[,2],  # T content
            contents[,3],  # G content
            contents[,4]   # C content
        )), c(
            paste0("A_content_", positions),
            paste0("T_content_", positions),
            paste0("G_content_", positions),
            paste0("C_content_", positions)
        ))
    )
    
    return(result_df)
}

# Set up parallel processing with 5 cores
num_cores <- 4
cl <- makeCluster(num_cores)

# Export required objects to the cluster
clusterExport(cl, c("fasta_sequences", "seq_names"))
clusterEvalQ(cl, library(Biostrings))

# Run parallel processing
all_results <- parLapply(cl, 1:length(fasta_sequences), process_sequence)

# Stop the cluster
stopCluster(cl)

# Combine all results
final_df <- do.call(rbind, all_results)

# Write to CSV
write.csv(final_df, "Base_content_25bp_windows_by_transcript.csv", row.names=FALSE)

# --- Part 3 ---
# Add these to your library imports at the top
library(rtracklayer)
library(GenomicRanges)

# Download CpG islands from UCSC
cpg_session <- browserSession("UCSC")
genome(cpg_session) <- "hg38"
cpg_query <- ucscTableQuery(cpg_session, table="cpgIslandExt")
cpg_islands <- track(cpg_query)

annotated_tss_data_updated <- read.csv("Results/TSS_annotated_transcripts_plusminus_100bp_with_genes_no_na.csv")

# Convert your TSS data to GRanges
tss_ranges <- GRanges(
    seqnames = annotated_tss_data_updated$seqnames,
    ranges = IRanges(
        start = annotated_tss_data_updated$start,
        end = annotated_tss_data_updated$end
    ),
    strand = annotated_tss_data_updated$strand
)

# Find overlaps between TSS regions and CpG islands
overlaps <- findOverlaps(tss_ranges, cpg_islands)

# Create a new column for CpG status
annotated_tss_data_updated$has_cpg <- FALSE
annotated_tss_data_updated$has_cpg[queryHits(overlaps)] <- TRUE

# Add CpG island details for overlapping regions
cpg_data <- as.data.frame(cpg_islands)[subjectHits(overlaps), ]
annotated_tss_data_updated$cpg_length <- NA
annotated_tss_data_updated$cpg_gc_percent <- NA
annotated_tss_data_updated$cpg_obs_exp_ratio <- NA

annotated_tss_data_updated$cpg_length[queryHits(overlaps)] <- cpg_data$length
annotated_tss_data_updated$cpg_gc_percent[queryHits(overlaps)] <- cpg_data$perCpg
annotated_tss_data_updated$cpg_obs_exp_ratio[queryHits(overlaps)] <- cpg_data$obsExp

# Save the updated data with CpG annotations
write.csv(annotated_tss_data_updated, 
          "TSS_annotated_transcripts_plusminus_100bp_with_genes_and_cpg.csv", 
          row.names = FALSE)