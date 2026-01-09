# ---------------------------------------------
# Title: Sequence logo plot
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Sequence logo plot.md
# ---------------------------------------------

library(ShortRead)
library(ggseqlogo)
library(ggplot2)

# Directory with your FASTQ files
fastq_dir <- "path/to/your/folder"
files <- list.files(fastq_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Output directory for plots
plot_dir <- "path/to/save/plots"
dir.create(plot_dir, showWarnings = FALSE)

# Parameters
fixed_length <- 30
max_reads <- 100000  # limit to speed up if needed

for (f in files) {
  # Get clean sample name
  sample_name <- tools::file_path_sans_ext(basename(f))
  sample_name <- gsub("\\.fastq$", "", sample_name)  # if double extension
  
  message("Processing: ", sample_name)
  
  fq <- readFastq(f)
  seqs <- as.character(sread(fq))
  
  # Sample to reduce memory use
  if (length(seqs) > max_reads) {
    seqs <- sample(seqs, max_reads)
  }
  
  # Trim and filter
  seqs_trimmed <- substr(seqs, 1, fixed_length)
  seqs_trimmed <- seqs_trimmed[nchar(seqs_trimmed) == fixed_length]
  
  # Skip if too few reads
  if (length(seqs_trimmed) < 100) {
    warning("Too few reads in ", sample_name, " after filtering")
    next
  }
  
  # Generate plot
  p <- ggseqlogo(seqs_trimmed) + ggtitle(sample_name)
  
  # Save
  ggsave(filename = file.path(plot_dir, paste0(sample_name, "_logo.png")),
         plot = p, width = 6, height = 4, dpi = 300)
}