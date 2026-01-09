#!/bin/bash
# ---------------------------------------------
# Title: Sequence random sampling
# Description: From: Source/1. Atlas/🛠️ Tools & Platforms/Sequence random sampling.md
# ---------------------------------------------

#!/bin/bash

# Loop through all .fastq.gz files in the current directory
for file in *.fastq.gz; do
    # Define the output file name
    output="sampled_${file}"

    # Run seqtk to sample 1% of the reads
    echo "Processing $file -> $output"
    zcat "$file" | seqtk sample -s42 - 0.01 | gzip > "$output"

    # Confirm completion
    echo "Finished sampling $file"
done

echo "All files processed!"