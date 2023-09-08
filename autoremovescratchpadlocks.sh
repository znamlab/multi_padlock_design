#!/bin/bash

# Navigate to the root directory
cd /nemo/lab/znamenskiyp/scratch/

# Find 'v1_genes_18_06_2023_chunk__*' directories
for dir in v1_genes_18_06_2023_chunk_*/; do
  if [ -d "$dir" ]; then
    # Find and delete 'TempFolder*' directories within 'olf_chunk_*' directories
    find "$dir" -type d -name 'TempFolder*' -exec rm -r '{}' \;
  fi
done
