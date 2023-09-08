awk 'BEGIN {n_seq=0; genes_in_chunk=0; max_genes_per_chunk=1; current_file=sprintf("<name_here>_chunk_%d.fasta", n_seq);} /^>/ { 
        split($0, header, ","); 
        gene_name = header[1]; 
        if (!(gene_name in gene_ids)) {
            gene_ids[gene_name] = 1;
            genes_in_chunk++;
            if (genes_in_chunk > max_genes_per_chunk) {
                close(current_file);
                n_seq++;
                current_file = sprintf("<name_here>_chunk_%d.fasta", n_seq);
                genes_in_chunk = 1;
            }
        }
        print >> current_file; 
        next; 
    } { 
        print >> current_file; 
    }' < <name_here>.fasta
