ncbi-blast-2.5.0+/bin/blastn -query ref_cdna.fasta -db ref_cdna.fasta.masked \
        -max_target_seqs 10000 -outfmt 6 \
        -evalue 1e-3 -lcase_masking \
        -num_threads  6 \
        -word_size 11  >  blast_pairs.outfmt6
