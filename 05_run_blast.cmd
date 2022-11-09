echo "Running blast against linkers" 
blastn -query ../reads_combined.fasta -db ../linkers -out ../linker_alignments.txt -task blastn-short -max_target_seqs 1 -outfmt 6 -num_threads 6

echo "Running blast against targets" 
blastn -query ../reads_combined.fasta -db ../genomic_targets -out ../genomic_target_alignments.txt -task blastn-short -max_target_seqs 1 -outfmt 6 -num_threads 6