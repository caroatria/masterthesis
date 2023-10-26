#!/bin/bash
database_file="/Users/carolinaatria/Desktop/ADesktop/Studium/Master/master_thesis/sponge/data_from_roger/python/data.nosync/mus_musculus/proteome/uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2023.06.03-16.55.21.10.fasta"

#QUERY
query_file=$1
outfile=$2

makeblastdb -in $database_file -dbtype prot

blastx -query $query_file -db $database_file -max_target_seqs 1 -outfmt "10 qseqid sseqid evalue" -evalue 0.001 > $outfile
