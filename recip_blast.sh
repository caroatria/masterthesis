#!/bin/bash
database_file="/Users/carolinaatria/Desktop/ADesktop/Studium/Master/master_thesis/sponge/data_from_roger/python/data.nosync/Sdomuncula_IsoseqStringtie/transcriptome/longest_isoform/transcriptome_updated.fa"

#QUERY
query_file=$1
outfile=$2
makeblastdb -in $database_file -dbtype nucl

tblastn -query $query_file -db $database_file -max_target_seqs 1 -outfmt "10 qseqid sseqid evalue" > $outfile
