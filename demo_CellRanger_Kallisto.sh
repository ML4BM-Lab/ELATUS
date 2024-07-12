# Guidelines to use Cell Ranger count and Kallisto preprocessing pipelines for generating the raw count matrices
# mm10 (Download annotation and fasta files from GENCODE)
mm39_gencode_genome_fasta
mm39_gencode_gtf
mm39_gencode_transcriptome_fasta

# 1. Cell Ranger (v.3.0.1)
# 1.1. Index generation: 
cellranger mkref --genome=$CellRanger_reference --fasta=$mm39_gencode_genome_fasta  --genes=$mm39_gencode_gtf
# 1.2. Generate raw count matrix
cellranger count --id=$id --transcriptome=$CellRanger_reference  --fastqs=$FASTQ_DIR --sample=$sample 

# 2. Kallisto (v.0.46.1)
# 2.1. Index generation: 
kallisto index -i $transcriptome_index_kallisto $mm10_gencode_transcriptome_fasta
# 2.2 Generate raw count matrix
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcriptome_index_kallisto.idx"
$TENxV3_WHITELIST (Download 10x whitelist from 10x website (https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist))
TRANSCRIPTS_TO_GENES="extdata/tr2g_mm10.tsv"
kallisto bus -i $transcriptome_index_kallisto -o $output -x $chemistry -t 8 $fasta_R1 $fasta_R2
cd $output
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
mkdir bustools_results
cd bustools_results
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -




