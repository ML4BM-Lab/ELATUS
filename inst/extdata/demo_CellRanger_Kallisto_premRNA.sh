# Guidelines to use Cell Ranger count, Kallisto and Kallisto multi preprocessing pipelines for generating the raw count matrices in snRNA-seq and generate the inputs needed for running ELATUS_premRNA
# hg38 (Download annotation and fasta files from GENCODE)
hg38_gencode_genome_fasta
hg38_gencode_gtf
hg38_gencode_transcriptome_fasta

# 1. Cell Ranger (v.5)
# 1.1. Index generation: 
cellranger mkref --genome=$CellRanger_reference --fasta=$hg38_gencode_genome_fasta  --genes=$hg38_gencode_gtf
# 1.2. Generate raw count matrix including introns 
cellranger count --id=$id --transcriptome=$CellRanger_reference  --fastqs=$FASTQ_DIR --sample=$sample --include-introns

# 2. Kallisto (v.0.46.1) 
# 2.1. Index generation (including introns. Generate the index according to: https://www.biostars.org/p/468180/) 
ml Python/3.9.6-GCCcore-11.2.0
pip install kb-python
kb ref -i kallisto_nuclei_introns_index.idx -g kallisto_nuclei_introns_t2g.txt -f1 kallisto_nuclei_introns_cdna.fa -f2 kallisto_nuclei_introns_intron.fa -c1 kallisto_nuclei_introns_cdna_t2c.txt -c2 kallisto_nuclei_introns_t2c.txt --workflow lamanno -n 8 hg38_gencode_genome_fasta hg38_gencode_gtf
cat kallisto_nuclei_introns_cdna.fa kallisto_nuclei_introns_intron.fa > kallisto_nuclei_cDNA_introns_ALL.fa
kallisto index -i /home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcriptome_index_kallisto.idx kallisto_nuclei_cDNA_introns_ALL.fa

# 2.2 Generate raw count matrix
KALLISTO_INDEX="/home/egonie/dato-activo/reference.genomes_kike/GRCm39/gencode/transcriptome_index_kallisto.idx"
$TENxV3_WHITELIST (Download 10x whitelist from 10x website (https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist))
TRANSCRIPTS_TO_GENES="extdata/tr2g_mm10.tsv"
kallisto bus -i $KALLISTO_INDEX -o output_intronic -x 10xv3 -t 8 $fasta_R1 $fasta_R2
cd $output_intronic
cat transcripts.txt | sed "s/|.*//" > transcripts_back.txt
rm transcripts.txt
mv transcripts_back.txt transcripts.txt
mkdir bustools_results
cd bustools_results
# 2.2.1. Kallisto without multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_NO_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts -

# 2.2.2. Kallisto with multimappers
bustools correct -w $TENxV3_WHITELIST -p ../output.bus | bustools sort -T temp/ -p - | bustools count -o cells_genes_multimapping -g $TRANSCRIPTS_TO_GENES -e ../matrix.ec -t ../transcripts.txt --genecounts --multimapping -

