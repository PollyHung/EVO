module load bwa
module load Picard
module load samtools

HG38_REF="$REF/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
KNOWN_INDELS="$REF/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILL_INDELS="$REF/vcf/Mills.indels.contig.adjusted.hg38.vcf.gz"
DBSNP="$REF/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz"

bwa mem -M -t 12 \
  "$HG38_REF" \
  "$SEQ_1" "$SEQ_2" \
  2> "$error" \
  > "$SAM"

java -jar /software/Picard/3.2.0/picard.jar SortSam \
  --INPUT "$SAM" \
  --OUTPUT "$SORTED_SAM" \
  --SORT_ORDER coordinate \
  --VALIDATION_STRINGENCY SILENT

java -jar /software/Picard/3.2.0/picard.jar AddOrReplaceReadGroups \
  --INPUT "$SORTED_SAM" \
  --OUTPUT "$TAGGED_SAM" \
  --RGLB "Short-Insert_Library" \
  --RGPL "DNBSEQ" \
  --RGPU "V350264599.2" \
  --RGSM "$sample_id"

java -jar /software/Picard/3.2.0/picard.jar MarkDuplicates \
  --INPUT "$TAGGED_SAM" \
  --OUTPUT "$MARKED_BAM" \
  --METRICS_FILE "$METRICS" \
  --ASSUME_SORTED true \
  --REMOVE_DUPLICATES true \
  --VALIDATION_STRINGENCY SILENT

cd "$RESULT"
samtools index "$MARKED_BAM"



