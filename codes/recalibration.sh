module load miniconda3
module load GenomeAnalysisTK/4.2.0.0

gatk BaseRecalibrator \
	-I "$MARKED_BAM" \
	-R "$HG38_REF" \
	--known-sites "$KNOWN_INDELS" \
	--known-sites "$MILL_INDELS" \
	--known-sites "$DBSNP" \
	-O "$RECAL_TABLE"

gatk ApplyBQSR \
	-R "$HG38_REF" \
	-I "$MARKED_BAM" \
	--bqsr-recal-file "$RECAL_TABLE" \
	-O "$CALIBRATED_BAM"
