library(QDNAseq)
library(dplyr)
library(tidyr)
library(magrittr)
library(QDNAseq.hg38)
library(tibble)

## Expected Folder Structure: 
## /Folder
##      /data/BAM
##      /data/FASTQ
##      /data/QDNAseq
##      /data/sample_id.txt
##      /data/reference_cn_unnormalised.txt
##      /scripts

# Set Overall Working Directory 
HOME="/Volumes/Fernyiges/hgsoc_cell-line_qdnaseq/YiHan/"
setwd(HOME)

# Get Sample IDS
samples <- read.delim("data/sample_id.txt", header = F) %>% unlist

# Get CN reference from FT190 
ref <- read.delim("data/reference_cn_unnormalised.txt", row.names = 1)

# Get hg38 ref
bins <- getBinAnnotations(binSize=1000, genome="hg38")

# Loop For Each Sample 
for(SAMPLE in samples){
  
  # set working directory to analysis 
  ANALYSIS=file.path(HOME, "data/QDNAseq", SAMPLE)
  if(!dir.exists(ANALYSIS)){
    dir.create(ANALYSIS)
    setwd(ANALYSIS)
  } else {
    setwd(ANALYSIS)
  }
  
  # Load bam files
  bam <- paste0(HOME, "data/BAM/", SAMPLE, ".sort.tag.dedup.cal.bam")
  readCounts <- binReadCounts(bins, bamfiles = bam)

  # Visualise raw copy numbers
  png("01_raw_copy_numbers.png", width = 8, height = 4, units = "in", res = 600)
  plot(readCounts, logTransform=FALSE)
  highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
  dev.off()

  # Filtering
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
  png("02_median_read_counts_per_bin.png", width = 5, height = 4, units = "in", res = 600)
  isobarPlot(readCountsFiltered)
  dev.off()

  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  png("03_sequence_depth_and_noise.png", width = 5, height = 4, units = "in", res = 600)
  noisePlot(readCountsFiltered)
  dev.off()

  # Correction for GC content
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  png("04_CN_profile_after_GC_correction.png", width = 8, height = 4, units = "in", res = 600)
  plot(copyNumbersSmooth)
  dev.off()

  # Export 1
  exportBins(copyNumbersSmooth, file = paste0(SAMPLE, ".txt"), logTransform = FALSE)
  exportBins(copyNumbersSmooth, file = paste0(SAMPLE, ".igv"), format = "igv")
  exportBins(copyNumbersSmooth, file = paste0(SAMPLE, ".bed"), format = "bed")

  # Segmentation
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "log2")
  png("05_CN_profile_after_segmentation.png", width = 8, height = 4, units = "in", res = 600)
  plot(copyNumbersSegmented)
  dev.off()

  # Call copy number
  copyNumbersCalled <- callBins(copyNumbersSegmented)
  png("06_CN_profile_after_calling_gains_and_losses.png", width = 8, height = 4, units = "in", res = 600)
  plot(copyNumbersCalled)
  dev.off()

  # Deduct CNV
  new_cn <- merge(copyNumbersSmooth@assayData$copynumber, ref, by = "row.names", all = TRUE) %>%
    column_to_rownames(var = "Row.names")
  new_cn[, 1] <- new_cn[, 1] / new_cn[, 2]
  new_cn <- new_cn[rownames(copyNumbersSmooth@assayData$copynumber), ]
  new_cn$reference <- NULL

  # Update CNV
  copyNumbersFitted <- copyNumbersSmooth
  unlockBinding("copynumber", copyNumbersFitted@assayData)
  copyNumbersFitted@assayData$copynumber <- new_cn
  lockBinding("copynumber", copyNumbersFitted@assayData)
  
  # Re-bin
  copyNumbersFittedSegmented <- segmentBins(copyNumbersFitted, transformFun = "log2")
  png("07_CN_calibrated_segmented_Bins.png", width = 8, height = 4, units = "in", res = 600)
  plot(copyNumbersFittedSegmented)
  dev.off()
  
  # Re-call 
  copyNumbersFittedBinned <- callBins(copyNumbersFittedSegmented, organism = "human", method = "cutoff")
  png("08_CN_calibrated_call_Bins.png", width = 8, height = 4, units = "in", res = 600)
  plot(copyNumbersFittedBinned)
  dev.off()
  
  # Export 2
  exportBins(copyNumbersCalled, file = paste0(SAMPLE, ".vcf"), format = "vcf")
  exportBins(copyNumbersCalled, file = paste0(SAMPLE, ".seg"), format = "seg")
  
  # Save 
  save.image("image.RData")
}  



