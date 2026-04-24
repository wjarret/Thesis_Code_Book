library(dada2)
library(openxlsx)
library(dplyr)
library(Biostrings)
library(optparse)


option_list <- list(
  make_option(c("--r1"), type="character", help="Forward FASTQ file"),
  make_option(c("--r2"), type="character", help="Reverse FASTQ file"),
  make_option(c("--out"), type="character", help="Output FASTA file")
)

opt <- parse_args(OptionParser(option_list=option_list))

r1_file <- opt$r1
r2_file <- opt$r2
output_file <- opt$out

filt_dir <- file.path(tempdir(), "filtered")
dir.create(filt_dir, showWarnings = FALSE)

# For a single sample, just define filtered file names directly
filtF <- file.path(filt_dir, "F_filt.fastq.gz")
filtR <- file.path(filt_dir, "R_filt.fastq.gz")
print(filtF)
print(filtR)

#Filtered dir should still be empty right before filtering. This is because dada2 is operating in append mode. This has led to a lot of issues
#and hopefully it runs through smoothly
# 
# cat("\nContents of filtered immediately before filterAndTrim:\n")
# print(list.files(filt_dir, full.names = TRUE))

# Filter and trim (removes reads with N because maxN=0)
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, compress = TRUE, multithread = TRUE)
out <- filterAndTrim(r1_file, filtF, r2_file, filtR, maxN = 0, compress = TRUE, multithread = TRUE)

print(out)

#learning error rate
errF <- learnErrors(filtF, multithread = TRUE)
errR <- learnErrors(filtR, multithread = TRUE)


# get folder where FASTA will be written
out_dir <- dirname(output_file) 

#write plot to PDF instead of screen
pdf(file = file.path(out_dir, "snRNA_error_plot_F.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(file = file.path(out_dir, "snRNA_error_plot_R.pdf"))
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Denoise sequences
dadaFs <- dada(filtF, err = errF, multithread = TRUE)
dadaRs <- dada(filtR, err = errR, multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, filtF, dadaRs, filtR, verbose = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
print(dim(seqtab))

# Inspect distribution of sequence lengths
len_tab <- table(nchar(getSequences(seqtab)))
print(len_tab)

asv_seqs    <- colnames(seqtab)
asv_lengths <- nchar(asv_seqs)


# Write to Excel (saved in the working dir where the job runs)
out_xlsx <- file.path(out_dir, "sequence_table_with_sequences_U1.xlsx")
write.xlsx(mergers, out_xlsx, rowNames = FALSE)
cat("Wrote Excel file to:", out_xlsx, "\n")

#Construct Fasta file
sequences <- mergers$sequence

#Make sequences column identifiable as DNA sequences
#Biostrings allows us to do that

DNA_sequences<- DNAStringSet(sequences)

#Names each ASV
names(DNA_sequences) <- paste0('ASV', 1:length(DNA_sequences))
DNA_sequences <- reverseComplement(DNA_sequences)
print(DNA_sequences)

#Outputs fasta to directory
writeXStringSet(x = DNA_sequences, filepath = output_file)
