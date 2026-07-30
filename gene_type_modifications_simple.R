args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript gene_type_modifications_simple.R <bedmethyl_file> <gtf_path> <output_dir> <modification> [cov_cutoff] [mod_cutoff] [output_prefix]")
}

bedmethyl_file <- args[1]
gtf_path <- args[2]
output_dir <- args[3]
modification <- args[4]
cov_cutoff <- if (length(args) >= 5) as.numeric(args[5]) else 10
mod_cutoff <- if (length(args) >= 6) as.numeric(args[6]) else 7
output_prefix <- if (length(args) >= 7) args[7] else tools::file_path_sans_ext(basename(bedmethyl_file))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

## 1. read bedmethyl + apply cutoffs to the whole dataset ---------------------
hdr <- c("chromosome", "start", "end", "modification", "score", "strand",
         "display_start", "display_end", "display_color", "Nvalid_cov", "percent_modified",
         "Nmod", "Ncanonical", "Nother_mod", "Ndeleted", "Nfail", "Ndiff_mod", "Nno_call")
bm <- read.table(bedmethyl_file, sep = "\t", h = FALSE)
colnames(bm) <- hdr
passed <- bm[bm$Nvalid_cov > cov_cutoff & bm$percent_modified >= mod_cutoff, ]
cat("Sites passing cutoffs (Nvalid_cov >", cov_cutoff, "& percent_modified >=", mod_cutoff, "):",
    nrow(passed), "\n")
write.table(passed[, c("chromosome", "start", "end")], "passed_sites.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
## 2. overlap passing sites with gene features, grab gene_biotype -------------
system(paste0("awk -F'\t' '$3==\"gene\"' ", gtf_path, " > genes_only.gtf"))
system("bedtools intersect -a passed_sites.bed -b genes_only.gtf -wb > hits.tsv")
hits <- read.table("hits.tsv", sep = "\t", quote = "", stringsAsFactors = FALSE)
attrs   <- hits[[ncol(hits)]]                                  # the GTF attribute blob
biotype <- sub('.*gene_biotype "([^"]+)".*', "\\1", attrs)     # pull the value
biotype[!grepl('gene_biotype "', attrs)] <- "unknown"          # genes with no biotype tag
## (optional) flag ribosomal protein genes by their gene name (rps/rpl/rpm) ---
gene_name <- sub('.*gene "([^"]+)".*', "\\1", attrs)
biotype[biotype == "protein_coding" & grepl("^rp[slm]", gene_name, ignore.case = TRUE)] <- "ribosomal_protein"
## 3. count + rank -----------------------------------------------------------
tab <- as.data.frame(table(biotype), stringsAsFactors = FALSE)
colnames(tab) <- c("gene_type", "n_modifications")
tab <- tab[order(-tab$n_modifications), ]
tab$pct <- round(100 * tab$n_modifications / sum(tab$n_modifications), 1)
print(tab, row.names = FALSE)
write.table(tab, paste0(output_prefix, ".", modification, ".gene_type_modifications.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("DONE:", format(Sys.time()), "\n")
