library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript general_graphing.R <bedmethyl_file> <color_metadata_file> <chromosome2species_file> <output_dir> <modification> [yeast_chromosome] [output_prefix] [cov_cutoff] [mod_cutoff]")
}

bedmethyl_file <- args[1]
color_metadata_file <- args[2]
chromosome2species_file <- args[3]
output_dir <- args[4]
modification <- args[5]
yeast_chromosome <- if (length(args) >= 6) args[6] else NA
output_prefix <- if (length(args) >= 7) args[7] else tools::file_path_sans_ext(basename(bedmethyl_file))
cov_cutoff <- if (length(args) >= 8) as.numeric(args[8]) else 10
mod_cutoff <- if (length(args) >= 9) as.numeric(args[9]) else 7

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

cat("Running bedmethyl analysis for modification:", modification, "\n")

####----------------------- colors metadata
mycolors_df <- read.table(color_metadata_file, sep = "\t", h = T, comment.char = "")
mycolors_scale <- mycolors_df$ColorsName
names(mycolors_scale) <- mycolors_df$Source_original

####----------------------- cutoffs and labeling helper
get_cutoff_counts <- function(df, cov_cutoff, mod_cutoff) {
  # hard cutoff on coverage: must be strictly greater than cov_cutoff to pass
  n_cov_below <- sum(df$Nvalid_cov <= cov_cutoff, na.rm = TRUE)
  n_cov_above <- sum(df$Nvalid_cov > cov_cutoff, na.rm = TRUE)
  n_mod_below <- sum(df$percent_modified < mod_cutoff, na.rm = TRUE)
  n_mod_above <- sum(df$percent_modified >= mod_cutoff, na.rm = TRUE)
  list(
    n_cov_below = n_cov_below,
    n_cov_above = n_cov_above,
    n_mod_below = n_mod_below,
    n_mod_above = n_mod_above
  )
}

get_cutoff_label <- function(df, cov_cutoff, mod_cutoff) {
  counts <- get_cutoff_counts(df, cov_cutoff, mod_cutoff)
  paste0(
    "Nvalid_cov <= ", cov_cutoff, ": n=", format(counts$n_cov_below, big.mark = ",", scientific = FALSE),
    "   |   Nvalid_cov > ", cov_cutoff, ": n=", format(counts$n_cov_above, big.mark = ",", scientific = FALSE),
    "\n% modified < ", mod_cutoff, ": n=", format(counts$n_mod_below, big.mark = ",", scientific = FALSE),
    "   |   % modified >= ", mod_cutoff, ": n=", format(counts$n_mod_above, big.mark = ",", scientific = FALSE)
  )
}

####----------------------- read bedmethyl
bedmethyl_data <- read.table(bedmethyl_file, sep = "\t", h = F)
dim(bedmethyl_data)

####----------------------- reference and metadata
bedmethyl_headers <- c("chromosome", "start", "end", "modification", "score", "strand", "display_start", "display_end", "display_color", "Nvalid_cov", "percent_modified")
bedmethyl_modkit_pileup_headers <- c(bedmethyl_headers, "Nmod", "Ncanonical", "Nother_mod", "Ndeleted", "Nfail", "Ndiff_mod", "Nno_call")
chromosme2species <- read.table(chromosome2species_file, h = T)

####--------------------------------------------------
colnames(bedmethyl_data) <- bedmethyl_modkit_pileup_headers
bedmethyl_data <- merge(bedmethyl_data, chromosme2species, by = "chromosome", all.x = T)
dim(bedmethyl_data)

####----------------------- cutoff counts: print to console/log
all_counts <- get_cutoff_counts(bedmethyl_data, cov_cutoff, mod_cutoff)
cat("---- Cutoff summary (all data, modification:", modification, ") ----\n")
cat("Nvalid_cov <=", cov_cutoff, ": n =", format(all_counts$n_cov_below, big.mark = ",", scientific = FALSE), "\n")
cat("Nvalid_cov >", cov_cutoff, ": n =", format(all_counts$n_cov_above, big.mark = ",", scientific = FALSE), "\n")
cat("percent_modified <", mod_cutoff, ": n =", format(all_counts$n_mod_below, big.mark = ",", scientific = FALSE), "\n")
cat("percent_modified >=", mod_cutoff, ": n =", format(all_counts$n_mod_above, big.mark = ",", scientific = FALSE), "\n")

####----------------------- scatter plot 1
myplot <- ggplot(data = bedmethyl_data, aes(x = Nvalid_cov, y = percent_modified)) +
	theme_bw() +
	geom_point(aes(colour = Source), alpha = 0.4) +
	geom_vline(xintercept = cov_cutoff, linetype = "dashed", colour = "black") +
	geom_hline(yintercept = mod_cutoff, linetype = "dashed", colour = "black") +
	scale_colour_manual(values = mycolors_scale) +
	ggtitle(paste0(output_prefix, ": ", modification)) +
	labs(subtitle = get_cutoff_label(bedmethyl_data, cov_cutoff, mod_cutoff))
myplot <- myplot + scale_x_log10()
png(paste0(output_prefix, ".", modification, ".plot1.png"), width = 8, height = 6, units = "in", res = 300)
print(myplot)
dev.off()

####----------------------- single-chromosome subset (e.g. yeast spike-in)
if (!is.na(yeast_chromosome)) {
	iYeast <- bedmethyl_data$chromosome == yeast_chromosome
	pseU_yeast <- bedmethyl_data[iYeast, ]
	dim(pseU_yeast)
	yeast_counts <- get_cutoff_counts(pseU_yeast, cov_cutoff, mod_cutoff)
	cat("---- Cutoff summary (", yeast_chromosome, " only, modification:", modification, ") ----\n")
	cat("Nvalid_cov <=", cov_cutoff, ": n =", format(yeast_counts$n_cov_below, big.mark = ",", scientific = FALSE), "\n")
	cat("Nvalid_cov >", cov_cutoff, ": n =", format(yeast_counts$n_cov_above, big.mark = ",", scientific = FALSE), "\n")
	cat("percent_modified <", mod_cutoff, ": n =", format(yeast_counts$n_mod_below, big.mark = ",", scientific = FALSE), "\n")
	cat("percent_modified >=", mod_cutoff, ": n =", format(yeast_counts$n_mod_above, big.mark = ",", scientific = FALSE), "\n")
	plotYeast <- ggplot(data = pseU_yeast, aes(x = Nvalid_cov, y = percent_modified)) +
		theme_bw() +
		geom_point(aes(colour = Source), alpha = 0.4) +
		geom_vline(xintercept = cov_cutoff, linetype = "dashed", colour = "black") +
		geom_hline(yintercept = mod_cutoff, linetype = "dashed", colour = "black") +
		scale_x_log10() +
		scale_colour_manual(values = mycolors_scale) +
		ggtitle(paste0(output_prefix, ": ", modification, ": ", yeast_chromosome, " only")) +
		labs(subtitle = get_cutoff_label(pseU_yeast, cov_cutoff, mod_cutoff))
	png(paste0(output_prefix, ".", modification, ".plot2.", yeast_chromosome, "_only.png"), width = 8, height = 6, units = "in", res = 300)
	print(plotYeast)
	dev.off()
}

####----------------------- boxplots
boxplot1 <- ggplot(data = bedmethyl_data, aes(x = Source, y = Nvalid_cov)) +
	geom_boxplot(aes(fill = Source)) +
	theme_bw() +
	scale_y_log10() +
	geom_hline(yintercept = cov_cutoff, linetype = "dashed", colour = "black") +
	scale_fill_manual(values = mycolors_scale) +
	ggtitle(paste0(output_prefix, ": ", modification, ": Coverage by Source"))
png(paste0(output_prefix, ".", modification, ".boxplot1.coverage.png"), width = 8, height = 6, units = "in", res = 300)
print(boxplot1)
dev.off()

boxplot2 <- ggplot(data = bedmethyl_data, aes(x = Source, y = percent_modified)) +
	geom_boxplot(aes(fill = Source)) +
	theme_bw() +
	scale_y_log10() +
	geom_hline(yintercept = mod_cutoff, linetype = "dashed", colour = "black") +
	scale_fill_manual(values = mycolors_scale) +
	ggtitle(paste0(output_prefix, ": ", modification, ": Percent Modified by Source"))
png(paste0(output_prefix, ".", modification, ".boxplot2.pct_modified.png"), width = 8, height = 6, units = "in", res = 300)
print(boxplot2)
dev.off()

####----------------------- summaries
summary(bedmethyl_data$percent_modified)
aggregate(x = bedmethyl_data, percent_modified ~ Source, FUN = summary)
aggregate(x = bedmethyl_data, Nvalid_cov ~ Source, FUN = summary)
cat("DONE:", format(Sys.time()), "\n")
