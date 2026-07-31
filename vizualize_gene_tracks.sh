#!/bin/bash
# visualize_gene_tracks.sh
# Picks a random subset of modified genes, then plots a Lollipop track view of selected genes

set -uo pipefail

# ---- input files
CSV="/path/to/analyze_genes_per_gene.csv"   # per-gene CSV from analyze_genes.sh
SITES="/path/to/tag_sites.bed"              # filtered sites BED from analyze_genes.sh (<tag>_sites.bed)
GBED="/path/to/tag_genes.bed"               # gene BED from analyze_genes.sh (<tag>_genes.bed)
GTF="/path/to/annotation.gtf"               # GTF for product lookup and exon extraction
OUTDIR="/path/to/outdir"                    # output directory
TAG="mytag"                                 # label for outputs
NGENES=10                                   # number of genes to randomly select
MAXSITES=100                                # max n_sites to consider

# Requires: bedtools, awk, Rscript (ggplot2, trackViewer + GenomicRanges optional)

# validate that all required files exist and tools are on PATH
for f in "$CSV" "$SITES" "$GBED" "$GTF"; do
  [ -s "$f" ] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done
for cmd in bedtools awk Rscript; do
  command -v "$cmd" >/dev/null || { echo "ERROR: '$cmd' not on PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

# -------- Step 1: randomly select NGENES genes below the MAXSITES cutoff
SELECTED="${OUTDIR}/${TAG}_selected_genes.tsv"
awk -F',' -v mx="$MAXSITES" '
  NR==1 { next }
  $3 > mx { next }
  { printf "%s\t%s\t%d\t%.3f\t%.1f\n", $1, $2, $3, $4, $3*$4 }
' "$CSV" | shuf -n "$NGENES" | sort -t$'\t' -k3,3n > "$SELECTED"

NSEL=$(wc -l < "$SELECTED")
[ "$NSEL" -gt 0 ] || { echo "ERROR: no genes selected" >&2; exit 1; }

# -------- Step 2: extract human-readable product names from GTF for each selected gene
PRODUCTS="${OUTDIR}/${TAG}_selected_products.tsv"
while IFS=$'\t' read -r gid gname nsites mean score; do
  # grep GTF for first line with this gene_id, pull the product="..." attribute
  prod=$(grep -m1 "\"$gid\"" "$GTF" | grep -oE 'product "[^"]+"' | head -1 | sed 's/^product "//; s/"$//')
  [ -z "$prod" ] && prod="$gname"
  # clean up noisy suffixes for shorter plot labels
  prod=$(echo "$prod" | sed 's/, transcript variant.*//; s/ homolog//; s/ isoform.*//')
  printf "%s\t%s\t%d\t%.3f\t%s\n" "$gid" "$gname" "$nsites" "$mean" "$prod"
done < "$SELECTED" > "$PRODUCTS"

# -------- Step 3: normalize mod-site and exon positions to 0-100% of gene body (0%=5', 100%=3')
TRACKS="${OUTDIR}/${TAG}_track_data.tsv"
: > "$TRACKS"
EXONS="${OUTDIR}/${TAG}_exon_structure.tsv"
: > "$EXONS"
# raw-coordinate versions for trackViewer (needs genomic positions, not 0-100%)
EXONS_RAW="${OUTDIR}/${TAG}_exon_raw.tsv"
: > "$EXONS_RAW"
GENE_COORDS="${OUTDIR}/${TAG}_gene_coords.tsv"
: > "$GENE_COORDS"

while IFS=$'\t' read -r gid gname nsites mean prod; do
  # look up this gene's genomic coordinates from the gene BED
  gene_line=$(awk -F'\t' -v g="$gid" 'index($4, g)==1' "$GBED" | head -1)
  [ -z "$gene_line" ] && continue
  gchr=$(echo "$gene_line" | cut -f1)
  gstart=$(echo "$gene_line" | cut -f2)
  gend=$(echo "$gene_line" | cut -f3)
  gstrand=$(echo "$gene_line" | cut -f6)
  glen=$((gend - gstart))
  [ "$glen" -le 0 ] && continue

  # save raw gene coordinates for trackViewer
  printf "%s\t%s\t%s\t%d\t%d\t%s\t%s\n" "$gid" "$prod" "$gchr" "$gstart" "$gend" "$gstrand" "$nsites" >> "$GENE_COORDS"

  # 3a: find all exon features for this gene in GTF, normalize start/end to 0-100%
  awk -F'\t' -v gid="$gid" '
    $3 == "exon" && index($9, "\"" gid "\"") > 0 {
      print $4, $5
    }
  ' "$GTF" | sort -k1,1n | \
  while read -r estart eend; do
    if [ "$gstrand" = "-" ]; then
      xstart_pct=$(awk -v ge="$gend" -v e="$eend"  -v gl="$glen" 'BEGIN{printf "%.2f", (ge - e)   * 100 / gl}')
      xend_pct=$(awk   -v ge="$gend" -v e="$estart" -v gl="$glen" 'BEGIN{printf "%.2f", (ge - e) * 100 / gl}')
    else
      xstart_pct=$(awk -v gs="$gstart" -v e="$estart" -v gl="$glen" 'BEGIN{printf "%.2f", (e - gs) * 100 / gl}')
      xend_pct=$(awk   -v gs="$gstart" -v e="$eend"   -v gl="$glen" 'BEGIN{printf "%.2f", (e - gs) * 100 / gl}')
    fi
    printf "%s\t%s\t%s\t%s\n" "$gid" "$prod" "$xstart_pct" "$xend_pct"
    # also save raw genomic coordinates for trackViewer
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$gid" "$prod" "$gchr" "$estart" "$eend" "$gstrand" >> "$EXONS_RAW"
  done >> "$EXONS"

  # 3b: find all mod sites overlapping this gene, normalize each position to 0-100%
  awk -F'\t' -v c="$gchr" -v s="$gstart" -v e="$gend" -v st="$gstrand" \
      '$1==c && $2>=s && $3<=e && $6==st {print $0}' "$SITES" | \
  while IFS=$'\t' read -r chr start end pct cov strand; do
    if [ "$gstrand" = "-" ]; then
      pos_norm=$(awk -v ge="$gend" -v s="$start" -v gl="$glen" 'BEGIN{printf "%.2f", (ge-s-1)*100/gl}')
    else
      pos_norm=$(awk -v gs="$gstart" -v s="$start" -v gl="$glen" 'BEGIN{printf "%.2f", (s-gs)*100/gl}')
    fi
    printf "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\n" "$gid" "$prod" "$nsites" "$glen" "$pos_norm" "$pct" "$cov" "$nsites"
  done >> "$TRACKS"
done < "$PRODUCTS"

NTRACK=$(wc -l < "$TRACKS")
[ "$NTRACK" -gt 0 ] || { echo "ERROR: no track data produced" >&2; exit 1; }

# add column headers to the data files for R
TRACKS_H="${OUTDIR}/${TAG}_track_data_h.tsv"
{ echo "gene_id	product	n_sites	gene_length	pos_pct	pct_mod	coverage	nsites_bin"; cat "$TRACKS"; } > "$TRACKS_H"

EXONS_H="${OUTDIR}/${TAG}_exon_structure_h.tsv"
{ echo "gene_id	product	xstart_pct	xend_pct"; cat "$EXONS"; } > "$EXONS_H"

# add headers to raw-coordinate files for trackViewer
EXONS_RAW_H="${OUTDIR}/${TAG}_exon_raw_h.tsv"
{ printf "gene_id\tproduct\tchr\texon_start\texon_end\tstrand\n"; cat "$EXONS_RAW"; } > "$EXONS_RAW_H"

GENE_COORDS_H="${OUTDIR}/${TAG}_gene_coords_h.tsv"
{ printf "gene_id\tproduct\tchr\tgene_start\tgene_end\tstrand\tn_sites\n"; cat "$GENE_COORDS"; } > "$GENE_COORDS_H"

# check whether GTF had exon features; if not, R will fall back to flat gene bars
NEXONS=$(wc -l < "$EXONS")
if [ "$NEXONS" -eq 0 ]; then
  echo "WARNING: no exon features found in GTF for selected genes." >&2
  echo "         The plot will fall back to uniform gene body bars." >&2
  echo "         Check that your GTF contains 'exon' feature lines." >&2
  HAS_EXONS="FALSE"
else
  HAS_EXONS="TRUE"
fi

# -------- Step 4: compute per-gene summary stats (mean/median/min/max of pct_mod and coverage)
STATS="${OUTDIR}/${TAG}_selected_gene_stats.csv"
{
  echo "gene_id,product,n_sites,gene_length_bp,mean_pct_mod,median_pct_mod,max_pct_mod,min_pct_mod,mean_cov,median_cov"
  awk -F'\t' '
    # accumulate per-gene sums, min/max, and value lists for median
    {
      g=$1; n[g]++; prod[g]=$2; glen[g]=$4
      sum_p[g]+=$6; sum_c[g]+=$7
      if ($6>max_p[g]) max_p[g]=$6
      if (!(g in min_p) || $6<min_p[g]) min_p[g]=$6
      plist[g]=(plist[g]==""? $6 : plist[g]" "$6)
      clist[g]=(clist[g]==""? $7 : clist[g]" "$7)
    }
    END {
      for (g in n) {
        # bubble-sort values then pick middle element(s) for median
        m=split(plist[g],v," ")
        for(i=1;i<=m;i++) for(j=i+1;j<=m;j++) if(v[i]>v[j]){t=v[i];v[i]=v[j];v[j]=t}
        if(m%2) medp=v[(m+1)/2]; else medp=(v[m/2]+v[m/2+1])/2
        m2=split(clist[g],w," ")
        for(i=1;i<=m2;i++) for(j=i+1;j<=m2;j++) if(w[i]>w[j]){t=w[i];w[i]=w[j];w[j]=t}
        if(m2%2) medc=w[(m2+1)/2]; else medc=(w[m2/2]+w[m2/2+1])/2
        printf "%s,%s,%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n", g, prod[g], n[g], glen[g], sum_p[g]/n[g], medp, max_p[g], min_p[g], sum_c[g]/n[g], medc
      }
    }
  ' "$TRACKS" | sort -t',' -k3,3nr
} > "$STATS"


# -------- Step 5: trackViewer lolliplot that shows intron/exon information, genomic coordinates of modifications, and degree of modification
if [ "$HAS_EXONS" = "TRUE" ]; then
  TV_PDF="${OUTDIR}/${TAG}_trackviewer_lollipops.pdf"
  Rscript - <<RSCRIPT3
suppressMessages({
  ok <- require(trackViewer, quietly = TRUE) &&
        require(GenomicRanges, quietly = TRUE)
})
if (!ok) {
  message("SKIP: trackViewer or GenomicRanges not installed.")
  message("      Install with: BiocManager::install('trackViewer')")
  quit(save = "no", status = 0)
}
suppressMessages(library(grid))

# ---- load data ----
genes  <- read.delim("${GENE_COORDS_H}", stringsAsFactors = FALSE)
exons  <- read.delim("${EXONS_RAW_H}",   stringsAsFactors = FALSE)
sites  <- read.delim("${SITES}", header = FALSE, sep = "\t",
                     col.names = c("chr","start","end","pct_mod","coverage","strand"),
                     stringsAsFactors = FALSE)

region_cols <- c("exon" = "#2c7bb6", "intron" = "#d7191c")
plotted <- 0

pdf("${TV_PDF}", width = 12, height = 5)

for (r in seq_len(nrow(genes))) {
  gid    <- genes\$gene_id[r]
  prod   <- genes\$product[r]
  gchr   <- genes\$chr[r]
  gstart <- genes\$gene_start[r]
  gend   <- genes\$gene_end[r]
  gstr   <- genes\$strand[r]
  ns     <- genes\$n_sites[r]

  # get exons and mod sites for this gene
  ex_g <- exons[exons\$gene_id == gid, , drop = FALSE]
  if (nrow(ex_g) == 0) next

  s_g <- sites[sites\$chr == gchr &
               sites\$start >= gstart &
               sites\$end   <= gend &
               sites\$strand == gstr, , drop = FALSE]
  if (nrow(s_g) == 0) next

  # ---- flip minus-strand coordinates so plot always reads 5'→3' left to right ----
  # For + strand: 5' = low genomic coord (left) — no change needed
  # For - strand: 5' = high genomic coord — mirror all positions around gene center
  if (gstr == "-") {
    s_g\$start_plot <- gstart + (gend - s_g\$start)
    s_g\$end_plot   <- gstart + (gend - s_g\$end)
    ex_g\$es_plot   <- gstart + (gend - ex_g\$exon_end)
    ex_g\$ee_plot   <- gstart + (gend - ex_g\$exon_start)
  } else {
    s_g\$start_plot <- s_g\$start
    s_g\$end_plot   <- s_g\$end
    ex_g\$es_plot   <- ex_g\$exon_start
    ex_g\$ee_plot   <- ex_g\$exon_end
  }

  # classify each site as exon or intron (using original coordinates)
  s_g\$region <- "intron"
  for (i in seq_len(nrow(s_g))) {
    if (any(s_g\$start[i] >= ex_g\$exon_start & s_g\$start[i] <= ex_g\$exon_end)) {
      s_g\$region[i] <- "exon"
    }
  }

  # ---- SNP GRanges: use flipped coordinates, set strand to + so trackViewer doesn't re-flip ----
  snp.gr <- GRanges(seqnames = gchr,
                    ranges   = IRanges(start = s_g\$start_plot, width = 1),
                    strand   = "+")
  snp.gr\$score  <- s_g\$pct_mod
  snp.gr\$color  <- ifelse(s_g\$region == "exon",
                           region_cols["exon"], region_cols["intron"])
  snp.gr\$border <- snp.gr\$color

  # ---- gene-model features: use flipped exon coordinates ----
  feat.gr <- GRanges(seqnames = gchr,
                     ranges   = IRanges(start = ex_g\$es_plot, end = ex_g\$ee_plot),
                     strand   = "+")
  feat.gr\$fill   <- rep("grey25", nrow(ex_g))
  feat.gr\$height <- rep(0.06, nrow(ex_g))

  feat.grl <- GRangesList(gene = feat.gr)
  names(feat.grl) <- prod

  # x-axis range = full gene span with 2% padding
  pad    <- max(1, round((gend - gstart) * 0.02))
  xrange <- GRanges(gchr, IRanges(gstart - pad, gend + pad))

  # ---- draw on a fresh page ----
  grid.newpage()

  # title with strand indicator
  strand_label <- ifelse(gstr == "+", "(+ strand)", "(\u2212 strand, flipped)")
  pushViewport(viewport(x = 0.5, y = 0.96, height = unit(1.2, "lines")))
  grid.text(paste0(prod, "  (", gid, ",  n=", ns, " sites)  ", strand_label),
            gp = gpar(fontsize = 13, fontface = "bold"))
  popViewport()

  # main plot area
  pushViewport(viewport(x = 0.5, y = 0.46, width = 0.92, height = 0.84))
  tryCatch(
    lolliplot(snp.gr, feat.grl, ranges = xrange, ylab = "% modified",
              newpage = FALSE),
    error = function(e) {
      tryCatch(
        lolliplot(snp.gr, feat.grl, ranges = xrange, ylab = "% modified"),
        error = function(e2) message("  plot failed for ", gid, ": ", e2\$message)
      )
    }
  )
  popViewport()

  # ---- 5' / 3' labels at bottom corners ----
  pushViewport(viewport(x = 0.08, y = 0.08, width = 0.05, height = 0.04))
  grid.text("5'", gp = gpar(fontsize = 11, fontface = "bold", col = "grey40"))
  popViewport()
  pushViewport(viewport(x = 0.92, y = 0.08, width = 0.05, height = 0.04))
  grid.text("3'", gp = gpar(fontsize = 11, fontface = "bold", col = "grey40"))
  popViewport()

  # ---- legend (grid-based, top-right corner) ----
  pushViewport(viewport(x = 0.90, y = 0.90, width = 0.10, height = 0.08,
                        just = c("centre","centre")))
  grid.rect(gp = gpar(fill = "white", col = "grey60", lwd = 0.5))
  # exon row
  grid.circle(x = 0.18, y = 0.72, r = unit(1.5, "mm"),
              gp = gpar(fill = region_cols["exon"], col = region_cols["exon"]))
  grid.text("exon",   x = 0.60, y = 0.72, gp = gpar(fontsize = 7), just = "centre")
  # intron row
  grid.circle(x = 0.18, y = 0.28, r = unit(1.5, "mm"),
              gp = gpar(fill = region_cols["intron"], col = region_cols["intron"]))
  grid.text("intron", x = 0.60, y = 0.28, gp = gpar(fontsize = 7), just = "centre")
  popViewport()

  plotted <- plotted + 1
}

dev.off()

if (plotted == 0) {
  message("WARNING: no genes had both exon and site data for trackViewer plot.")
  file.remove("${TV_PDF}")
}
message("trackViewer: plotted ", plotted, " genes")
RSCRIPT3
  [ -f "$TV_PDF" ] && echo "Step 7: trackViewer lollipop plot → $TV_PDF" >&2 \
                   || echo "WARNING: trackViewer plot skipped (library not available or no data)." >&2
else
  echo "Step 7: skipped trackViewer plot (no exon data)." >&2
fi
