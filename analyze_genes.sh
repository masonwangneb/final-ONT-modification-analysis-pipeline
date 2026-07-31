#!/bin/bash
# analyze_genes.sh
#
# Per-gene modification analysis:
# -Assigns each modified site to a gene
# -Computes number of modified sites and mean percent modification (for genes with >0 modified sites, sites with >10 coverage)
# -Creates csv table + scatter plot of modified genes
#
set -uo pipefail

# ---- inputs
BEDM="/path/to/modkit.bedmethyl"    # modkit bedmethyl (any codes; filtered by MOD_CODE)
GTF="/path/to/annotation.gtf"       # GTF/GFF3 with 'gene' features and gene_id attribute
CODE="a"                            # modkit code in column 4 (a, 17596, 17802, m, ...)
TAG="mytag"                         # label for outputs (e.g. HG002_m6A, mosquito_inosine)
OUTDIR="/path/to/outdir"            # output directory
MINCOV=10                           # Nvalid_cov must be > this
TOPN=30                             # number of genes to label on the plot
CHROM_FILTER=""                     # only sites on matching chromosomes (e.g. '^NC_'), empty = no filter

# Requires: bedtools, awk, sort, Rscript (with ggplot2, ggrepel optional)

[ -s "$BEDM" ] || { echo "ERROR: bedmethyl not found: $BEDM" >&2; exit 1; }
[ -s "$GTF" ]  || { echo "ERROR: annotation not found: $GTF" >&2; exit 1; }
for cmd in bedtools awk sort Rscript; do
  command -v "$cmd" >/dev/null || { echo "ERROR: '$cmd' not on PATH" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

# -------- extract gene BED with gene_id + gene_name from GTF (works for
# NCBI-style and Ensembl-style attribute quoting)
GENE_BED="${OUTDIR}/${TAG}_genes.bed"
echo "[genes] parsing gene features from GTF..."
grep -v '^#' "$GTF" | awk -F'\t' '$3=="gene"' | awk 'BEGIN{OFS="\t"}
  {
    gid=""; gname=""
    if (match($0, /gene_id "[^"]+"/))   gid   = substr($0, RSTART+9,  RLENGTH-10)
    if (match($0, /gene_name "[^"]+"/)) gname = substr($0, RSTART+11, RLENGTH-12)
    if (gname=="") gname=gid
    print $1, $4-1, $5, gid"|"gname, ".", $7
  }' | sort -k1,1 -k2,2n > "$GENE_BED"
NGENES=$(wc -l < "$GENE_BED")
[ "$NGENES" -gt 0 ] || { echo "ERROR: no 'gene' features found in GTF" >&2; exit 1; }
echo "[genes] $NGENES gene features"

# -------- filter bedmethyl to modified + cov-passing sites for this mod code
SITES_BED="${OUTDIR}/${TAG}_sites.bed"
echo "[genes] filtering bedmethyl (mod=${CODE}, cov>${MINCOV}, pct>0)..."
awk -F'\t' -v c="$CODE" -v m="$MINCOV" \
    '$4==c && $11>0 && $10>m {print $1"\t"$2"\t"$3"\t"$11"\t"$10"\t"$6}' "$BEDM" \
  | ( [ -n "$CHROM_FILTER" ] && awk -F'\t' -v r="$CHROM_FILTER" '$1 ~ r' || cat ) \
  | sort -k1,1 -k2,2n > "$SITES_BED"
NSITES=$(wc -l < "$SITES_BED")
[ "$NSITES" -gt 0 ] || { echo "ERROR: 0 sites passed filter -- check CODE and MINCOV" >&2; exit 1; }
echo "[genes] ${NSITES} modified sites passing filter"

# -------- assign each site to its overlapping gene(s), strand-aware
# bedtools intersect: sites -a (BED6-ish), genes -b (BED6). '-wa -wb' keeps both.
# columns after intersect: 1-6 site (chr start end pct cov strand), 7-12 gene (chr start end gid|gname . strand)
ASSIGN="${OUTDIR}/${TAG}_site_gene_assignments.tsv"
bedtools intersect -s -a "$SITES_BED" -b "$GENE_BED" -wa -wb \
  | awk -F'\t' 'BEGIN{OFS="\t"} {print $10, $4, $5}' \
  > "$ASSIGN"
# columns: gene_id|gene_name, percent_modified, Nvalid_cov
NASSIGN=$(wc -l < "$ASSIGN")
echo "[genes] ${NASSIGN} site-to-gene assignments (sites in >1 gene appear >1x)"

# -------- aggregate per gene: site count, mean %modified, mean coverage
PER_GENE="${OUTDIR}/${TAG}_per_gene.csv"
{
  echo "gene_id,gene_name,n_sites,mean_pct_mod,median_pct_mod,max_pct_mod,mean_cov"
  awk -F'\t' '
    {
      split($1, a, "|"); gid=a[1]; gname=(a[2]?a[2]:gid)
      n[gid]++; name[gid]=gname
      sum_pct[gid]+=$2; sum_cov[gid]+=$3
      if ($2 > max_pct[gid]) max_pct[gid]=$2
      pct_list[gid] = (pct_list[gid]=="" ? $2 : pct_list[gid]" "$2)
    }
    END {
      for (g in n) {
        # median
        m = split(pct_list[g], v, " ")
        # simple sort
        for (i=1;i<=m;i++) for (j=i+1;j<=m;j++) if (v[i]>v[j]) { t=v[i]; v[i]=v[j]; v[j]=t }
        if (m%2) med=v[(m+1)/2]; else med=(v[m/2]+v[m/2+1])/2
        printf "%s,%s,%d,%.3f,%.3f,%.2f,%.2f\n", g, name[g], n[g], sum_pct[g]/n[g], med, max_pct[g], sum_cov[g]/n[g]
      }
    }
  ' "$ASSIGN" | sort -t',' -k3,3nr
} > "$PER_GENE"
NGENES_MOD=$(($(wc -l < "$PER_GENE") - 1))
echo "[genes] ${NGENES_MOD} genes have ≥1 modified site -> ${PER_GENE}"

# -------- plot: scatter of n_sites vs mean_pct_mod, top-N labeled
PLOT="${OUTDIR}/${TAG}_genes_scatter.png"
Rscript - <<RSCRIPT
suppressMessages({ library(ggplot2) })
have_repel <- requireNamespace("ggrepel", quietly = TRUE)

d <- read.csv("${PER_GENE}", stringsAsFactors = FALSE)
cat("genes in table:", nrow(d), "\n")
if (nrow(d) == 0) { message("no genes to plot"); quit(save="no") }

# label top-N by n_sites * mean_pct_mod (combined ranking)
d\$score <- d\$n_sites * d\$mean_pct_mod
d_top <- head(d[order(-d\$score),], ${TOPN})

p <- ggplot(d, aes(x = n_sites, y = mean_pct_mod)) +
  theme_bw() +
  geom_point(alpha = 0.4, size = 1) +
  scale_x_log10() +
  labs(x = "number of modified sites (log10)",
       y = "mean percent_modified per site",
       title = paste0("${TAG}: per-gene modification"),
       subtitle = paste0(nrow(d), " genes with >=1 modified site (cov > ${MINCOV})"))

if (have_repel) {
  p <- p + ggrepel::geom_text_repel(data = d_top,
    aes(label = gene_name), size = 3, max.overlaps = Inf,
    min.segment.length = 0, box.padding = 0.3)
} else {
  p <- p + geom_text(data = d_top, aes(label = gene_name),
    size = 3, vjust = -0.5, check_overlap = TRUE)
}

png("${PLOT}", width = 9, height = 7, units = "in", res = 300)
print(p)
dev.off()
cat("wrote plot: ${PLOT}\n")
RSCRIPT

echo ""
echo "outputs:"
echo "  ${PER_GENE}     (all genes, ranked by n_sites)"
echo "  ${PLOT}         (scatter with top ${TOPN} labeled)"
echo "  ${ASSIGN}       (raw site-gene assignments, for custom analysis)"
