#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l m_mem_free=8G
#$ -pe smp 4
#$ -m e
#$ -M masonwang@neb.com
#$ -V

eval "$(micromamba shell hook --shell bash)"
micromamba activate base_env
set -ue

NUMCPU=4

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 <ref_fasta> <gtf> <sites_bed> <outdir>"
  exit 1
fi

REF_FASTA="$1"
GTF="$2"
SITES_BED="$3"
OUTDIR="$4"

mkdir -p "$OUTDIR"

# STEP 0: genome.txt (FASTA record order) - all downstream sorts use this via -g,
# not plain `sort`, so bedtools never hits chromosome-order mismatches
[ -s "${REF_FASTA}.fai" ] || samtools faidx "$REF_FASTA"
cut -f1,2 "${REF_FASTA}.fai" > "${OUTDIR}/genome.txt"

# STEP 1: find all genomic DRACH motifs, strand-aware, write as BED
python - "$REF_FASTA" "$OUTDIR" <<'EOF'
import sys, re
from Bio import SeqIO

fasta, outdir = sys.argv[1], sys.argv[2]
# DRACH: D=[AGT] R=[AG] A C H=[ACT], modified A is at motif position 3 (0-indexed 2)
pattern = re.compile(r'(?=([AGT][AG]AC[ACT]))', re.IGNORECASE)

out_bed = f"{outdir}/all_drach_motifs.bed"
n_plus = 0
n_minus = 0

with open(out_bed, "w") as out:
    for record in SeqIO.parse(fasta, "fasta"):
        chrom = record.id
        seq = str(record.seq).upper()
        seqlen = len(seq)

        # plus strand: motif start = match position, A (modified base) at match+2
        for m in pattern.finditer(seq):
            motif_start = m.start()
            a_pos = motif_start + 2   # 0-based position of the modified A
            out.write(f"{chrom}\t{a_pos}\t{a_pos+1}\t.\t.\t+\n")
            n_plus += 1

        # minus strand: scan the reverse complement, then map coordinates back
        # to the original forward-strand coordinate system
        rev = str(record.seq.reverse_complement()).upper()
        for m in pattern.finditer(rev):
            motif_start = m.start()
            a_pos_rev = motif_start + 2
            a_pos_fwd = seqlen - 1 - a_pos_rev
            out.write(f"{chrom}\t{a_pos_fwd}\t{a_pos_fwd+1}\t.\t.\t-\n")
            n_minus += 1

print(f"[drach] plus-strand motifs: {n_plus}", file=sys.stderr)
print(f"[drach] minus-strand motifs: {n_minus}", file=sys.stderr)
print(f"[drach] total: {n_plus + n_minus}", file=sys.stderr)
EOF

bedtools sort -g "${OUTDIR}/genome.txt" -i "${OUTDIR}/all_drach_motifs.bed" \
  > "${OUTDIR}/all_drach_motifs.sorted.bed"
mv "${OUTDIR}/all_drach_motifs.sorted.bed" "${OUTDIR}/all_drach_motifs.bed"
N_ALL_DRACH_RAW=$(wc -l < "${OUTDIR}/all_drach_motifs.bed")
echo "[drach] total genomic DRACH motifs, both strands, unfiltered: $N_ALL_DRACH_RAW"

# STEP 2: build region BEDs (exon, intron-by-subtraction, unannotated)
echo "[regions] extracting exon and gene features from GTF..."

grep -v '^#' "$GTF" | awk -F'\t' '$3=="exon"' | \
  awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5, ".", ".", $7}' | \
  bedtools sort -g "${OUTDIR}/genome.txt" -i - | \
  bedtools merge -s -i - -c 6 -o distinct | \
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,".",".",$4}' > "${OUTDIR}/exons.merged.bed"

grep -v '^#' "$GTF" | awk -F'\t' '$3=="gene"' | \
  awk 'BEGIN{OFS="\t"}{print $1, $4-1, $5, ".", ".", $7}' | \
  bedtools sort -g "${OUTDIR}/genome.txt" -i - > "${OUTDIR}/genes.bed"

# introns = gene span minus exons (strand-aware subtraction)
bedtools subtract -s -a "${OUTDIR}/genes.bed" -b "${OUTDIR}/exons.merged.bed" \
  | bedtools sort -g "${OUTDIR}/genome.txt" -i - > "${OUTDIR}/introns.bed"

# unannotated = whole genome minus (union of all gene spans, either strand)
awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' "${OUTDIR}/genes.bed" | \
  bedtools sort -g "${OUTDIR}/genome.txt" -i - | \
  bedtools merge -i - > "${OUTDIR}/genes.merged.noStrand.bed"

bedtools complement -i "${OUTDIR}/genes.merged.noStrand.bed" -g "${OUTDIR}/genome.txt" \
  > "${OUTDIR}/unannotated.bed"

N_EXON_BP=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/exons.merged.bed")
N_INTRON_BP=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/introns.bed")
N_UNANN_BP=$(awk '{sum+=$3-$2} END{print sum+0}' "${OUTDIR}/unannotated.bed")
echo "[regions] exon bp: $N_EXON_BP | intron bp: $N_INTRON_BP | unannotated bp: $N_UNANN_BP"

# STEP 3: exclude antisense (non-transcribed-strand) DRACH motifs.
# dRNA-seq only ever reads the transcribed strand, so a DRACH motif on the
# opposite strand of a gene body is genomic background, not a real candidate
# site. -S (capital) finds antisense-overlapping motifs; we subtract those
# from the master set. Motifs outside any gene body are left untouched since
# there's no annotation to judge sense vs antisense against.
echo "[strand] identifying antisense (non-transcribed-strand) DRACH motifs within genes..."

bedtools intersect -S -u -a "${OUTDIR}/all_drach_motifs.bed" -b "${OUTDIR}/genes.bed" \
  > "${OUTDIR}/antisense_drach_motifs.bed"
N_ANTISENSE=$(wc -l < "${OUTDIR}/antisense_drach_motifs.bed")
echo "[strand] antisense DRACH motifs excluded: $N_ANTISENSE"

# transcribable DRACH universe = all motifs minus antisense-in-gene motifs.
# -s is critical: without it, a sense-strand motif sharing coordinates with
# an antisense motif at the same locus would be incorrectly removed.
bedtools intersect -s -v -a "${OUTDIR}/all_drach_motifs.bed" -b "${OUTDIR}/antisense_drach_motifs.bed" \
  > "${OUTDIR}/all_drach_motifs.transcribable.bed"
N_ALL_DRACH=$(wc -l < "${OUTDIR}/all_drach_motifs.transcribable.bed")
echo "[strand] transcribable DRACH universe (antisense excluded): $N_ALL_DRACH"

DRACH_BED="${OUTDIR}/all_drach_motifs.transcribable.bed"

# STEP 4: intersect transcribable DRACH motifs with each region
# (-s for exon/intron: strand must match the annotation; unannotated has no
# strand to match against)
echo "[intersect] assigning DRACH motifs to regions..."

bedtools intersect -s -u -a "$DRACH_BED" -b "${OUTDIR}/exons.merged.bed" \
  > "${OUTDIR}/drach_in_exon.bed"
bedtools intersect -s -u -a "$DRACH_BED" -b "${OUTDIR}/introns.bed" \
  > "${OUTDIR}/drach_in_intron.bed"
bedtools intersect -u -a "$DRACH_BED" -b "${OUTDIR}/unannotated.bed" \
  > "${OUTDIR}/drach_in_unannotated.bed"

N_DRACH_EXON=$(wc -l < "${OUTDIR}/drach_in_exon.bed")
N_DRACH_INTRON=$(wc -l < "${OUTDIR}/drach_in_intron.bed")
N_DRACH_UNANN=$(wc -l < "${OUTDIR}/drach_in_unannotated.bed")
echo "[intersect] DRACH in exon: $N_DRACH_EXON | intron: $N_DRACH_INTRON | unannotated: $N_DRACH_UNANN"

# STEP 5: intersect modified sites (called m6A) with each region. These come
# from actual Nanopore reads, which are inherently strand-specific already,
# so no antisense-exclusion step is needed - only region assignment with
# strand-matching (-s) for exon/intron.
echo "[intersect] assigning modified sites to regions..."

bedtools intersect -s -u -a "$SITES_BED" -b "${OUTDIR}/exons.merged.bed" \
  > "${OUTDIR}/modsites_in_exon.bed"
bedtools intersect -s -u -a "$SITES_BED" -b "${OUTDIR}/introns.bed" \
  > "${OUTDIR}/modsites_in_intron.bed"
bedtools intersect -u -a "$SITES_BED" -b "${OUTDIR}/unannotated.bed" \
  > "${OUTDIR}/modsites_in_unannotated.bed"

N_MOD_EXON=$(wc -l < "${OUTDIR}/modsites_in_exon.bed")
N_MOD_INTRON=$(wc -l < "${OUTDIR}/modsites_in_intron.bed")
N_MOD_UNANN=$(wc -l < "${OUTDIR}/modsites_in_unannotated.bed")
echo "[intersect] modified sites in exon: $N_MOD_EXON | intron: $N_MOD_INTRON | unannotated: $N_MOD_UNANN"

# STEP 6: summary table - base length per region, DRACH site density
# (transcribable motifs/kb), modified site density (modified sites/kb),
# and modified/potential ratio (modified sites / transcribable motifs)
SUMMARY="${OUTDIR}/region_summary.csv"

{
  echo "region,region_bp,n_drach_motifs,drach_density_per_kb,n_modified_sites,modified_density_per_kb,modified_over_potential_pct"

  awk -v bp="$N_EXON_BP" -v ndrach="$N_DRACH_EXON" -v nmod="$N_MOD_EXON" 'BEGIN{
    drach_dens = (bp>0) ? ndrach/(bp/1000) : 0
    mod_dens   = (bp>0) ? nmod/(bp/1000)   : 0
    ratio      = (ndrach>0) ? (nmod/ndrach)*100 : 0
    printf "exon,%d,%d,%.4f,%d,%.4f,%.4f\n", bp, ndrach, drach_dens, nmod, mod_dens, ratio
  }'

  awk -v bp="$N_INTRON_BP" -v ndrach="$N_DRACH_INTRON" -v nmod="$N_MOD_INTRON" 'BEGIN{
    drach_dens = (bp>0) ? ndrach/(bp/1000) : 0
    mod_dens   = (bp>0) ? nmod/(bp/1000)   : 0
    ratio      = (ndrach>0) ? (nmod/ndrach)*100 : 0
    printf "intron,%d,%d,%.4f,%d,%.4f,%.4f\n", bp, ndrach, drach_dens, nmod, mod_dens, ratio
  }'

  awk -v bp="$N_UNANN_BP" -v ndrach="$N_DRACH_UNANN" -v nmod="$N_MOD_UNANN" 'BEGIN{
    drach_dens = (bp>0) ? ndrach/(bp/1000) : 0
    mod_dens   = (bp>0) ? nmod/(bp/1000)   : 0
    ratio      = (ndrach>0) ? (nmod/ndrach)*100 : 0
    printf "unannotated,%d,%d,%.4f,%d,%.4f,%.4f\n", bp, ndrach, drach_dens, nmod, mod_dens, ratio
  }'
} > "$SUMMARY"

echo ""
echo "=== SUMMARY ==="
column -s',' -t "$SUMMARY"
echo ""
echo "wrote: $SUMMARY"
echo "DONE: $(date)"
