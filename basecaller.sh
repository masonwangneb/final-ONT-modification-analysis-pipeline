#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l m_mem_free=8G
#$ -pe smp 4
#$ -l gpu=1
##$ -l chipset=intel
#$ -q all.q
#$ -m e
#$ -M masonwang@neb.com
#$ -V

################################################
set -ue;

if [ -z "${1:-}" ]; then
  echo "Usage: $0 <modification>"; exit 1;
fi

modification="$1";
echo "Running dorado + modkit for modification: ${modification}";

micromamba activate dorado_env;
################################################
NUMCPU=4;
################################################

which dorado;
echo "DORADO_MODELS_DIRECTORY = $DORADO_MODELS_DIRECTORY";

base_dir="YOUR_BASE_DIR";

in_pod5_dir="PATH_TO_YOUR_POD5";
out_dir="YOUR_OUT_DIR"

########--------------------------------------

[ -d $out_dir ] || mkdir $out_dir;

########--------------------------------------

# reference genome to map your data to
ref_genome="PATH_TO_YOUR_FASTA"
echo "ref_genome = $ref_genome";
ls -lh $ref_genome;

########--------------------------------------

# obtained from https://software-docs.nanoporetech.com/dorado/latest/models/list/
rna_model="YOUR_RNA_MODEL";

########--------------------------------------

my_out_bam="${out_dir}/polyA_mRNA_dorado_mrna_mods.refV7.${modification}.bam";
my_sorted_bam="${out_dir}/polyA_mRNA_dorado_mrna_mods.refV7.${modification}.sorted.bam";
my_out_bedMethyl="${out_dir}/polyA_mRNA_dorado_mrna_mods.refV7.${modification}.modkit.conf99.bedmethyl";

dorado basecaller \
  rna_model \
  $in_pod5_dir/ \
  --models-directory $DORADO_MODELS_DIRECTORY \
  --reference $ref_genome \
  --modified-bases $modification \
  --emit-moves  \
  --device cuda:0 \
  > $my_out_bam;

CMD="samtools sort --threads $NUMCPU -m 4G $my_out_bam -o $my_sorted_bam";
echo;echo "Running: $CMD [`date`]";eval ${CMD};

CMD="samtools index $my_sorted_bam";
echo;echo "Running: $CMD [`date`]";eval ${CMD};

modkit pileup \
	--mod-threshold a:0.99 --mod-threshold 17802:0.99 --mod-threshold 17596:0.99 \
	--mod-threshold 69426:0.99 --mod-threshold 19228:0.99 --mod-threshold m:0.99 \
	--mod-threshold 19229:0.99 --mod-threshold 19227:0.99 \
	--reference $ref_genome \
	$my_sorted_bam \
	$my_out_bedMethyl;

echo "DONE: `date`";
############### END OF SCRIPT #################################
