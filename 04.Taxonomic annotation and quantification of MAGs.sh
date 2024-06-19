# 1. Taxonomic annotation of strain-level MAGs

mkdir mkdir /mnt/z/xb/metagenome/finalbins/allsamples/MAGclass
time gtdbtk classify_wf \
    --cpus 200 \
    --pplacer_cpus 64 \
    -x fa \
    --genome_dir /mnt/z/xb/metagenome/finalbins/allsamples/ \
    --out_dir /mnt/z/xb/metagenome/finalbins/allsamples/MAGclass \
    --skip_ani_screen && echo "**MAGS Taxonomy  Annoation done **"
	

# 2. The abundances of MAGs in each contig

for i in hostsname; do
outdir=/mnt/z/xb/metagenome/finalbins/${i}/metawrap
for b in samples ; do
metaWRAP quant_bins \
	-t 120 \
	-b ${outdir}/bins/${b} \
	-o ${outdir}/${i}.${b} \
	-a /mnt/z/xb/metagenome/finalbins/${i}/L1000contig/${i}.${b}.contig.fasta /mnt/z/xb/metagenome/finalbins/${i}/hostfreefq/${i}.${b}.hostFree_1.fastq /mnt/z/xb/metagenome/finalbins/${i}/hostfreefq/${i}.${b}.hostFree_2.fastq
done
done
# -a 输入每个物种每个样本的过滤后的contig和每个样本的去宿主污染过后的原始双端宏基因组reads.

# 3. The correlation (Spearman's r) between the abundance and completeness of MAGs 

library(ggplot2)
library(ggpubr)
library(ggExtra)

td <-read.table("MAGS_abundances_vs_completeness.txt", header=T)
p1=ggplot(td, aes(abundances,completeness)) + 
  xlab("abundances")+ylab("completeness")+
  geom_point(shape = 21, colour = "red3", fill = "lightgray", size = 1.5, stroke = .5,alpha=0.8)+ geom_smooth(method="lm",formula = y ~ x,linetype=2,color="black",fill="lightgray") + theme_bw()+
  stat_cor(method = 'spearman', aes(x =abundances, y =completeness))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "gray"),yparams = list(fill = "gray"))
p2


