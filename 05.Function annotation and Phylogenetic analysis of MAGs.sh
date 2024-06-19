# 1. Function annotation of MAGs
for b in hostsample; do
for i in /mnt/z/xb/metagenome/finalbins/${b}/straindRep/dereplicated_genomes/*.fa; do
cd /mnt/z/xb/metagenome/finalbins/${b}
file=${i##*/}
base=${file%.fa}
bakta ${i} \
	--db /mnt/z/xb/metagenome/database/db \
	--meta --threads 350 --prefix $base --skip-trn --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-crispr --skip-pseudo --skip-sorf --skip-gap --skip-ori --skip-plot \
	--output /mnt/z/xb/metagenome/finalbins/${b}/$base
done
done

# 2. CAZYmes in the MAGs
mkdir /mnt/z/xb/metagenome/finalbins/faa
cp -r  /mnt/z/xb/metagenome/finalbins/*/$base/*.faa  /mnt/z/xb/metagenome/finalbins/faa
for i in /mnt/z/xb/metagenome/finalbins/faa/*.faa; do
file=${i##*/}
base=${file%.faa}
run_dbcan  \
	--db_dir /mnt/z/xb/metagenome/database/dbCAN/db \
	--hmm_cov 0.35 --hmm_eval 1e-15 --hmm_cpu 20 \
	--dia_eval 1e-102 --dia_cpu 20 \
	--out_dir /mnt/z/xb/metagenome/finalbins/faa/$base --out_pre $base ${i} protein
done

## Sumary of the CAZYmes

mkdir /mnt/z/xb/metagenome/SV/ALLbins
mkdir /mnt/z/xb/metagenome/SV/ALLbins/CAZY/
mkdir /mnt/z/xb/metagenome/SV/ALLbins/CAZY/CYAZmes
cp -r /mnt/z/xb/metagenome/finalbins/allsamples/SGBdRep /mnt/z/xb/metagenome/SV/ALLbins/
for i in $(ls /mnt/z/xb/metagenome/finalbins/faa/*.faa);do
m=${i##*/}
n=${m%.faa}
cd /mnt/z/xb/metagenome/SV/ALLbins/bins
volum=`seqkit stat -T /mnt/z/xb/metagenome/SV/ALLbins/bins/${n}.fa | sed '1d' |cut -f 5`
cd /mnt/z/xb/metagenome/SV/ALLbins/CAZY/CYAZmes
cp -r /mnt/z/xb/metagenome/finalbins/faa/*/*.overview.txt /mnt/z/xb/metagenome/SV/ALLbins/CAZY/CYAZmes
num1=`wc -l ${n}overview.txt | awk '{print$1-1}'`
echo -e "${n}\t${volum}\t${num1}" >> /mnt/z/xb/metagenome/SV/ALLbins/CYAZ/0.stat.tsv
done &
for i in $(ls /mnt/z/xb/metagenome/finalbins/faa/*.faa);do
m=${i##*/}
n=${m%.faa}
cd /mnt/z/xb/metagenome/SV/ALLbins/CYAZ/CYAZmes
for b in $(echo `sed '1d' ${n}overview.txt | awk '{print $1}' `) ; do
cd /mnt/z/xb/metagenome/SV/ALLbins/CYAZ/binstsv
num2=`sed '1d' ${n}.gff3 |grep ${b} | awk '{print $1}'`
num3=`sed '1d' ${n}.gff3 |grep ${b} | awk '{print $4}'`
num4=`sed '1d' ${n}.gff3 |grep ${b} | awk '{print $5}'`
echo -e "${n}\t$num2\t$num3\t$num4" >> /mnt/z/xb/metagenome/SV/ALLbins/CYAZ/1.stat.tsv
done
done

# 3. Polygenetic analysis of MAGs
# method-1
mkdir /mnt/z/xb/metagenome/finalbins/allsamples/phylophlantree
phylophlan \
     --input_folder /mnt/z/xb/metagenome/finalbins/allsamples/SGBdRep/dereplicated_genomes \
     --output_folder /mnt/z/xb/metagenome/finalbins/allsamples/phylophlantree \
     -d phylophlan \
     --diversity high \
     -f /mnt/z/xb/metagenome/finalbins/allsamples/phylophlan/faatree.cfg \
     --nproc 250 \
     -t a \
    --fast \
    --force_nucleotides \
    --genome_extension ".fa" \
     -i wgt 2>&1 | tee phylophlan.log


# method-2
mkdir /mnt/z/xb/metagenome/finalbins/allsamples/GTDBTKtree
gtdbtk classify_wf \
	--genome_dir /mnt/z/xb/metagenome/finalbins/allsamples/SGBdRep/dereplicated_genomes \
	--out_dir mnt/z/xb/metagenome/finalbins/allsamples/GTDBTKtree \
	--extension .fa \
	--cpus 200

# 4. The primary metabolic gene clusters (MGCs) in SGBs 
for i in /mnt/z/xb/metagenome/finalbins/allsamples/SGBdRep/dereplicated_genomes/*.fa
python3 gutsmash/run_gutsmash.py --minimal --cb-knownclusters --enable-genefunctions --genefinding-tool prodigal -c 36  ${i}
done