# 1. SVs prediction

mkdir /mnt/z/xb/metagenome/SV/
mkdir /mnt/z/xb/metagenome/SV/MUMandCo
OUTPUT_DIR="/mnt/z/xb/metagenome/SV/MUMandCo"

# filter the clusters present in >1 samples 

for i in `sed '1d' /mnt/z/xb/metagenome/SV/ALLbins/SGBdRep/data_tables/Cdb.csv |cut -f2 -d',' |sort |uniq -c | sort -nr | awk '($1>1) {print}' |awk '{print $2}'`;do
  mkdir $OUTPUT_DIR/${i}
  for j in `sed '1d' /mnt/z/xb/metagenome/SV/ALLbins/SGBdRep/data_tables/Cdb.csv | awk -v var="${i}" 'BEGIN {FS=","} ($2==var) {print $1}'`;do
    cp /mnt/z/xb/metagenome/SV/ALLbins/bins/$j $OUTPUT_DIR/${i}
  done
  ref=`sed '1d' /mnt/z/xb/metagenome/SV/ALLbins/SGBdRep/data_tables/Widb.csv |awk -v var="${i}" 'BEGIN {FS=","} ($8==var) {print $1}'`
  rm $OUTPUT_DIR/${i}/${ref}
  for m in $(ls $OUTPUT_DIR/${i}/*fa);do 
    m=`echo ${m} |cut -d'/' -f 9`
    gsize=`seqkit stat -T /mnt/z/xb/metagenome/SV/ALLbins/bins/${ref} | sed '1d' |cut -f 5`
    cd $OUTPUT_DIR/${i}/
    sh /mnt/z/xb/metagenome/SV/mumandco_v3.8.sh -r /mnt/z/xb/metagenome/SV/ALLbins/bins/${ref} \
    -q $OUTPUT_DIR/${i}/${m} -g ${gsize} -o ${ref}-${m}- -t 24
wait
    cd /mnt/z/xb/metagenome/SV/MUMandCo
  done
  rm $OUTPUT_DIR/${i}/*fa
wait
done

# 2. SVs extraction
#plateau=BOG BOM EQK PRP PSN CEA CEC GAS OVAM
#plain=BOT CAP CAH CEN EQC OVA
mkdir /mnt/z/xb/metagenome/SV/SVfilter/
mkdir /mnt/z/xb/metagenome/SV/SVfilter/plateau
mkdir /mnt/z/xb/metagenome/SV/SVfilter/plain

for m in BOG BOM EQK PRP PSN CEA CEC GAS OVAM ; do
for n in BOT CAP CAH CEN EQC OVA; do
for i in $(ls /mnt/z/xb/metagenome/SV/MUMandCo/*/*_output/${m}.*.fa-${n}.*.fa-.SVs_all.withfragment.tsv);do
cp -r ${i} /mnt/z/xb/metagenome/SV/SVfilter/plain
done
done
done
wait
for m in BOT CAP CAH CEN EQC OVA ; do
for n in BOG BOM EQK PRP PSN CEA CEC GAS OVAM; do
for j in $(ls /mnt/z/xb/metagenome/SV/MUMandCo/*/*_output/${m}.*.fa-${n}.*.fa-.SVs_all.withfragment.tsv);do
cp -r ${j} /mnt/z/xb/metagenome/SV/SVfilter/plateau
done
done
done

# 3. statistic of numbers of predicted SVs
cat /mnt/z/xb/metagenome/SV/SVfilter/plateau/*.tsv /mnt/z/xb/metagenome/SV/SVfilter/plain/*.tsv > /mnt/z/xb/metagenome/SV/SVfilter/all.tsv

ref=`echo /mnt/z/xb/metagenome/SV/SVfilter/all.tsv |cut -d'/' -f8 |sed 's/_output//' |cut -d'-' -f1`
query=`echo /mnt/z/xb/metagenome/SV/SVfilter/all.tsv |cut -d'/' -f8 |sed 's/_output//' |cut -d'-' -f2`

num1=`sed '1d' /mnt/z/xb/metagenome/SV/SVfilter/all.tsv |grep "deletion"  | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

num2=`sed '1d' /mnt/z/xb/metagenome/SV/SVfilter/all.tsv |grep "insertion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3}' |wc -l `

num3=`sed '1d' /mnt/z/xb/metagenome/SV/SVfilter/all.tsv |grep "duplication" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l `

num4=`sed '1d' /mnt/z/xb/metagenome/SV/SVfilter/all.tsv |grep transloc | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`
  
num5=`sed '1d' /mnt/z/xb/metagenome/SV/SVfilter/all.tsv |grep "inversion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 |wc -l`

total=`expr $num1 + $num2 + $num3 + $num4 + $num5`

echo -e "$ref\t$query\t$num1\t$num2\t$num3\t$num4\t$num5\t$total" >>/mnt/z/xb/metagenome/SV/SVsummary/MUMandCo_SVs_num_all.summary

sed -i '1 i ref\tquery\tdeletion_num\tinsertion_num\tduplication_num\ttransloc_num\tinversion_num\ttotal' /mnt/z/xb/metagenome/SV/SVsummary/MUMandCo_SVs_num_all.summary

# 4. statistic of length of all predicted SVs
cat /mnt/z/xb/metagenome/SV/SVfilter/all.tsv | grep "deletion" |awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 >>/mnt/z/xb/metagenome/SV/SVsummary/MUMandCo_SVs_num_all.summary
cat /mnt/z/xb/metagenome/SV/SVfilter/all.tsv | grep "insertion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3}' >>/mnt/z/xb/metagenome/SV/SVsummary/MUMandCo_SVs_num_all.summary
cat /mnt/z/xb/metagenome/SV/SVfilter/all.tsv | grep "inversion" | awk '($3>10) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6 >>/mnt/z/xb/metagenome/SV/SVsummary/MUMandCo_SVs_num_all.summary
sed -i '1 i contigs\tstart\tend\tlength\ttype' /mnt/z/xb/metagenome/SV/SVsummary/MUMandCo_SVs_num_all.summary


# 5. SVs fillting
cd /mnt/z/xb/metagenome/SV/SVfilter/plateau/
for i in *.SVs_all.withfragment.tsv ;do
m=`echo ${i} |cut -d'-' -f 1`
n=`echo ${i} |cut -d'-' -f 2`
mkdir /mnt/z/xb/metagenome/SV/SVfilter/deletion/${m}-${n}

cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "deletion" |awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6,9 >>/mnt/z/xb/metagenome/SV/SVfilter/deletion/${m}-${n}/${m}-${n}.deletion.tsv

cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "deletion" | awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f1,3,4,5,6,9 |awk '{print ">"$1"\n"$6"\n"}' >>/mnt/z/xb/metagenome/SV/SVfilter/deletion/${m}-${n}/${m}.d.fa

mkdir /mnt/z/xb/metagenome/SV/SVfilter/insertion/${m}-${n}

cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "insertion" |awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f2,5,6,7,8,9 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$6}' >>/mnt/z/xb/metagenome/SV/SVfilter/insertion/${m}-${n}/${m}-${n}.insertion.tsv

cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "insertion" | awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8,9 |awk '{print ">"$1"\n"$6"\n"}' >>/mnt/z/xb/metagenome/SV/SVfilter/insertion/${m}-${n}/${n}.i.fa

done

wait

cd /mnt/z/xb/metagenome/SV/SVfilter/plain/
for i in *.SVs_all.withfragment.tsv ;do
m=`echo ${i} |cut -d'-' -f 1`
n=`echo ${i} |cut -d'-' -f 2`
mkdir /mnt/z/xb/metagenome/SV/SVfilter/deletion/${m}-${n}
cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "insertion" |awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f2,5,6,7,8,9 |awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$6}' >>/mnt/z/xb/metagenome/SV/SVfilter/deletion/${m}-${n}/${m}-${n}.insertion.tsv

cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "insertion" | awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f2,5,6,7,8,9 |awk '{print ">"$1"\n"$6"\n"}' >>/mnt/z/xb/metagenome/SV/SVfilter/deletion/${m}-${n}/${n}.i.fa

mkdir /mnt/z/xb/metagenome/SV/SVfilter/insertion/${m}-${n}

cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "deletion" |awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' |cut -f1,3,4,5,6,9 >>/mnt/z/xb/metagenome/SV/SVfilter/insertion/${m}-${n}/${m}-${n}.deletion.tsv

cat ${m}-${n}-.SVs_all.withfragment.tsv | grep "deletion" | awk '($3>10) {print}' |awk '($5>200) {print}' |awk '($7>10) {print}' | awk '{split($1,a,"_");if($4<a[4]-10){print}}' | \
    awk '{split($2,a,"_");if($8<a[4]-10){print}}' | cut -f1,3,4,5,6,9 |awk '{print ">"$1"\n"$6"\n"}' >>/mnt/z/xb/metagenome/SV/SVfilter/insertion/${m}-${n}/${m}.d.fa

done

# 6. Prediction the genes in SVs
mkdir /mnt/z/xb/metagenome/SV/SVfilter/deletion/annotation
for i in /mnt/z/xb/metagenome/SV/SVfilter/deletion/*/*.fa; do
file=${i##*/}
base=${file%.fa}
bakta $i \
	--db /mnt/z/xb/metagenome/database/db \
	--meta \
	--threads 200 \
	--prefix $base \
	--output /mnt/z/xb/metagenome/SV/SVfilter/deletion/annotation/$base --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-ori --skip-plot
echo -e "bakta annotation Done"
done

mkdir /mnt/z/xb/metagenome/SV/SVfilter/insertion/annotation
for i in /mnt/z/xb/metagenome/SV/SVfilter/insertion/*/*.fa; do
file=${i##*/}
base=${file%.fa}
bakta $i \
	--db /mnt/z/xb/metagenome/database/db \
	--meta \
	--threads 200 \
	--prefix $base \
	--output /mnt/z/xb/metagenome/SV/SVfilter/insertion/annotation/$base --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-ori --skip-plot
echo -e "bakta annotation Done"
done

# 7. The CYAZmes overlapping with SVs
preparing the coordinate of CAZYmes according to the gff file generated by the bakta
preparing the coordinate of SVs according to the tsv file generated by the mumandco_v3.8.sh
bedtools intersect -a CAZYmes.bed -b SVs.bed -wb > SVs_Overlap.tsv