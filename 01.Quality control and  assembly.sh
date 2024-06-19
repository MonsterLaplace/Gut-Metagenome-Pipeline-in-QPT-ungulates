# requirments (fastp, fastQC, multiQC, bwa-mem2, metaSPAdes)
#!/bin/sh
CONDA_BASE=$(conda info --base)
 source $CONDA_BASE/etc/profile.d/conda.sh

# 1. 指定物种
SP=hosts

# 2. 指定是谁的样品, 比如李斌的样品则为LB
BL=XB

# 3.指定原始rawData数据存放目录
OD=/mnt/z/xb/metagenome/00.rawdata/XB.hosts

# 指定宿主基因组存放路径
HG=/mnt/z/xb/metagenome/00.rawdata/XB.hosts.genome

# 指定宿主基因组名称, 普氏原羚用的是自己组装的基因组.
GF="hosts.fa"

# 4. 指定结果存放的根目录
WP=/mnt/z/xb/metagenome

# 5. metaSPAdes核心及内存控制
CPU=32
MEM=600

### 以下内容勿动 ###
###################
###################

# 01. 进行原始数据的质控

# 激活qcENV环境, 进行原始数据的质控处理
conda activate qcENV;

# 01.01 质控原始数据并生成报告

# 建立qc工作目录
LOG=logFile
mkdir -p ${WP}/01.QC/${SP}.${BL}.report.rawData/report.fastQC/${LOG} && \
mkdir -p ${WP}/01.QC/${SP}.${BL}.report.rawData/report.multiQC/${LOG} && \
mkdir -p ${WP}/01.QC/${SP}.${BL}.report.cleanData/report.fastp/${LOG} && \
mkdir -p ${WP}/01.QC/${SP}.${BL}.cleanData && \
mkdir -p ${WP}/01.QC/${SP}.${BL}.report.cleanData/report.fastQC/${LOG}  && \
mkdir -p ${WP}/01.QC/${SP}.${BL}.report.cleanData/report.multiQC/${LOG}

# 替换路径
fastQC=${WP}/01.QC/${SP}.${BL}.report.rawData/report.fastQC
multiQC=${WP}/01.QC/${SP}.${BL}.report.rawData/report.multiQC
fastP=${WP}/01.QC/${SP}.${BL}.report.cleanData/report.fastp
CLEAN=${WP}/01.QC/${SP}.${BL}.cleanData
DATE=$(date -d today +%Y%m%d)

# 生成原始数据的质控报告
filenames=$(for FILE in ${OD}/*
do
  if [ -f $FILE ]  # 如果是文件
  then
    # 使用 cut 命令和点号分隔符截取文件名
    filename=$(basename $FILE .raw.fastq.gz | cut -d. -f1-3)
    echo $filename
  fi
done | uniq)
wait
for FID in $filenames
do
	fastqc \
	${OD}/${FID}.*.gz \
	-t 6 \
	-o ${fastQC} \
	1>${fastQC}/${LOG}/${FID}.rawData.${DATE}.log \
	2>${fastQC}/${LOG}/${FID}.rawData.${DATE}.err &
done

# 使用fastp进行质控
for FID in $filenames
do
fastp \
	-i ${OD}/${FID}.R1.raw.fastq.gz \
	-I ${OD}/${FID}.R2.raw.fastq.gz \
	-o ${CLEAN}/${FID}.R1.clean.fastq \
	-O ${CLEAN}/${FID}.R2.clean.fastq \
	-W 5 \
	-M 20 \
	-5 \
	-q 15 \
	-u 40 \
	-n 0 \
	-l 75 \
	-w 4 \
	-j ${fastP}/${FID}.fastp.json \
	-h ${fastP}/${FID}.fastp.html \
	-R ${fastP}/${FID}.report \
	1>${fastP}/${LOG}/${FID}.${DATE}.log \
	2>${fastP}/${LOG}/${FID}.${DATE}.err &
done
wait

# 对原始数据的fastqc进行multiqc整合
multiqc ${fastQC} \
	-o ${multiQC} \
	1>${multiQC}/${LOG}/${SP}.${BL}.rawData.multiQC.${DATE}.log 2>${multiQC}/${LOG}/${SP}.${BL}.rawData.multiQC.${DATE}.err

fastQC=${WP}/01.QC/${SP}.${BL}.report.cleanData/report.fastQC
multiQC=${WP}/01.QC/${SP}.${BL}.report.cleanData/report.multiQC

# 对质控后的数据进行fastqc
for FID in $filenames
do
fastqc \
	${CLEAN}/${FID}.*.clean.fastq \
	-t 6 \
	-o ${fastQC} \
	1>${fastQC}/${LOG}/${FID}.cleanData.${DATE}.log \
	2>${fastQC}/${LOG}/${FID}.cleanData.${DATE}.err &
done
wait

# 对质控后fastqc报告进行multiqc整合
multiqc ${fastQC} \
	-o ${multiQC} \
	1>${multiQC}/${LOG}/${SP}.${BL}.cleanData.${DATE}.log \
	2>${multiQC}/${LOG}/${SP}.${BL}.cleanData.${DATE}.err

# 2. 构建宿主基因组bwa index文件
# 激活bwa环境
conda deactivate && \
conda activate bwa && \
mkdir -p ${HG}/${LOG}
cd ${HG};
bwa-mem2 index \
	-p ${SP} ${HG}/${SP}.fa \
	1>${HG}/${LOG}/buildBWAindex.${SP}.${BL}.${DATE}.log \
	2>${HG}/${LOG}/buildBWAindex.${SP}.${BL}.${DATE}.err &&

# 3. bwa比对,剔除宿主基因组污染
mkdir -p ${WP}/02.bwa.mapping/02.01.bwa.samFile/ && \
mkdir -p ${WP}/02.bwa.mapping/02.02.bwa.sortedBamFile/ && \
BSAM=${WP}/02.bwa.mapping/02.01.bwa.samFile
BBAM=${WP}/02.bwa.mapping/02.02.bwa.sortedBamFile
for FID in $filenames
do
	bwa-mem2 mem \
	-t 200 \
	${HG}/${SP} \
	${CLEAN}/${FID}.R1.clean.fastq \
	${CLEAN}/${FID}.R2.clean.fastq > ${BSAM}/${FID}.bwa.sam &
done
wait
echo "bwa mapping process finished"

# 使用samtools将sam文件转换为bam文件, 并排序.
for FID in $filenames
do
	samtools view -@ 4 -bS ${BSAM}/${FID}.bwa.sam | samtools sort -@ 4  -o ${BBAM}/${FID}.sorted.bam &
done
wait

# 使用samtools对bam文件进行索引.
for FID in $filenames
do
	samtools index ${BBAM}/${FID}.sorted.bam &
done
wait
echo "sam2bam and sorted process finished"

# 删除无用的过程数据
rm -rf ${BSAM}/*.bwa.sam
# 使用samtools提取未必对上的所有序列, 用于后续宏基因组分析
mkdir -p ${WP}/02.bwa.mapping/02.03.hostFree/ && \
HF=${WP}/02.bwa.mapping/02.03.hostFree
for FID in $filenames
do
	samtools view -@ 2 -b -f 12 -F 256 -o ${HF}/${FID}.bothUnmapped.bam ${BBAM}/${FID}.sorted.bam &
	samtools view -@ 2 -b -f 4 -F 264 -o ${HF}/${FID}.R1Unmapped.bam ${BBAM}/${FID}.sorted.bam &
	samtools view -@ 2 -b -f 8 -F 260 -o ${HF}/${FID}.R2Unmapped.bam ${BBAM}/${FID}.sorted.bam &
done
wait
echo "extract unmapped data finished, processing merge"

# 使用samtools对未比对上的所有序列进行合并, 排序, 并提取fastq序列
for FID in $filenames
do
	samtools merge -@ 4 ${HF}/${FID}.hostFree.bam ${HF}/${FID}.bothUnmapped.bam ${HF}/${FID}.R1Unmapped.bam ${HF}/${FID}.R2Unmapped.bam &
done
wait
for FID in $filenames
do
	samtools sort -@ 4  -o ${HF}/${FID}.hostFree.sorted.bam ${HF}/${FID}.hostFree.bam &
done
wait

# 删除无用的过程数据
rm -rf ${HF}/*.bothUnmapped.bam ${HF}/*.R1Unmapped.bam ${HF}/*.R2Unmapped.bam
使用samtools对bam文件建立索引
for FID in $filenames
do
	samtools index ${HF}/${FID}.hostFree.sorted.bam &
done
wait

# 使用bedtools从bam文件中提取fastq序列
mkdir -p ${WP}/03.hostFreeData/
HFFAST=${WP}/03.hostFreeData
for FID in $filenames
do
	bamToFastq \
	-i ${HF}/${FID}.hostFree.sorted.bam \
	-fq ${HFFAST}/${FID}.hostFree.R1.fq \
	-fq2 ${HFFAST}/${FID}.hostFree.R2.fq &
done
wait
echo "extract hostFree fq file finished"

提取能够比对上宿主的fq序列, 如果需要, 请将以下注释删除.
####################################################
mkdir -p ${WP}/04.onlyHostData/
OH=${WP}/04.onlyHostData
for FID in $filenames
do
	samtools view -@ 4 -b -f 1 -F 12 ${BBAM}/${FID}.sorted.bam | samtools sort -@ 4  -o ${OH}/${FID}.onlyMapToHost.sorted.bam &
done
wait

# 使用samtools对bam文件建立索引
for FID in $filenames
do
	samtools index ${OH}/${FID}.onlyMapToHost.sorted.bam &
done
wait

# 使用bedtools从bam文件中提取fastq序列
for FID in $filenames
do
	bamToFastq \
	-i ${OH}/${FID}.onlyMapToHost.sorted.bam \
	-fq ${OH}/${FID}.onlyMapToHost.R1.fq \
	-fq2 ${OH}/${FID}.onlyMapToHost.R2.fq &
done
wait
echo "extract onlyMapToHost fq file finished"
#####################################################
conda deactivate

# 4. 使用metaSPAdes基于上述hostFree reads组装宏基因组
conda activate metaSPAdes &&
mkdir -p ${WP}/05.assembly.metaSPAdes/logFile
ASP=${WP}/05.assembly.metaSPAdes
# 一次最多仅拼接8个样品, 在这里需要控制数量
# 控制同时运行的进程数
MAX_PROCESSES=8
CURRENT_PROCESSES=0

for FID in $filenames; do
  # 检查输出文件夹是否已存在
  if [[ ! -d "${ASP}/${FID}" ]]; then
    mkdir -p "${ASP}/${FID}"
  fi

  # 启动后台进程
  metaspades.py \
    -1 "${HFFAST}/${FID}.hostFree.R1.fq" \
    -2 "${HFFAST}/${FID}.hostFree.R2.fq" \
    -o "${ASP}/${FID}/" \
    -m ${MEM} \
    -t ${CPU} &
  CURRENT_PROCESSES=$((CURRENT_PROCESSES + 1))

  # 判断是否达到最大并发数，如果是则等待所有进程完成
  if [[ "$CURRENT_PROCESSES" -ge "$MAX_PROCESSES" ]]; then
    wait
    CURRENT_PROCESSES=0
  fi
done

# 等待所有进程完成
wait
# 使用metaSPAdes组装宏基因组MAGs完成.
echo "all process done, lol!"