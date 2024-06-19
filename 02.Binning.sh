#!/bin/bash

# 1. Rename ContigFile

# 定义目标路径
target_dir="/mnt/p/spf/projectMeta/05.assembly.metaSPAdes"
out_dir="/mnt/p/spf/projectMeta/06.contigsFile/01.rawContigsFile"
log_dir="/mnt/p/spf/projectMeta/06.contigsFile/metaSPAdes.logFile"
missing_file="/mnt/p/spf/projectMeta/check.missingContigsFile.txt"

# 创建目标路径
mkdir -p "${out_dir}"
mkdir -p "${log_dir}"

# 清空可能存在的输出文件
> "${missing_file}"

# 进入目标路径
cd "${target_dir}"

# 为每个子目录执行操作
for dir in */; do
    # 如果存在contigs.fasta和spades.log文件
    if [ -f "${dir}contigs.fasta" ] && [ -f "${dir}spades.log" ]; then
        # 重命名contigs.fasta和spades.log文件
        new_contigs_name="$(basename ${dir}).metaSPAdes.contigs.fasta"
        new_log_name="$(basename ${dir}).metaSPAdes.log"
        mv "${dir}contigs.fasta" "${dir}${new_contigs_name}"
        mv "${dir}spades.log" "${dir}${new_log_name}"
        
        # 并行复制文件
        cp "${dir}${new_contigs_name}" "${out_dir}/${new_contigs_name}" &
        cp "${dir}${new_log_name}" "${log_dir}/${new_log_name}" &
    else
        # 如果不存在contigs.fasta文件，将样品名称输出到指定的文件中
        echo $(basename ${dir}) >> "${missing_file}"
    fi
done

# 等待所有后台任务完成
wait

# 2. Filter Contigs > 1000 bp

#!/bin/bash


conda_BASE=$(conda info --base)
source $conda_BASE/etc/profile.d/conda.sh


input_dir="/mnt/p/spf/projectMeta/06.contigsFile/01.rawContigsFile"
renamed_output_dir="/mnt/p/spf/projectMeta/06.contigsFile/02.renamedContigsFile"
removeL000_output_dir="/mnt/p/spf/projectMeta/06.contigsFile/03.L1000"

# 创建输出目录，如果不存在的话
mkdir -p "$renamed_output_dir"
mkdir -p "$removeL000_output_dir"


cd "$input_dir"

for fasta in *.fasta; do
    prefix="${fasta%.metaSPAdes.contigs.fasta}"
    awk -v prefix="$prefix" '{ if ($1 ~ /^>/) { print ">" prefix ".metaSPAdes." substr($0, 2) } else { print } }' "$fasta" > "${renamed_output_dir}/${prefix}.metaSPAdes.contigs.renamed.fasta" &

done

# 等待上述脚本执行完成
wait

conda activate seqkit

for fasta_file in "$renamed_output_dir"/*.renamed.fasta; do
    base_name=$(basename "$fasta_file" .renamed.fasta)
    output_file="${removeL000_output_dir}/${base_name}.L1000.fasta"
    seqkit seq -m 1000 "$fasta_file" > "$output_file"
done

conda deactivate

# 3. bowtie2 Indexing and Mapping

#!/bin/bash
conda_BASE=$(conda info --base)
source $conda_BASE/etc/profile.d/conda.sh

conda activate bowtie2

# 注意修改以下路径
input_dir="/mnt/p/spf/projectMeta/06.contigsFile/03.L1000"
bowtie2Index_output_dir="/mnt/p/spf/projectMeta/06.contigsFile/04.bowtie2Index"
fq_dir="/mnt/p/spf/projectMeta/04.hostFreeData"
bowtie2mapping_output_dir="/mnt/p/spf/projectMeta/06.contigsFile/05.bowtie2alignedData"


mkdir -p "$bowtie2Index_output_dir"
mkdir -p "${bowtie2mapping_output_dir}"

for fasta_file in "$input_dir"/*.L1000.fasta; do
    base_name=$(basename "$fasta_file" .L1000.fasta)
    index_prefix="${bowtie2Index_output_dir}/${base_name}.bowtie2.index"
    bowtie2-build -f "$fasta_file" "$index_prefix" --threads 2 &
done

# 等待上述的bowtie2全部建立完毕
 wait



for index in "${bowtie2Index_output_dir}"/*.bowtie2.index.1.bt2; do
    sample_id=$(basename "${index}" | sed 's/.metaSPAdes.contigs.bowtie2.index.1.bt2//')
    fq1="${fq_dir}/${sample_id}.hostFree.R1.fq.gz"
    fq2="${fq_dir}/${sample_id}.hostFree.R2.fq.gz"
    output="${bowtie2mapping_output_dir}/${sample_id}.metaSPAdes.contigs.bowtie2.sam"

    bowtie2 \
        -x "${bowtie2Index_output_dir}/${sample_id}.metaSPAdes.contigs.bowtie2.index" \
        -1 "${fq1}" \
        -2 "${fq2}" \
        -S "${output}" \
        -p 2 \
        --very-sensitive-local &
done
wait
done

# 4. sam to bam and sorting

#!/bin/sh
conda_BASE=$(conda info --base)
source $conda_BASE/etc/profile.d/conda.sh

conda activate samtools

input_dir="/mnt/p/spf/projectMeta/06.contigsFile/05.bowtie2alignedData"
output_dir="/mnt/p/spf/projectMeta/07.binning/bowtie2sorted"

mkdir -p "${output_dir}"

for sam_file in "${input_dir}"/*.sam; do
    sample_id=$(basename "${sam_file}" .metaSPAdes.contigs.bowtie2.sam)
    sorted_bam_file="${output_dir}/${sample_id}.metaSPAdes.contigs.bowtie2.sorted.bam"

    # Convert SAM to BAM, sort BAM file, and index sorted BAM file using a pipeline
    samtools view -@ 4 -bS "${sam_file}" | samtools sort -@ 4 -o --write-index -o "${sorted_bam_file}" &
done
wait
rm -rf ${input_dir}/*.metaSPAdes.contigs.bowtie2.sam

# 5. Call Depth

#!/bin/sh
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

input_dir="/mnt/p/spf/projectMeta/07.binning/bowtie2sorted"
output_dir="/mnt/p/spf/projectMeta/07.binning/metabat2_depth"

# 创建输出目录
mkdir -p "$output_dir"

# 使用conda激活metabat2环境
conda activate metaBat2

# 遍历所有BAM文件并行计算contig深度
find "$input_dir" -type f -name "*.bowtie2.sorted.bam" | while read -r bam_file; do
  base_name=$(basename "$bam_file" .bowtie2.sorted.bam)
  depth_file="${output_dir}/${base_name}.depth.txt"
  echo "Processing: $base_name"
  jgi_summarize_bam_contig_depths --outputDepth "$depth_file" "$bam_file" &
done
wait

conda deactivate

# 6. metaBat2 Binning

#!/bin/bash
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate metaBat2

fasta_dir="/mnt/p/spf/projectMeta/06.contigsFile/03.L1000"
depth_dir="/mnt/p/spf/projectMeta/07.binning/metabat2_depth"
output_dir="/mnt/p/spf/projectMeta/07.binning/binningResult"

# 创建输出目录，如果不存在
mkdir -p "$output_dir"

# 对每个.fasta文件和对应的_depth.txt文件执行Metabat2
for fasta_file in "$fasta_dir"/*.metaSPAdes.contigs.L1000.fasta; do
    base_name=$(basename "$fasta_file" .metaSPAdes.contigs.L1000.fasta)
    depth_file="$depth_dir/${base_name}.metaSPAdes.contigs.depth.txt"

    # 检查_depth.txt文件是否存在
    if [ -f "$depth_file" ]; then
        output_prefix="$output_dir/${base_name}"
        metabat2 -i "$fasta_file" -a "$depth_file" -o "${output_prefix}_bin" -m 1500 -t 2 &
    else
        echo "Error: Depth file not found for $fasta_file"
    fi
done

conda deactivate

# 7. Quality control of MAGs

#!/bin/bash
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate checkm

bins_dir="/mnt/p/spf/projectMeta/07.binning/binningResult"
output_dir="/mnt/p/spf/projectMeta/07.binning/checkmResult"

# 创建输出目录，如果不存在
mkdir -p "$output_dir"

# 检查当前目录中的bins文件数量
bins_count=$(find "$bins_dir" -name "*.fa" | wc -l)

# 计算需要的批次数量
num_batches=$(($bins_count / 1000 + 1))

# 对每个批次的 bins 文件执行 CheckM
for i in $(seq 1 $num_batches); do
    batch_bins_dir="$bins_dir/batch_$i"
    mkdir -p "$batch_bins_dir"
    find "$bins_dir" -maxdepth 1 -name "*.fa" -print | head -n 1000 | xargs -I '{}' mv '{}' "$batch_bins_dir"
    output_batch_dir="$output_dir/batch_$i"
    mkdir -p "$output_batch_dir"
    checkm lineage_wf -t 128 -x fa --pplacer_threads 32 "$batch_bins_dir" "$output_batch_dir"
done