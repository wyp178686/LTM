from Bio import SeqIO
import pandas as pd


def parse_fasta(fasta_file):
    """ 解析FASTA文件，返回一个字典 {染色体: 序列} """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = record.seq
    return sequences


def parse_csv(csv_file):
    """ 解析CSV文件，返回一个字典，键是基因ID，值是基因的位置信息 """
    genes = {}
    df = pd.read_csv(csv_file)

    for _, row in df.iterrows():
        gene_id = row['gene_id']
        seqid = row['seqid']
        start = row['start']
        end = row['end']
        strand = row['strand']
        #if strand == '+':
        # 只处理具有gene_id的记录，并且只提取相同gene_id的第一行
        genes[gene_id] = {
                "seqid": seqid,
                "start": start,
                "end": end,
                "strand": strand
            }

    return genes

def extract_gene_sequences(fasta_sequences, gff_genes):
    """ 根据基因位置信息从FASTA文件中提取基因的DNA序列，处理重叠基因 """
    gene_sequences = {}
    chromosome_boundaries = {}

    # 将字典转换为列表以便索引
    gene_list = list(gff_genes.items())

    for i, (gene_id, gene_info) in enumerate(gene_list):
        seqid = gene_info["seqid"]
        start = gene_info["start"]
        end = gene_info["end"]
        strand = gene_info["strand"]

        if seqid not in fasta_sequences:
            continue  # 跳过没有该染色体的基因

        chromosome_seq = fasta_sequences[seqid]

        # 检查是否有下一个基因
        if i < len(gene_list) - 1:
            next_gene_id, next_gene_info = gene_list[i + 1]
            next_seqid = next_gene_info['seqid']

            # 检查是否在同一染色体上
            if seqid == next_seqid:
                # 检查条件
                if start < next_gene_info['end']:
                    # 调整起始位置
                    start = next_gene_info['end'] + 1

        # 确保起始和结束位置在有效范围内
        if start < 1:
            start = 1
        if end > len(chromosome_seq):
            end = len(chromosome_seq)

        # 提取基因序列
        gene_seq = chromosome_seq[start - 1:end]  # GFF文件是1基准，Python是0基准

        gene_sequences[f"{seqid}_{gene_id}"] = {
            "sequence": str(gene_seq),
            "strand": strand,
            "start": start,
            "end": end
        }

        # 记录染色体的起始和结束位置
        if seqid not in chromosome_boundaries:
            chromosome_boundaries[seqid] = {"start": start, "end": end}
        else:
            chromosome_boundaries[seqid]["start"] = min(chromosome_boundaries[seqid]["start"], start)
            chromosome_boundaries[seqid]["end"] = max(chromosome_boundaries[seqid]["end"], end)

    # 提取染色体开头和结尾的非编码区域
    for seqid, boundaries in chromosome_boundaries.items():
        chromosome_seq = fasta_sequences[seqid]
        # 开头部分
        if boundaries["start"] > 1:
            start_noncoding_seq = chromosome_seq[:boundaries["start"] - 1]
            gene_sequences[f"{seqid}_start_noncoding"] = {
                "sequence": str(start_noncoding_seq),
                "strand": "+",
                "start": 1,
                "end": boundaries["start"] - 1
            }
        # 结尾部分
        if boundaries["end"] < len(chromosome_seq):
            end_noncoding_seq = chromosome_seq[boundaries["end"]:]
            gene_sequences[f"{seqid}_end_noncoding"] = {
                "sequence": str(end_noncoding_seq),
                "strand": "+",
                "start": boundaries["end"] + 1,
                "end": len(chromosome_seq)
            }

    return gene_sequences

def extract_non_coding_regions_dict(fasta_sequences, gff_genes):
    """ 提取非编码区序列并存储到字典 """
    non_coding_dict = {}
    data = []
    for gene_id, info in gff_genes.items():
        data.append({
            "gene_id": gene_id,
            "seqid": info["seqid"],
            "start": info["start"],
            "end": info["end"],
            "strand": info["strand"]
        })

        # 创建DataFrame并按染色体降序排序
    gene_df = pd.DataFrame(data)
    gene_df_sorted = gene_df.sort_values(['seqid', 'end'], ascending=[True, False])

    print("排序后的基因数据:")
    print(gene_df_sorted.head())
    print(f"总基因数: {len(gene_df_sorted)}")

    for seqid, group in gene_df_sorted.groupby("seqid"):  # 按染色体分组
        if seqid not in fasta_sequences:
            print(f"警告：染色体 {seqid} 不在参考序列中")
            continue

        chromosome_seq = fasta_sequences[seqid]
        previous_start = None  # 记录前一个基因的起始位置
        previous_gene_id = None  # 记录前一个基因的名称

        # 重置计数器
        non_coding_count = 0

        for _, row in group.iterrows():
            start, end = row["start"], row["end"]
            gene_id = row["gene_id"]
            strand = row["strand"]

            # 第一次迭代时跳过
            if previous_start is None:
                previous_start = start
                previous_gene_id = gene_id
                continue

                # 检查是否有非编码区
            if previous_start > end:
                non_coding_start = end + 1
                non_coding_end = previous_start - 1

                # 额外检查非编码区长度
                if non_coding_end >= non_coding_start:
                    try:
                        non_coding_seq = chromosome_seq[non_coding_start - 1:non_coding_end]

                        # 构造键名：染色体名 + 上一个基因名称 + 非编码区计数 + "_non_coding"
                        key = f"{seqid}_{previous_gene_id}_non_coding_{non_coding_count}"
                        non_coding_dict[key] = {
                            "sequence": str(non_coding_seq),
                            "strand": strand,
                            "start": non_coding_start,
                            "end": non_coding_end
                        }

                        # 打印非编码区信息
                        print(f"非编码区: {key}")
                        print(f"  长度: {len(non_coding_seq)}")
                        print(f"  起始: {non_coding_start}")
                        print(f"  结束: {non_coding_end}")

                        non_coding_count += 1
                    except IndexError:
                        print(f"警告：提取非编码区时出错 {seqid}: {non_coding_start} - {non_coding_end}")

                        # 更新前一个基因的起始位置和基因名称
            previous_start = start
            previous_gene_id = gene_id

    print(f"总非编码区数: {len(non_coding_dict)}")
    return non_coding_dict



def save_to_csv(gene_sequences, output_file):
    """ 将基因ID和序列保存到CSV文件 """
    data = []
    for gene_id, info in gene_sequences.items():
        data.append({
            "Gene_ID": gene_id,
            "DNA_Sequence": info["sequence"],
            "Strand": info["strand"],
            "Start": info["start"],
            "End": info["end"]
        })
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)


# 主程序
def main():
    fasta_file = r"./hg19.fa"
    csv_file = r"./hg19_zhengfu_sorted_output_file_end.csv"
    output_file = r"./hg19_gene_annotation.csv"

    # 解析FASTA和CSV文件
    fasta_sequences = parse_fasta(fasta_file)
    gff_genes = parse_csv(csv_file)
    print(gff_genes)

    # 提取基因序列和非编码区域
    gene_sequences = extract_gene_sequences(fasta_sequences, gff_genes)
    non_coding_regions = extract_non_coding_regions_dict(fasta_sequences, gff_genes)

    # 将结果保存为CSV文件
    all_sequences = {**gene_sequences, **non_coding_regions}
    save_to_csv(all_sequences, output_file)


if __name__ == "__main__":
    main()