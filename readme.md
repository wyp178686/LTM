好的，我已经根据您的要求，将 `Merge_the_multiplicities.py` 脚本作为一个独立的步骤整合到您现有的 `README.md` 文件中。

新的工作流程现在包含五个步骤，逻辑更加清晰。

-----

# Genomic Language Model Pre-training Pipeline with Ladderpath Tokenization

## 项目概述

本项目提供了一套完整的端到端工作流程，用于预训练一个基于 **Ladderpath** 分词算法的基因组语言模型。整个流程从原始的基因组序列（FASTA 格式）和基因注释文件开始，通过数据预处理、并行化 Ladderpath 程序计算、结果合并、词汇表构建，最终使用 Hugging Face Transformers 库训练一个 BERT 架构的掩码语言模型。

## 工作流程

整个流程分为五个主要步骤，每个步骤由一个独立的 Python 脚本完成。请务必按照以下顺序执行：

1.  **`split_chrom.py`**: **数据预处理**

      * 读取一个 FASTA 格式的参考基因组文件和一个 CSV 格式的基因注释文件。
      * 提取所有基因序列和基因间的非编码区序列。
      * 处理重叠基因，并提取染色体两端的序列。
      * 将所有提取出的 DNA 序列整合成一个统一的 CSV 文件，作为下一步的输入。

2.  **`Ladderpath_multiprocess_DNASequence.py`**: **Ladderpath 并行计算**

      * 读取上一步生成的包含所有 DNA 序列的 CSV 文件。
      * 将长序列分割成适合处理的批次（默认总字符数不超过 100 万）。
      * 利用多进程并行运行 Ladderpath 算法，对每个批次的序列进行层级化结构分析。
      * 将每个批次分析得到的偏序多重集 (Poset-Multiset, POM) 结果保存为独立的 JSON 文件。

3.  **`Merge_the_multiplicities.py`**: **合并并行结果**

      * 读取第二步并行计算生成的所有 JSON 文件。
      * 将所有文件中每个“梯元”(Token) 的出现次数（重数）进行累加。
      * 将最终聚合的结果保存到一个单独的、包含全局统计信息的 JSON 文件中 (`merge_multiplicities.json`)。

4.  **`build_vocab.py`**: **统计与词汇表构建**

      * 读取上一步生成的合并后的 JSON 文件 (`merge_multiplicities.json`)。
      * 根据重数对所有梯元进行降序排序，并选择排名最高的 N 个梯元（例如 1000 个）构建最终的词汇表文件 `vocab.txt`。

5.  **`LPT_pretrain.py`**: **模型预训练**

      * 这是一个基于 Hugging Face Transformers 的标准模型训练脚本。
      * 它使用第四步生成的 `vocab.txt` 文件和预先分词好的训练数据来预训练一个 BERT 模型。
      * **注意**: 此脚本要求训练数据 (`train_data_file`) 是一个文本文件，其中每行是使用 `vocab.txt` 转换后的 Token ID，并以空格分隔。

## 环境准备

### 1\. 克隆项目

```bash
git clone <your-repo-url>
cd <your-repo-name>
```

### 2\. 安装依赖

建议使用 `conda` 创建一个独立的虚拟环境。

```bash
conda create -n genomics_lm python=3.9
conda activate genomics_lm
```

然后通过 `pip` 安装所有必要的库。

```bash
pip install pandas = 2.2.3
biopython = 1.84
psutil = 6.1.0
torch = 2.6.0
transformers = 4.47.1
tqdm = 4.67.1
```

## 使用方法

### 步骤 1: 执行数据预处理

首先，确保你的参考基因组文件 (如 `hg19.fa`) 和基因注释文件 (如 `hg19_zhengfu_sorted_output_file_end.csv`) 已准备好。然后修改 `split_chrom.py` 脚本中的文件路径：

```python
# split_chrom.py
...
def main():
    fasta_file = r"./hg19.fa"  # <-- 修改为你的 FASTA 文件路径
    csv_file = r"./hg19_zhengfu_sorted_output_file_end.csv"  # <-- 修改为你的注释文件路径
    output_file = r"./hg19_gene_annotation.csv" # <-- 输出文件路径
...
```

然后运行脚本：

```bash
python split_chrom.py
```

执行成功后，将生成 `hg19_gene_annotation.csv` 文件。

### 步骤 2: 执行 Ladderpath 分析

修改 `Ladderpath_multiprocess_DNASequence.py` 脚本中的输入文件路径、输出文件夹路径以及并行工作进程数。

```python
# Ladderpath_multiprocess_DNASequence.py
...
if __name__ == "__main__":
    max_length = 1000000  # 每个批次的最大总字符数
    output_folder = r'./hg19_json_file'  # <-- POM JSON 文件的输出目录
    file_path = r"./hg19_gene_annotation.csv"  # <-- 使用上一步生成的 CSV 文件
    
    num_workers = 50 # <-- 根据你的 CPU 核心数和内存进行调整
...
```

运行脚本：

```bash
python Ladderpath_multiprocess_DNASequence.py
```

该脚本会耗时较长，执行完毕后，`hg19_json_file` 目录下会生成大量的 JSON 文件。

### 步骤 3: 合并并行结果

修改 `Merge_the_multiplicities.py` 脚本中的输入目录和输出文件路径。

```python
# Merge_the_multiplicities.py
...
# 文件路径定义
file_paths = [
    './hg19_json_file', # <-- 确保这是上一步的输出目录
]
...
# 写出结果到 JSON 文件
output_file = './merge_multiplicities.json'  # <-- 输出文件名
...
```

运行脚本：

```bash
python Merge_the_multiplicities.py
```

执行成功后，将生成 `merge_multiplicities.json` 文件，其中包含了所有 Token 的全局频率。

### 步骤 4: 构建词汇表

修改 `build_vocab.py` 脚本，使其**直接读取上一步生成的 `merge_multiplicities.json` 文件**，而不是重新扫描整个目录。

```python
# build_vocab.py (建议修改)
import json

# 直接读取合并后的文件
with open('./merge_multiplicities.json', 'r') as f:
    data_dict = json.load(f)

# 按频率降序排序
sorted_items = sorted(data_dict.items(), key=lambda item: item[1], reverse=True)

# 提取 Top K 词元 (可修改 K 的值)
top_k_keys = [key for key, value in sorted_items[:1000]]

# 输出词汇表
with open('./vocab.txt', 'w') as f:
    for key in top_k_keys:
        f.write(key + '\n')

print("Finished writing top keys to vocab.txt")
```

运行修改后的脚本：

```bash
python build_vocab.py
```

执行成功后，将生成最终的 `vocab.txt` 文件。

### 步骤 5: 训练模型

这是最后一步。你需要准备好以下文件：

1.  **`vocab.txt`**: 上一步生成的词汇表。
2.  **`bert_config.json`**: 模型的配置文件（例如，层数、隐藏层大小等）。
3.  **分词后的训练数据**: 一个 `.txt` 文件，其中每一行都是由空格分隔的 Token ID。你需要自己写脚本，使用 `vocab.txt` 将 `hg19_gene_annotation.csv` 中的序列转换为这个格式。

然后通过命令行运行训练脚本。以下是一个示例命令：

```bash
python LPT_pretrain.py \
    --output_dir ./lpt_model_output \
    --model_type bert \
    --do_train \
    --train_data_file ./path/to/your/tokenized_data.txt \
    --tokenizer_name ./ \
    --config_name ./bert_config.json \
    --mlm \
    --mlm_probability 0.15 \
    --per_gpu_train_batch_size 16 \
    --gradient_accumulation_steps 4 \
    --num_train_epochs 3 \
    --learning_rate 5e-5 \
    --save_steps 10000 \
    --logging_steps 500
```

训练完成后，模型检查点和相关文件将保存在 `--output_dir` 指定的目录中。

## 文件结构（示例）

```
.
├── split_chrom.py
├── Ladderpath_multiprocess_DNASequence.py
├── Merge_the_multiplicities.py
├── build_vocab.py
├── LPT_pretrain.py
├── README.md
│
├── data/
│   ├── hg19.fa
│   └── hg19_zhengfu_sorted_output_file_end.csv
│
├── hg19_gene_annotation.csv        # 步骤 1 输出
│
├── hg19_json_file/                   # 步骤 2 输出
│   ├── pom_data_batch_1.json
│   ├── pom_data_batch_2.json
│   └── ...
│
├── merge_multiplicities.json       # 步骤 3 输出
├── vocab.txt                       # 步骤 4 输出
│
├── tokenized_data.txt              # 步骤 5 输入 
├── bert_config.json                # 步骤 5 输入
│
└── lpt_model_output/                 # 步骤 5 输出
    ├── checkpoint-10000/
    ├── checkpoint-20000/
    └── ...
```