import json
import os
from collections import defaultdict

# 文件路径定义
file_paths = [
    './hg19_json_file_1000000',
]

# 收集所有 JSON 文件路径
files = []
for file_path in file_paths:
    files.extend([os.path.join(file_path, f) for f in os.listdir(file_path) if f.endswith('.json')])

# 使用 defaultdict 来统计重数
token_counts = defaultdict(int)

# 读取文件并统计重数
for file in files:
    with open(file, 'r', encoding='utf-8') as f:
        data = json.load(f)
        for layer in data.values():  # 遍历每一层
            for token, count in layer.items():
                # 仅当键长度大于等于 2 时才进行统计
                if len(token) >= 1:
                    token_counts[token] += count  # 累加重数

# 打印字典信息
print(f'字典的长度为：{len(token_counts)}')

# 计算所有键的长度累积
total_length = sum(len(key) for key in token_counts)
print("字典中所有键的长度累积是：", total_length)

# 写出结果到 JSON 文件
output_file = './merge_multiplicities.json'  # 输出文件名
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(token_counts, f, ensure_ascii=False, indent=4)