import json
import os
import gc
import ladderpath as lp
import time
from multiprocessing import Process
import psutil
import pandas as pd

# 获取当前进程
process = psutil.Process(os.getpid())

# 定义合法碱基集合
valid_bases = {"A", "T", "C", "G"}

def log_memory_usage(stage):
    """
    打印当前内存使用情况，方便检查内存增长点。
    """
    mem_info = process.memory_info()
    print(f"[{stage}] 当前内存使用：{mem_info.rss / 1024 ** 2:.2f} MB")

def process_batch(batch_index, current_sequences, output_folder):
    """
    处理一个批次的数据。
    """
    time.sleep(2)
    print(
        f"处理批次 {batch_index}，总长度：{sum(len(seq) for seq in current_sequences)}, 列表长度：{len(current_sequences)}"
    )
    # 调用 ladderpath 函数
    lpjson = lp.get_ladderpath(
        current_sequences,
        info='V1.0.0.20240910_Alpha',
        estimate_eta=False,
        estimate_eta_para=[10, 'global'],
        save_file_name=f'./ladder_{batch_index}.json',
        show_version=True
    )

    # 检查返回值是否有效
    if lpjson is None:
        print(f"批次 {batch_index} 的输入无效，跳过处理。")
    else:
        # 显示梯径的 3 个指标
        index3 = lp.disp3index(lpjson)
        print(f"批次 {batch_index} 的梯径指标：{index3}")

        # 得到梯径的偏序多重集（POM）表示
        pom, pom_str = lp.POM_from_JSON(lpjson, display_str=False)

        pom_file_path = f"{output_folder}/pom_data_batch_{batch_index}.json"
        with open(pom_file_path, 'w') as json_file:
            json.dump(pom, json_file, indent=4)
        print(f"批次 {batch_index} 的 POM 字典已保存到 {pom_file_path}")

    # 在子进程处理完任务后进行垃圾回收
    print(f"批次 {batch_index} 处理完，开始垃圾回收")
    gc.collect()  # 手动触发垃圾回收
    log_memory_usage(f"处理批次 {batch_index} 后")  # 内存检查点

# 定义生成器函数
def sequence_generator(sequences, max_length):
    current_length = 0
    current_sequences = []
    for sequence in sequences:
        sequence = sequence.upper()  # 转为大写
        # 分割序列：在非标准字符处断开
        fragments = []
        fragment = []
        for char in sequence:
            if char in valid_bases:
                fragment.append(char)  # 构建有效片段
            else:
                if fragment:
                    fragments.append("".join(fragment))  # 保存有效片段
                    fragment = []  # 重置片段
        if fragment:
            fragments.append("".join(fragment))  # 保存最后的有效片段

        # 将分割后的片段加入 current_sequences
        for frag in fragments:
            frag_length = len(frag)

            # 如果片段长度小于等于1，则跳过
            if frag_length <= 1:
                continue

            # 如果片段长度超过 max_length，则切割该片段
            if frag_length > max_length:
                split_fragments = [frag[i:i + max_length] for i in range(0, len(frag), max_length)]
                for split_frag in split_fragments:
                    if current_length + len(split_frag) > max_length:
                        yield current_sequences
                        current_sequences = []
                        current_length = 0
                    current_sequences.append(split_frag)
                    current_length += len(split_frag)
            else:
                if current_length + frag_length > max_length:
                    yield current_sequences
                    current_sequences = []
                    current_length = 0
                current_sequences.append(frag)
                current_length += frag_length

    # 生成最后一个批次（如果有剩余）
    if current_sequences:
        yield current_sequences

# 主处理逻辑
if __name__ == "__main__":
    max_length = 1000000  # 最大长度限制
    output_folder = r'./hg19_json_file'
    file_path = r"./hg19_gene_annotation.csv"
    data = pd.read_csv(file_path)
    sequences = data['DNA_Sequence']

    num_workers = 50 # 控制并发数，根据 CPU 和内存调整
    processes = []  # 存放进程的列表
    batch_index = 1
    # 使用生成器创建参数
    for current_sequences in sequence_generator(sequences, max_length):
        # 创建并启动每个子进程
        process = Process(target=process_batch, args=(batch_index, current_sequences, output_folder))
        process.start()
        processes.append(process)

        # 控制并发数，等待部分进程结束
        if len(processes) >= num_workers:
            for p in processes:
                p.join()  # 等待子进程结束
            processes = []  # 清空进程列表
        batch_index += 1

    # 等待所有进程结束
    for p in processes:
        p.join()

    print("主程序执行完毕")

