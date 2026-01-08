import matplotlib.pyplot as plt
import collections
import os

def verify_heterogeneity():
    fastq_file = "sim_data_R1.fastq"
    if not os.path.exists(fastq_file): return

    print(f"正在分析 {fastq_file}...")
    counts = collections.defaultdict(int)
    
    with open(fastq_file, 'r') as f:
        while True:
            line = f.readline()
            if not line: break
            seq = f.readline().strip()
            f.readline(); f.readline() 
            
            # 滑动窗口查找最长 CAG
            best_run = 0
            for frame in range(3):
                current_run = 0
                max_run_in_frame = 0
                for i in range(frame, len(seq)-2, 3):
                    if seq[i:i+3] == "CAG":
                        current_run += 1
                    else:
                        max_run_in_frame = max(max_run_in_frame, current_run)
                        current_run = 0
                max_run_in_frame = max(max_run_in_frame, current_run)
                best_run = max(best_run, max_run_in_frame)
            
            if best_run >= 10: 
                counts[best_run] += 1

    x = sorted(counts.keys())
    y = [counts[k] for k in x]
    
    plt.figure(figsize=(12, 6))
    bars = plt.bar(x, y, color='#A0CBE8', edgecolor='gray', alpha=0.9)
    
    # 【预期目标】
    # 15, 18: 正常二倍体，应该很高
    # 40: 扩增目标。由于 Read 150bp, 40repeats=120bp, 只有 30bp 空间留给侧翼。
    # 所以大部分 Reads 会在 30-38 之间饱和，只有极少数完美覆盖的 Reads 能达到 40。
    expected_peaks = [15, 18, 40]
    
    for rect, val in zip(bars, y):
        height = rect.get_height()
        repeat = int(rect.get_x() + rect.get_width()/2.0)
        
        # 标记数值
        if height > max(y)*0.1:
            plt.text(rect.get_x() + rect.get_width()/2.0, height, f'{repeat}', 
                     ha='center', va='bottom', fontsize=8)
            
        # 高亮逻辑
        # 15, 18 用橙色
        if abs(repeat - 15) <= 1 or abs(repeat - 18) <= 1:
            rect.set_color('#F28E2B')
        # 30-40 都是扩增信号区间 (红色)
        # 因为 150bp Read 测 40 repeats 很容易被截断在 30多
        elif 30 <= repeat <= 40:
            rect.set_color('#E15759')

    plt.title("Allele Distribution: 15/18 (Normal) + 40 (Tumor Expansion)")
    plt.xlabel("CAG Repeats visible in Read")
    plt.ylabel("Read Count")
    
    # 说明文字
    note = "Note:\nPeak at 40 is rare (requires perfect centering).\nSignals in 30-40 range indicate expansion > 30."
    plt.text(0.95, 0.95, note, transform=plt.gca().transAxes, 
             va='top', ha='right', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.savefig("output_verification.png")
    print(f"分析完成。请查看 output_verification.png")
    print(f"预期现象：15/18 双峰耸立；30-40 区间有连续的信号分布（代表扩增等位基因 40）。")

if __name__ == "__main__":
    verify_heterogeneity()