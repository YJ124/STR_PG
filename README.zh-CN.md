# STR-PG

STR-PG 使用固定 pointer node 表示每个 STR 位点，并把等位基因序列、重复数和
种群频率保存到外部 registry，从而将 STR 等位多样性与图拓扑解耦。

本版本的生产默认 likelihood backend 是 Smith–Waterman（SW）。修复后的
motif-aware pHMM 仍被保留，但仅作为可选、实验性和可复现比较后端。两个后端
读取完全相同的、与 backend 无关的 informative-read cache。

## 安装

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[test]"
```

## 快速运行

默认 SW：

```powershell
PowerShell -ExecutionPolicy Bypass -File scripts\run_demo.ps1
```

显式 pHMM：

```powershell
PowerShell -ExecutionPolicy Bypass -File scripts\run_phmm_demo.ps1
```

CLI 中不指定 backend 等价于：

```text
strpg genotype ... --likelihood-backend sw
```

pHMM 必须显式启用：

```text
strpg genotype ... --likelihood-backend phmm
```

## 分型流程

```text
FASTQ/BAM + GAF
  -> backend-independent read selection
  -> immutable informative reads
  -> SW or pHMM likelihood matrix
  -> diploid genotype likelihood
  -> population/HWE/BN/length prior
  -> posterior, GT and GQ
  -> optional mixture diagnostics
```

read selection 不允许调用 SW、pHMM、prior、candidate ranking、posterior 或
truth。cache 保存 read ID、mate、序列、质量、MAPQ、选择原因和输入文件
SHA256。

## 科学边界

- SW 是当前默认生产后端。
- pHMM 是实验后端；现有证据不支持其优于 SW。
- pointer topology 是 representation 设计，不声明其自身提高分型准确率。
- 当前 targeted-locus 结果不能表述为 genome-wide validation。
- 源项目没有许可证；公开发布前必须由版权持有人决定许可证，详见
  `LICENSE_REQUIRED.md`。

完整英文说明、输入输出、架构、复现和迁移说明请阅读 `README.md` 与 `docs/`。
