# PanSnpEff: 泛基因组变异注释工具（支持 GFA 分支数统计）
<div align="center">
 <img width="400" height="400" alt="logo" src="https://github.com/user-attachments/assets/d6f00e15-4456-4924-9c7b-f61a108f90d8" />
</div>

A SnpEff-like variant annotation tool tailored for pan-genome analysis, supporting branch count statistics from GFA files and standard VCF output.
# 项目简介
PanSnpEff 是一款基于 Python 的变异注释工具，核心功能类似经典工具 SnpEff，但额外支持从泛基因组组装图（GFA 文件）中提取参考序列并统计基因区域的分支数，最终输出包含功能影响分级和分支结构信息的标准 VCF 文件，适用于泛基因组变异功能解析场景。
# 核心价值
兼容传统基因组（FASTA）和泛基因组（GFA）两种参考格式

提供与 SnpEff 一致的变异影响程度分级（HIGH/MODERATE/LOW/MODIFIER）

创新添加基因区域分支数（BRANCH_COUNT），关联变异与泛基因组结构多样性

输出符合 VCFv4.2 规范，可直接用于下游分析（如 bcftools、IGV）

# 核心功能
- 参考序列提取	优先从 GFA 的 P/S 行提取泛基因组参考序列，FASTA 作为备选补充
- GFF 基因结构解析	解析 gene/exon/CDS/transcript 等区域，构建染色体 - 基因 - 区域的层级索引
- 变异功能注释	预测变异所在区域（基因间区 / 内含子 / 外显子 / CDS）及对蛋白质的影响（错义 / 同义 / 移码等）
- SnpEff 兼容分级	基于变异类型自动映射影响程度（HIGH/MODERATE/LOW/MODIFIER）
- 泛基因组分支数统计	基于 GFA 的 L 行（片段连接），统计变异所在基因区域的分支数
- 标准 VCF 输出	在原始 VCF 中添加 7 个注释字段，保留所有原始信息
# 安装依赖
环境要求:
Python 3.7+

依赖库：biopython（用于序列处理和密码子翻译）
# 推荐使用conda创建独立环境
conda create -n pansnpeff python=3.9 -y

conda activate pansnpeff

# 安装依赖库
pip install biopython
# 使用方法
输入文件

需提供 3 个必需文件 + 1 个可选文件：
1.VCF 文件：待注释的变异文件（VCFv4.2+）
2.GFF 文件：基因结构注释文件（需包含 gene/exon/CDS 记录）
3.GFA 文件：泛基因组组装图（需包含 P/S/L 行，用于提取参考和统计分支）
4.FASTA 文件（可选）：传统参考基因组，当 GFA 参考提取失败时使用（填 none 表示不使用）
# 基本用法（使用GFA作为参考，不使用FASTA）
`python main.py input.vcf genomic.gff pangenome.gfa none`
# 备选用法（GFA优先，FASTA作为补充）
`python main.py input.vcf genomic.gff pangenome.gfa ref_genome.fasta`
# 输出文件
输出文件名为 annotated_with_impact_branches_<输入VCF名>.vcf，在原始 VCF 基础上新增以下 INFO 字段：
|字段名	|类型	|描述|
|------|-----|-----|
|GENE_ID	|String	|变异所在的基因 ID（. 表示基因间区）|
|REGION	|String	|变异所在区域（intergenic/intron/exon/CDS）|
|EFFECT	|String	|变异的生物学影响|
|CODON_CHANGE	|String	|CDS 区域 SNV 的密码子变化（非 CDS 区域为 .）|
|AA_CHANGE	|String	|CDS 区域 SNV 的氨基酸变化（非 CDS 区域为 .）|
|BRANCH_COUNT	|Integer	|变异所在基因区域的分支数（反映泛基因组结构多样性）|
|IMPACT	|String	|变异影响程度（HIGH/MODERATE/LOW/MODIFIER）|
# 注意事项
文件格式兼容性：
GFF 文件需为 GFF3 格式，染色体 ID（seqid）需与 GFA/FASTA 的参考名完全一致

GFA 文件需包含完整的 P 行（路径，定义参考序列）、S 行（片段序列）和 L 行（片段连接，用于分支统计）

VCF 文件需为 VCFv4.2 及以上版本，CHROM 字段需与参考序列名匹配

## 分支数统计逻辑：
分支数定义为 GFA 中与该基因区域重叠的片段连接数（L 行数量）

若 GFA 中无 L 行或片段无位置信息（PO 标签），分支数将显示为 0
