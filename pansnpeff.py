import sys
import re
from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO

# --------------------------
# 变异影响程度映射表
# --------------------------
IMPACT_MAP = {
    # 高影响：直接破坏蛋白质功能
    'nonsense_mutation': 'HIGH',
    'frameshift_variant': 'HIGH',
    'nonstop_mutation': 'HIGH',
    'start_lost': 'HIGH',
    'stop_lost': 'HIGH',

    # 中影响：可能改变蛋白质功能
    'missense_mutation': 'MODERATE',
    'inframe_insertion': 'MODERATE',
    'inframe_deletion': 'MODERATE',
    'splice_acceptor_variant': 'MODERATE',
    'splice_donor_variant': 'MODERATE',

    # 低影响：对蛋白质功能影响较小
    'synonymous_variant': 'LOW',
    'coding_sequence_variant': 'LOW',

    # 修饰性：不直接影响蛋白质序列
    'intergenic_variant': 'MODIFIER',
    'intron_variant': 'MODIFIER',
    'exon_variant': 'MODIFIER',
    'genic_variant': 'MODIFIER',
    'reference_mismatch': 'MODIFIER',
    'incomplete_codon_variant': 'MODIFIER'
}


# --------------------------
# 1. 参考序列提取
# --------------------------
def extract_reference(gfa_file, fasta_file=None):
    """从GFA优先提取参考序列，GFA不足时使用FASTA补充"""
    ref_from_gfa, chrom_map = _extract_from_gfa(gfa_file)

    # 若GFA未提取到序列，尝试从FASTA提取
    if not ref_from_gfa and fasta_file:
        ref_from_gfa = _extract_from_fasta(fasta_file)
        chrom_map = {k: k for k in ref_from_gfa.keys()}

    if not ref_from_gfa:
        raise ValueError("无法从GFA或FASTA提取参考序列")

    return ref_from_gfa, chrom_map


def _extract_from_gfa(gfa_file):
    """从GFA的P行和S行提取参考序列"""
    seg_sequences = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if fields[0] == 'S' and len(fields) >= 3:
                seg_id = fields[1]
                seg_seq = fields[2].upper()
                seg_sequences[seg_id] = seg_seq

    ref_genome = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if fields[0] == 'P' and len(fields) >= 3:
                ref_name = fields[1]
                seg_list = fields[2].split(',')
                full_seq = []
                for seg in seg_list:
                    if seg.endswith(('+', '-')):
                        seg_id = seg[:-1]
                        orientation = seg[-1]
                    else:
                        seg_id = seg
                        orientation = '+'

                    if seg_id not in seg_sequences:
                        continue  # 跳过未找到的片段

                    seg_seq = seg_sequences[seg_id]
                    if orientation == '-':
                        seg_seq = str(Seq(seg_seq).reverse_complement())

                    full_seq.append(seg_seq)

                ref_genome[ref_name] = ''.join(full_seq)
                print(f"从GFA P行提取参考序列：{ref_name}（长度：{len(ref_genome[ref_name])}bp）")

    return ref_genome, {k: k for k in ref_genome.keys()}


def _extract_from_fasta(fasta_file):
    """从FASTA文件提取参考序列"""
    ref_genome = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        ref_genome[record.id] = str(record.seq).upper()
    return ref_genome


# --------------------------
# 2. GFF文件解析
# --------------------------
def parse_gff(gff_file, target_chroms):
    """解析GFF文件，提取基因、外显子、CDS等结构信息"""
    gff_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    transcript_info = {}

    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != 9:
                continue

            seqid, _, feature, start, end, _, strand, _, attr = fields
            if seqid not in target_chroms:
                continue  # 只保留目标染色体的信息
            start = int(start)
            end = int(end)

            attr_dict = _parse_attributes(attr)

            if feature == 'gene':
                _process_gene(gff_data, seqid, attr_dict, start, end, strand)
            elif feature == 'transcript':
                _process_transcript(transcript_info, attr_dict)
            elif feature == 'exon':
                _process_exon(gff_data, seqid, transcript_info, attr_dict, start, end, strand)
            elif feature == 'CDS':
                _process_cds(gff_data, seqid, transcript_info, attr_dict, start, end, strand)

    return gff_data, transcript_info


def _parse_attributes(attr_str):
    """解析GFF属性字段（如ID=gene123;Parent=transcript456）"""
    attr_dict = {}
    for item in attr_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attr_dict[k.strip()] = v.strip()
    return attr_dict


def _process_gene(gff_data, seqid, attr_dict, start, end, strand):
    gene_id = attr_dict.get('ID')
    if gene_id:
        gff_data[seqid][gene_id]['gene'].append((start, end, strand))


def _process_transcript(transcript_info, attr_dict):
    transcript_id = attr_dict.get('ID')
    gene_id = attr_dict.get('Parent')
    if transcript_id and gene_id:
        transcript_info[transcript_id] = {
            'gene_id': gene_id,
            'strand': '+',  # 后续会被CDS/exon的strand覆盖
            'cds_list': []
        }


def _process_exon(gff_data, seqid, transcript_info, attr_dict, start, end, strand):
    transcript_id = attr_dict.get('Parent')
    if not transcript_id or transcript_id not in transcript_info:
        return
    gene_id = transcript_info[transcript_id]['gene_id']
    gff_data[seqid][gene_id]['exon'].append((start, end, strand, transcript_id))
    transcript_info[transcript_id]['strand'] = strand  # 记录链方向


def _process_cds(gff_data, seqid, transcript_info, attr_dict, start, end, strand):
    transcript_id = attr_dict.get('Parent')
    if not transcript_id or transcript_id not in transcript_info:
        return
    gene_id = transcript_info[transcript_id]['gene_id']
    phase = attr_dict.get('Phase', '0')
    gff_data[seqid][gene_id]['CDS'].append((start, end, strand, transcript_id, phase))
    transcript_info[transcript_id]['cds_list'].append((start, end, phase))
    transcript_info[transcript_id]['cds_list'].sort(key=lambda x: x[0])  # 按位置排序
    transcript_info[transcript_id]['strand'] = strand  # 记录链方向


# --------------------------
# 3. 分支数统计
# --------------------------
def count_branches(gfa_file, gff_data, ref_names):
    """统计每个基因区域的分支数"""
    path_segments, seg_coords, links, seg_to_ref = _parse_gfa(gfa_file)
    region_segs = _map_regions_to_segments(gff_data, seg_coords, seg_to_ref, ref_names)
    return _count_branch_numbers(region_segs, links)


def _parse_gfa(gfa_file):
    """解析GFA获取片段、路径和连接信息"""
    path_segments = defaultdict(set)  # 参考序列→包含的片段
    seg_coords = {}  # 片段ID→(参考起始位置, 参考结束位置)
    links = []  # 片段连接关系 (from_seg, to_seg)
    seg_to_ref = {}  # 片段ID→所属参考序列

    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if fields[0] == 'P' and len(fields) >= 3:
                ref_name = fields[1]
                segments = fields[2].split(',')
                for seg in segments:
                    seg_id = seg[:-1] if seg.endswith(('+', '-')) else seg
                    path_segments[ref_name].add(seg_id)
                    seg_to_ref[seg_id] = ref_name
            elif fields[0] == 'S' and len(fields) >= 3:
                seg_id = fields[1]
                seg_seq = fields[2]
                # 提取PO标签（片段在参考序列中的位置）
                po_tag = next((t for t in fields[3:] if t.startswith('PO:f:')), None)
                if po_tag:
                    ref_start = int(float(po_tag.split(':')[2]))
                    ref_end = ref_start + len(seg_seq) - 1
                    seg_coords[seg_id] = (ref_start, ref_end)
            elif fields[0] == 'L' and len(fields) >= 5:
                from_seg = fields[1]
                to_seg = fields[3]
                links.append((from_seg, to_seg))

    return path_segments, seg_coords, links, seg_to_ref


def _map_regions_to_segments(gff_data, seg_coords, seg_to_ref, ref_names):
    """将基因区域（如exon）映射到GFA片段"""
    region_segs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    for ref_name in ref_names:
        if ref_name not in gff_data:
            continue
        genes = gff_data[ref_name]
        for gene_id, features in genes.items():
            for feature_type, regions in features.items():
                for region in regions:
                    # 只提取start和end（前两个值），兼容不同长度的区域记录
                    if len(region) < 2:
                        continue
                    reg_start, reg_end = region[0], region[1]

                    # 找到与该区域重叠的片段
                    for seg_id, (seg_start, seg_end) in seg_coords.items():
                        if seg_to_ref.get(seg_id) != ref_name:
                            continue
                        # 判断片段与区域是否重叠
                        if not (seg_end < reg_start or seg_start > reg_end):
                            region_segs[ref_name][gene_id][feature_type].add(seg_id)
    return region_segs


def _count_branch_numbers(region_segs, links):
    """统计每个区域的分支数（片段连接数）"""
    branch_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for ref_name, genes in region_segs.items():
        for gene_id, features in genes.items():
            for feature_type, segs in features.items():
                count = 0
                for (from_seg, to_seg) in links:
                    if from_seg in segs or to_seg in segs:
                        count += 1
                branch_counts[ref_name][gene_id][feature_type] = count
    return branch_counts


# --------------------------
# 4. 变异注释核心逻辑
# --------------------------
def annotate_variant(chrom, pos, ref, alt, gff_data, transcript_info, ref_genome, branch_counts):
    """注释单个变异，返回包含功能影响、分支数和影响程度的结果"""
    annotation = {
        'gene_id': '.',
        'region': 'intergenic',
        'effect': 'intergenic_variant',
        'codon_change': '.',
        'aa_change': '.',
        'branch_count': 0,
        'impact': 'MODIFIER'  # 变异影响程度
    }

    if chrom not in gff_data:
        annotation['impact'] = IMPACT_MAP.get(annotation['effect'], 'MODIFIER')
        return annotation

    # 遍历基因，判断变异是否在基因区域内
    for gene_id, features in gff_data[chrom].items():
        for (gene_start, gene_end, _) in features.get('gene', []):
            if gene_start <= pos <= gene_end:
                annotation['gene_id'] = gene_id
                annotation['region'] = 'genic'
                annotation = _check_regions(
                    chrom, pos, ref, alt, gene_id, features, transcript_info,
                    ref_genome, branch_counts, annotation
                )
                # 计算影响程度
                annotation['impact'] = IMPACT_MAP.get(annotation['effect'], 'MODIFIER')
                return annotation

    # 基因间区的影响程度
    annotation['impact'] = IMPACT_MAP.get(annotation['effect'], 'MODIFIER')
    return annotation


def _check_regions(chrom, pos, ref, alt, gene_id, features, transcript_info, ref_genome, branch_counts, annotation):
    """检查变异具体位于哪个基因区域（exon/intron/CDS）"""
    # 检查外显子
    for region in features.get('exon', []):
        if len(region) < 3:
            continue  # 确保有足够的信息
        exon_start, exon_end, strand, transcript_id = region[:4]
        if exon_start <= pos <= exon_end:
            annotation['region'] = 'exon'
            annotation['effect'] = 'exon_variant'
            annotation['branch_count'] = branch_counts.get(chrom, {}).get(gene_id, {}).get('exon', 0)
            return annotation

    # 检查CDS
    for region in features.get('CDS', []):
        if len(region) < 4:
            continue  # 确保有足够的信息
        cds_start, cds_end, strand, transcript_id = region[:4]
        if cds_start <= pos <= cds_end:
            annotation['region'] = 'CDS'
            annotation['branch_count'] = branch_counts.get(chrom, {}).get(gene_id, {}).get('CDS', 0)
            return _predict_cds_effect(
                chrom, pos, ref, alt, strand, transcript_id, transcript_info, ref_genome, annotation
            )

    # 内含子（非exon/CDS的基因内区域）
    annotation['region'] = 'intron'
    annotation['effect'] = 'intron_variant'
    annotation['branch_count'] = branch_counts.get(chrom, {}).get(gene_id, {}).get('intron', 0)
    return annotation


def _predict_cds_effect(chrom, pos, ref, alt, strand, transcript_id, transcript_info, ref_genome, annotation):
    """预测CDS区域变异对氨基酸的影响"""
    transcript = transcript_info.get(transcript_id)
    if not transcript:
        return annotation

    cds_list = transcript['cds_list']
    gene_strand = transcript['strand']
    full_cds = _get_full_cds_sequence(chrom, cds_list, ref_genome)
    variant_in_cds_pos = _get_variant_cds_position(pos, cds_list)

    if variant_in_cds_pos is None:
        return annotation

    # 处理负链（反向互补）
    if gene_strand == '-':
        variant_in_cds_pos = len(full_cds) - 1 - variant_in_cds_pos
        ref = str(Seq(ref).reverse_complement())
        alt = str(Seq(alt).reverse_complement())

    # 处理不同变异类型
    if len(ref) == 1 and len(alt) == 1:
        return _annotate_snv(full_cds, variant_in_cds_pos, ref, alt, annotation)
    elif len(ref) < len(alt):
        annotation['effect'] = 'frameshift_variant' if (len(alt) - len(ref)) % 3 != 0 else 'inframe_insertion'
    elif len(ref) > len(alt):
        annotation['effect'] = 'frameshift_variant' if (len(ref) - len(alt)) % 3 != 0 else 'inframe_deletion'

    return annotation


def _get_full_cds_sequence(chrom, cds_list, ref_genome):
    """拼接完整的CDS序列"""
    cds_sequence = []
    for (start, end, _) in cds_list:
        cds_seq = ref_genome[chrom][start - 1:end]  # GFF是1-based，序列是0-based
        cds_sequence.append(cds_seq)
    return ''.join(cds_sequence)


def _get_variant_cds_position(pos, cds_list):
    """计算变异在CDS序列中的位置（0-based）"""
    cds_offset = 0
    for (start, end, _) in cds_list:
        if start <= pos <= end:
            return cds_offset + (pos - start)
        cds_offset += (end - start + 1)
    return None


def _annotate_snv(full_cds, pos, ref, alt, annotation):
    """注释单核苷酸变异（SNV）"""
    if full_cds[pos] != ref:
        annotation['effect'] = 'reference_mismatch'
        return annotation

    codon_index = pos // 3
    codon_pos = pos % 3
    codon_start = codon_index * 3
    codon_end = codon_start + 3
    original_codon = full_cds[codon_start:codon_end]

    if len(original_codon) < 3:
        annotation['effect'] = 'incomplete_codon_variant'
        return annotation

    # 计算变异后的密码子和氨基酸
    mutated_codon = list(original_codon)
    mutated_codon[codon_pos] = alt
    mutated_codon = ''.join(mutated_codon)

    original_aa = str(Seq(original_codon).translate())
    mutated_aa = str(Seq(mutated_codon).translate())

    # 确定变异类型
    if mutated_aa == '*':
        annotation['effect'] = 'nonsense_mutation'
    elif original_aa == '*':
        annotation['effect'] = 'nonstop_mutation'
    elif original_aa != mutated_aa:
        annotation['effect'] = 'missense_mutation'
    else:
        annotation['effect'] = 'synonymous_variant'

    annotation['codon_change'] = f"{original_codon}→{mutated_codon}"
    annotation['aa_change'] = f"{original_aa}{codon_index + 1}{mutated_aa}"
    return annotation


# --------------------------
# 5. VCF文件处理
# --------------------------
def process_vcf(input_vcf, output_vcf, gff_data, transcript_info, ref_genome, branch_counts):
    """处理VCF文件，添加注释并输出"""
    with open(input_vcf, 'r') as in_f, open(output_vcf, 'w') as out_f:
        # 处理头部
        _write_vcf_header(in_f, out_f)

        # 处理数据行
        _process_vcf_data(in_f, out_f, gff_data, transcript_info, ref_genome, branch_counts)


def _write_vcf_header(in_f, out_f):
    """写入VCF头部，添加新注释字段定义"""
    for line in in_f:
        line = line.strip()
        if line.startswith('##'):
            out_f.write(f"{line}\n")
        elif line.startswith('#CHROM'):
            # 添加新增INFO字段的定义
            out_f.write('##INFO=<ID=GENE_ID,Number=1,Type=String,Description="变异所在的基因ID">\n')
            out_f.write(
                '##INFO=<ID=REGION,Number=1,Type=String,Description="变异所在的基因区域（intergenic/intron/exon/CDS）">\n')
            out_f.write('##INFO=<ID=EFFECT,Number=1,Type=String,Description="变异的生物学影响">\n')
            out_f.write('##INFO=<ID=CODON_CHANGE,Number=1,Type=String,Description="密码子变化（仅CDS区域）">\n')
            out_f.write('##INFO=<ID=AA_CHANGE,Number=1,Type=String,Description="氨基酸变化（仅CDS区域）">\n')
            out_f.write('##INFO=<ID=BRANCH_COUNT,Number=1,Type=Integer,Description="变异所在基因区域的分支数">\n')
            out_f.write(
                '##INFO=<ID=IMPACT,Number=1,Type=String,Description="变异影响程度（HIGH/MODERATE/LOW/MODIFIER）">\n')
            out_f.write(f"{line}\n")
            break


def _process_vcf_data(in_f, out_f, gff_data, transcript_info, ref_genome, branch_counts):
    """处理VCF数据行，添加注释"""
    for line in in_f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        fields = line.split('\t')
        if len(fields) < 8:
            out_f.write(f"{line}\n")
            continue

        chrom, pos, id_, ref, alt, qual, filter_, info = fields[:8]
        other_fields = fields[8:]

        # 注释当前变异
        ann = annotate_variant(
            chrom, int(pos), ref, alt.split(',')[0],  # 取第一个等位基因
            gff_data, transcript_info, ref_genome, branch_counts
        )

        # 构建新的INFO字段
        new_info = [
            f"GENE_ID={ann['gene_id']}",
            f"REGION={ann['region']}",
            f"EFFECT={ann['effect']}",
            f"CODON_CHANGE={ann['codon_change']}",
            f"AA_CHANGE={ann['aa_change']}",
            f"BRANCH_COUNT={ann['branch_count']}",
            f"IMPACT={ann['impact']}"
        ]
        if info != '.':
            new_info = [info] + new_info

        # 写入注释后的行
        output_fields = fields[:7] + [';'.join(new_info)] + other_fields
        out_f.write('\t'.join(output_fields) + '\n')


# --------------------------
# 主函数
# --------------------------
def main():
    if len(sys.argv) != 5:
        print(
            "用法：python complete_variant_annotator.py <输入vcf文件> <gff文件> <gfa文件> <fasta文件(可选，填none则不使用)>")
        print("示例：python complete_variant_annotator.py input.vcf genes.gff assembly.gfa ref_genome.fasta")
        sys.exit(1)

    input_vcf = sys.argv[1]
    gff_file = sys.argv[2]
    gfa_file = sys.argv[3]
    fasta_file = sys.argv[4] if sys.argv[4].lower() != 'none' else None
    output_vcf = f"annotated_with_impact_branches_{input_vcf}"

    try:
        # 1. 提取参考序列
        print("提取参考序列...")
        ref_genome, chrom_map = extract_reference(gfa_file, fasta_file)
        target_chroms = set(chrom_map.values())

        # 2. 解析GFF
        print("解析GFF文件...")
        gff_data, transcript_info = parse_gff(gff_file, target_chroms)

        # 3. 统计分支数
        print("统计基因区域分支数...")
        branch_counts = count_branches(gfa_file, gff_data, target_chroms)

        # 4. 处理VCF并输出
        print("处理VCF并添加注释...")
        process_vcf(input_vcf, output_vcf, gff_data, transcript_info, ref_genome, branch_counts)

        print(f"完成！结果保存至：{output_vcf}")

    except Exception as e:
        print(f"运行出错：{str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
