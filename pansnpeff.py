import sys
import re
from collections import defaultdict
from Bio.Seq import Seq  # 用于序列操作（如反向互补、翻译）
from Bio import SeqIO   # 用于读取FASTA格式文件

# --------------------------
# 变异影响程度映射表（参考SnpEff标准）
# 说明：根据变异类型将其划分为4个影响等级，用于快速评估变异的潜在生物学意义
# --------------------------
IMPACT_MAP = {
    # 高影响：直接破坏蛋白质功能（可能导致蛋白失活）
    'nonsense_mutation': 'HIGH',         # 无义突变（提前终止密码子）
    'frameshift_variant': 'HIGH',        # 移码突变（改变阅读框）
    'nonstop_mutation': 'HIGH',          # 无终止突变（终止密码子变为氨基酸）
    'start_lost': 'HIGH',                # 起始密码子丢失
    'stop_lost': 'HIGH',                 # 终止密码子丢失
    
    # 中影响：可能改变蛋白质功能（氨基酸序列改变）
    'missense_mutation': 'MODERATE',     # 错义突变（氨基酸改变）
    'inframe_insertion': 'MODERATE',     # 框内插入（插入长度为3的倍数，不改变阅读框）
    'inframe_deletion': 'MODERATE',      # 框内缺失（缺失长度为3的倍数，不改变阅读框）
    'splice_acceptor_variant': 'MODERATE',# 剪接受体变异（影响mRNA剪接）
    'splice_donor_variant': 'MODERATE',  # 剪接供体变异（影响mRNA剪接）
    
    # 低影响：对蛋白质功能影响较小（氨基酸序列不变）
    'synonymous_variant': 'LOW',         # 同义突变（密码子改变但氨基酸不变）
    'coding_sequence_variant': 'LOW',    # 编码区变异（不改变氨基酸序列）
    
    # 修饰性：不直接影响蛋白质序列（可能影响基因表达调控）
    'intergenic_variant': 'MODIFIER',    # 基因间区变异
    'intron_variant': 'MODIFIER',        # 内含子变异
    'exon_variant': 'MODIFIER',          # 外显子变异（非编码区）
    'genic_variant': 'MODIFIER',         # 基因区变异（非编码区）
    'reference_mismatch': 'MODIFIER',    # 参考序列不匹配（变异与参考序列不一致）
    'incomplete_codon_variant': 'MODIFIER'# 不完整密码子变异（位于密码子不完整区域）
}

# --------------------------
# 1. 参考序列提取模块
# 功能：从GFA或FASTA文件中提取参考基因组序列，优先使用GFA
# --------------------------
def extract_reference(gfa_file, fasta_file=None):
    """
    从GFA文件提取参考序列，若提取失败且提供了FASTA文件，则从FASTA提取
    
    参数：
        gfa_file (str): GFA文件路径（泛基因组组装图）
        fasta_file (str, optional): FASTA文件路径（传统参考基因组），默认为None
    
    返回：
        tuple: (ref_genome, chrom_map)
            - ref_genome: 字典，键为参考序列名，值为对应的DNA序列（字符串）
            - chrom_map: 字典，记录参考序列名的映射关系（通常为{name: name}）
    """
    # 从GFA提取参考序列
    ref_from_gfa, chrom_map = _extract_from_gfa(gfa_file)
    
    # 若GFA未提取到序列，且提供了FASTA，则从FASTA提取
    if not ref_from_gfa and fasta_file:
        ref_from_gfa = _extract_from_fasta(fasta_file)
        chrom_map = {k: k for k in ref_from_gfa.keys()}  # 序列名映射（自身映射）
    
    # 若仍未提取到序列，抛出错误
    if not ref_from_gfa:
        raise ValueError("无法从GFA或FASTA提取参考序列，请检查文件格式")
    
    return ref_from_gfa, chrom_map

def _extract_from_gfa(gfa_file):
    """
    从GFA文件的S行（片段序列）和P行（路径）提取参考序列
    
    逻辑：
        1. 先从S行收集所有片段的序列
        2. 再从P行解析参考序列的组成（片段顺序和方向）
        3. 按顺序拼接片段序列（负链片段需反向互补）
    
    参数：
        gfa_file (str): GFA文件路径
    
    返回：
        tuple: (ref_genome, chrom_map)，同extract_reference
    """
    seg_sequences = {}  # 存储片段ID到序列的映射（{seg_id: sequence}）
    
    # 第一步：读取S行，收集所有片段的序列
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # 跳过空行
            fields = line.split('\t')  # GFA每行以制表符分隔
            # S行格式：S\tseg_id\tsequence\t[标签...]
            if fields[0] == 'S' and len(fields) >= 3:
                seg_id = fields[1]          # 片段ID
                seg_seq = fields[2].upper() # 片段序列（转为大写）
                seg_sequences[seg_id] = seg_seq
    
    ref_genome = {}  # 存储参考序列（{ref_name: full_sequence}）
    
    # 第二步：读取P行，拼接参考序列
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            # P行格式：P\tref_name\tseg1+,seg2-,...\t[标签...]
            if fields[0] == 'P' and len(fields) >= 3:
                ref_name = fields[1]              # 参考序列名（如染色体ID）
                seg_list = fields[2].split(',')   # 片段列表（含方向，如seg1+）
                full_seq = []                     # 用于拼接完整参考序列
                
                for seg in seg_list:
                    # 解析片段ID和方向（+表示正向，-表示反向）
                    if seg.endswith(('+', '-')):
                        seg_id = seg[:-1]         # 去掉方向符号，保留片段ID
                        orientation = seg[-1]     # 方向（+/-）
                    else:
                        seg_id = seg              # 无方向符号，默认为正向
                        orientation = '+'
                    
                    # 跳过未找到序列的片段
                    if seg_id not in seg_sequences:
                        continue
                    
                    # 获取片段序列，负链需做反向互补
                    seg_seq = seg_sequences[seg_id]
                    if orientation == '-':
                        # 反向互补（Seq对象的reverse_complement方法）
                        seg_seq = str(Seq(seg_seq).reverse_complement())
                    
                    full_seq.append(seg_seq)  # 拼接序列
                
                # 保存完整参考序列
                ref_genome[ref_name] = ''.join(full_seq)
                print(f"从GFA P行提取参考序列：{ref_name}（长度：{len(ref_genome[ref_name])}bp）")
    
    return ref_genome, {k: k for k in ref_genome.keys()}

def _extract_from_fasta(fasta_file):
    """
    从FASTA文件提取参考序列（当GFA提取失败时使用）
    
    参数：
        fasta_file (str): FASTA文件路径
    
    返回：
        dict: 键为序列ID，值为对应的DNA序列（字符串）
    """
    ref_genome = {}
    # 使用Bio.SeqIO解析FASTA，record.id为序列ID，record.seq为序列
    for record in SeqIO.parse(fasta_file, 'fasta'):
        ref_genome[record.id] = str(record.seq).upper()  # 转为大写字符串
    return ref_genome

# --------------------------
# 2. GFF文件解析模块
# 功能：解析GFF文件，提取基因、外显子、CDS等结构的位置信息
# --------------------------
def parse_gff(gff_file, target_chroms):
    """
    解析GFF文件，提取基因、转录本、外显子、CDS的位置和关联信息
    
    参数：
        gff_file (str): GFF文件路径
        target_chroms (set): 需要保留的染色体ID集合（与参考序列匹配）
    
    返回：
        tuple: (gff_data, transcript_info)
            - gff_data: 嵌套字典，结构为{gff_seqid: {gene_id: {feature_type: [位置信息...]}}}
            - transcript_info: 字典，存储转录本信息{transcript_id: {gene_id, strand, cds_list}}
    """
    # gff_data结构说明：
    # 第一层键：染色体ID（与参考序列匹配）
    # 第二层键：基因ID
    # 第三层键：特征类型（gene/exon/CDS）
    # 值：列表，每个元素为该特征的位置信息（元组）
    gff_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    
    # transcript_info结构说明：
    # 键：转录本ID
    # 值：字典，包含基因ID（gene_id）、链方向（strand）、CDS列表（cds_list）
    transcript_info = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue  # 跳过注释行（#开头）
            fields = line.split('\t')
            if len(fields) != 9:
                continue  # GFF文件每行应包含9个字段，跳过格式错误的行
            
            # 解析GFF的9个字段
            seqid = fields[0]       # 染色体ID（需与参考序列匹配）
            source = fields[1]      # 注释来源（如ENSEMBL）
            feature = fields[2]     # 特征类型（gene/exon/CDS等）
            start = int(fields[3])  # 起始位置（1-based）
            end = int(fields[4])    # 结束位置（1-based）
            score = fields[5]       # 得分（通常为.）
            strand = fields[6]      # 链方向（+/-）
            phase = fields[7]       # 相位（仅CDS有效，0/1/2/.）
            attr = fields[8]        # 属性字段（如ID=gene123;Parent=transcript456）
            
            # 只保留目标染色体的信息（避免处理无关序列）
            if seqid not in target_chroms:
                continue
            
            # 解析属性字段（转为字典，如{'ID': 'gene123', 'Parent': 'transcript456'}）
            attr_dict = _parse_attributes(attr)
            
            # 根据特征类型分别处理
            if feature == 'gene':
                _process_gene(gff_data, seqid, attr_dict, start, end, strand)
            elif feature == 'transcript':
                _process_transcript(transcript_info, attr_dict)
            elif feature == 'exon':
                _process_exon(gff_data, seqid, transcript_info, attr_dict, start, end, strand)
            elif feature == 'CDS':
                _process_cds(gff_data, seqid, transcript_info, attr_dict, start, end, strand, phase)
    
    return gff_data, transcript_info

def _parse_attributes(attr_str):
    """
    解析GFF的属性字段（如"ID=gene123;Parent=transcript456"）
    
    参数：
        attr_str (str): 属性字段字符串
    
    返回：
        dict: 解析后的属性字典（键为属性名，值为属性值）
    """
    attr_dict = {}
    # 按分号分隔属性，处理每个键值对
    for item in attr_str.split(';'):
        item = item.strip()
        if '=' in item:
            # 按第一个等号分割（避免值中包含等号）
            k, v = item.split('=', 1)
            attr_dict[k.strip()] = v.strip()
    return attr_dict

def _process_gene(gff_data, seqid, attr_dict, start, end, strand):
    """处理gene类型特征，记录基因位置"""
    gene_id = attr_dict.get('ID')  # 基因ID（从属性字段的ID获取）
    if gene_id:
        # 存储基因位置信息：(起始位置, 结束位置, 链方向)
        gff_data[seqid][gene_id]['gene'].append((start, end, strand))

def _process_transcript(transcript_info, attr_dict):
    """处理transcript类型特征，记录转录本与基因的关联"""
    transcript_id = attr_dict.get('ID')  # 转录本ID
    gene_id = attr_dict.get('Parent')    # 所属基因ID（从Parent属性获取）
    if transcript_id and gene_id:
        # 初始化转录本信息：基因ID、链方向（默认+，后续会被CDS/exon覆盖）、CDS列表
        transcript_info[transcript_id] = {
            'gene_id': gene_id,
            'strand': '+',
            'cds_list': []  # 存储CDS片段的位置和相位
        }

def _process_exon(gff_data, seqid, transcript_info, attr_dict, start, end, strand):
    """处理exon类型特征，记录外显子位置及所属转录本"""
    transcript_id = attr_dict.get('Parent')  # 所属转录本ID（从Parent获取）
    # 跳过无转录本ID或转录本未记录的外显子
    if not transcript_id or transcript_id not in transcript_info:
        return
    gene_id = transcript_info[transcript_id]['gene_id']  # 所属基因ID
    # 存储外显子信息：(起始位置, 结束位置, 链方向, 转录本ID)
    gff_data[seqid][gene_id]['exon'].append((start, end, strand, transcript_id))
    transcript_info[transcript_id]['strand'] = strand  # 更新转录本的链方向

def _process_cds(gff_data, seqid, transcript_info, attr_dict, start, end, strand, phase):
    """处理CDS类型特征，记录CDS位置、相位及所属转录本"""
    transcript_id = attr_dict.get('Parent')  # 所属转录本ID
    if not transcript_id or transcript_id not in transcript_info:
        return
    gene_id = transcript_info[transcript_id]['gene_id']  # 所属基因ID
    # 存储CDS信息：(起始位置, 结束位置, 链方向, 转录本ID, 相位)
    gff_data[seqid][gene_id]['CDS'].append((start, end, strand, transcript_id, phase))
    # 记录CDS到转录本的CDS列表（用于后续拼接完整CDS序列）
    transcript_info[transcript_id]['cds_list'].append((start, end, phase))
    # 按起始位置排序CDS（确保拼接顺序正确）
    transcript_info[transcript_id]['cds_list'].sort(key=lambda x: x[0])
    transcript_info[transcript_id]['strand'] = strand  # 更新转录本的链方向

# --------------------------
# 3. 分支数统计模块（修改后不依赖PO标签）
# 功能：基于GFA的片段连接关系，统计基因区域的分支数（反映泛基因组结构多样性）
# --------------------------
def count_branches(gfa_file, gff_data, ref_names):
    """
    统计每个基因区域（gene/exon/CDS）的分支数（片段连接数）
    
    参数：
        gfa_file (str): GFA文件路径
        gff_data (dict): 从GFF解析的基因结构数据（见parse_gff返回值）
        ref_names (set): 参考序列名集合（染色体ID）
    
    返回：
        dict: 分支数统计结果，结构为{ref_name: {gene_id: {feature_type: count}}}
    """
    # 第一步：提取所有片段的长度（用于计算位置）
    seg_lengths = _parse_segment_lengths(gfa_file)
    
    # 第二步：解析GFA，计算片段位置、收集连接关系
    path_segments, seg_coords, links, seg_to_ref = _parse_gfa(gfa_file, seg_lengths)
    
    # 第三步：将基因区域映射到对应的GFA片段
    region_segs = _map_regions_to_segments(gff_data, seg_coords, seg_to_ref, ref_names)
    
    # 第四步：统计每个区域的分支数
    return _count_branch_numbers(region_segs, links)

def _parse_segment_lengths(gfa_file):
    """
    从GFA的S行提取片段长度（优先用LN标签，无标签则用序列长度）
    
    参数：
        gfa_file (str): GFA文件路径
    
    返回：
        dict: {seg_id: length}，片段ID到长度的映射
    """
    seg_lengths = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if fields[0] == 'S' and len(fields) >= 3:
                seg_id = fields[1]      # 片段ID
                seg_seq = fields[2]     # 片段序列
                
                # 优先从LN标签获取长度（格式：LN:i:100）
                ln_tag = next((t for t in fields[3:] if t.startswith('LN:i:')), None)
                if ln_tag:
                    # 提取长度值（如LN:i:100 → 100）
                    seg_lengths[seg_id] = int(ln_tag.split(':')[2])
                else:
                    # 无LN标签时，用序列长度作为片段长度
                    seg_lengths[seg_id] = len(seg_seq)
    return seg_lengths

def _parse_gfa(gfa_file, seg_lengths):
    """
    解析GFA文件，不依赖PO标签，通过P行和片段长度计算片段位置
    
    核心逻辑：
        1. 从P行解析参考序列的片段组成（顺序和方向）
        2. 按片段顺序累加长度，计算每个片段在参考序列中的start和end（0-based）
        3. 从L行收集片段间的连接关系（用于统计分支数）
    
    参数：
        gfa_file (str): GFA文件路径
        seg_lengths (dict): 片段ID到长度的映射（见_parse_segment_lengths）
    
    返回：
        tuple: (path_segments, seg_coords, links, seg_to_ref)
            - path_segments: {ref_name: {seg_id1, seg_id2, ...}} 参考序列包含的片段
            - seg_coords: {seg_id: (start, end)} 片段在参考序列中的位置（0-based）
            - links: [(from_seg, to_seg), ...] 片段间的连接关系
            - seg_to_ref: {seg_id: ref_name} 片段所属的参考序列
    """
    path_segments = defaultdict(set)  # 参考序列包含的片段
    seg_coords = {}                   # 片段位置（start, end）
    links = []                        # 片段连接关系
    seg_to_ref = {}                   # 片段所属参考序列
    
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            
            # 处理P行：解析参考序列的片段组成，并计算位置
            if fields[0] == 'P' and len(fields) >= 3:
                ref_name = fields[1]              # 参考序列名
                segments = fields[2].split(',')   # 片段列表（如seg1+,seg2-）
                current_pos = 0                   # 当前位置（0-based，初始为0）
                
                for seg in segments:
                    # 解析片段ID和方向（+/-）
                    if seg.endswith(('+', '-')):
                        seg_id = seg[:-1]  # 片段ID（去掉方向符号）
                        orientation = seg[-1]  # 方向（不影响位置计算，仅影响序列方向）
                    else:
                        seg_id = seg
                        orientation = '+'
                    
                    # 跳过无长度信息的片段
                    if seg_id not in seg_lengths:
                        continue
                    
                    # 计算片段在参考序列中的位置（0-based）
                    seg_len = seg_lengths[seg_id]       # 片段长度
                    seg_start = current_pos              # 片段起始位置
                    seg_end = current_pos + seg_len - 1  # 片段结束位置（闭区间）
                    
                    # 记录片段信息
                    seg_coords[seg_id] = (seg_start, seg_end)  # 位置
                    seg_to_ref[seg_id] = ref_name              # 所属参考序列
                    path_segments[ref_name].add(seg_id)         # 参考序列包含的片段
                    
                    # 更新当前位置（累加片段长度，为下一个片段计算做准备）
                    current_pos += seg_len
            
            # 处理L行：记录片段间的连接关系（用于统计分支数）
            elif fields[0] == 'L' and len(fields) >= 5:
                # L行格式：L\tfrom_seg\tfrom_orient\tto_seg\tto_orient\toverlap
                from_seg = fields[1]  # 起始片段
                to_seg = fields[3]    # 目标片段
                links.append((from_seg, to_seg))  # 记录连接关系
    
    return path_segments, seg_coords, links, seg_to_ref

def _map_regions_to_segments(gff_data, seg_coords, seg_to_ref, ref_names):
    """
    将基因区域（gene/exon/CDS）映射到对应的GFA片段（判断位置重叠）
    
    参数：
        gff_data (dict): 基因结构数据
        seg_coords (dict): 片段位置（start, end）
        seg_to_ref (dict): 片段所属参考序列
        ref_names (set): 参考序列名集合
    
    返回：
        dict: 区域-片段映射，结构为{ref_name: {gene_id: {feature_type: {seg_id1, seg_id2, ...}}}}
    """
    # 初始化映射关系（嵌套字典）
    region_segs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    
    for ref_name in ref_names:
        if ref_name not in gff_data:
            continue  # 跳过无基因结构的参考序列
        
        # 遍历该参考序列下的所有基因
        genes = gff_data[ref_name]
        for gene_id, features in genes.items():
            # 遍历基因的所有特征（gene/exon/CDS）
            for feature_type, regions in features.items():
                # 遍历每个特征的位置区域
                for region in regions:
                    if len(region) < 2:
                        continue  # 跳过位置信息不完整的区域
                    
                    # 提取区域的起始和结束位置（GFF是1-based）
                    reg_start, reg_end = region[0], region[1]
                    
                    # 寻找与该区域重叠的片段
                    for seg_id, (seg_start, seg_end) in seg_coords.items():
                        # 片段必须属于当前参考序列
                        if seg_to_ref.get(seg_id) != ref_name:
                            continue
                        
                        # 判断区域（1-based）与片段（0-based）是否重叠
                        # 转换逻辑：1-based的[reg_start, reg_end] → 0-based的[reg_start-1, reg_end-1]
                        # 重叠条件：片段结束 >= 区域起始-1 且 片段起始 <= 区域结束-1
                        if not (seg_end < reg_start - 1 or seg_start > reg_end - 1):
                            # 记录重叠的片段
                            region_segs[ref_name][gene_id][feature_type].add(seg_id)
    
    return region_segs

def _count_branch_numbers(region_segs, links):
    """
    统计每个基因区域的分支数（与区域重叠的片段所涉及的连接数）
    
    参数：
        region_segs (dict): 区域-片段映射（见_map_regions_to_segments）
        links (list): 片段连接关系[(from_seg, to_seg), ...]
    
    返回：
        dict: 分支数统计结果，结构同count_branches返回值
    """
    branch_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    
    # 遍历每个参考序列的基因区域
    for ref_name, genes in region_segs.items():
        for gene_id, features in genes.items():
            for feature_type, segs in features.items():
                count = 0
                # 统计与该区域片段相关的连接数
                for (from_seg, to_seg) in links:
                    # 若连接的起始或目标片段属于该区域，则计数+1
                    if from_seg in segs or to_seg in segs:
                        count += 1
                # 记录统计结果
                branch_counts[ref_name][gene_id][feature_type] = count
    
    return branch_counts

# --------------------------
# 4. 变异注释核心逻辑
# 功能：注释变异的基因区域、功能影响、分支数和影响程度
# --------------------------
def annotate_variant(chrom, pos, ref, alt, gff_data, transcript_info, ref_genome, branch_counts):
    """
    注释单个变异的详细信息（基因区域、功能影响、分支数等）
    
    参数：
        chrom (str): 变异所在的染色体ID
        pos (int): 变异位置（1-based）
        ref (str): 参考等位基因
        alt (str): 替代等位基因（取第一个）
        gff_data (dict): 基因结构数据
        transcript_info (dict): 转录本信息
        ref_genome (dict): 参考序列
        branch_counts (dict): 分支数统计结果
    
    返回：
        dict: 变异注释结果，包含gene_id/region/effect等信息
    """
    # 初始化注释结果
    annotation = {
        'gene_id': '.',               # 基因ID（.表示无）
        'region': 'intergenic',       # 所在区域（基因间区）
        'effect': 'intergenic_variant',# 变异影响类型
        'codon_change': '.',          # 密码子变化（仅CDS区域）
        'aa_change': '.',             # 氨基酸变化（仅CDS区域）
        'branch_count': 0,            # 分支数
        'impact': 'MODIFIER'          # 影响程度
    }
    
    # 若变异所在染色体无基因结构数据，直接返回基因间区注释
    if chrom not in gff_data:
        annotation['impact'] = IMPACT_MAP.get(annotation['effect'], 'MODIFIER')
        return annotation
    
    # 遍历基因，判断变异是否在基因区域内
    for gene_id, features in gff_data[chrom].items():
        # 检查变异是否在当前基因的范围内
        for (gene_start, gene_end, _) in features.get('gene', []):
            if gene_start <= pos <= gene_end:
                # 变异在基因区域内，更新注释
                annotation['gene_id'] = gene_id
                annotation['region'] = 'genic'
                # 进一步判断具体在基因的哪个子区域（exon/intron/CDS）
                annotation = _check_regions(
                    chrom, pos, ref, alt, gene_id, features, transcript_info,
                    ref_genome, branch_counts, annotation
                )
                # 计算影响程度
                annotation['impact'] = IMPACT_MAP.get(annotation['effect'], 'MODIFIER')
                return annotation
    
    # 若未匹配到任何基因，保持基因间区注释
    annotation['impact'] = IMPACT_MAP.get(annotation['effect'], 'MODIFIER')
    return annotation

def _check_regions(chrom, pos, ref, alt, gene_id, features, transcript_info, ref_genome, branch_counts, annotation):
    """
    检查变异在基因内的具体区域（exon/intron/CDS）并更新注释
    
    参数：
        同上，增加annotation（当前注释结果）
    
    返回：
        dict: 更新后的注释结果
    """
    # 第一步：检查是否在外显子区域
    for region in features.get('exon', []):
        if len(region) < 3:
            continue  # 跳过信息不完整的外显子
        exon_start, exon_end, strand, transcript_id = region[:4]
        if exon_start <= pos <= exon_end:
            annotation['region'] = 'exon'               # 更新区域为外显子
            annotation['effect'] = 'exon_variant'       # 更新影响类型
            # 获取外显子区域的分支数
            annotation['branch_count'] = branch_counts.get(chrom, {}).get(gene_id, {}).get('exon', 0)
            return annotation
    
    # 第二步：检查是否在CDS区域（编码区）
    for region in features.get('CDS', []):
        if len(region) < 4:
            continue  # 跳过信息不完整的CDS
        cds_start, cds_end, strand, transcript_id = region[:4]
        if cds_start <= pos <= cds_end:
            annotation['region'] = 'CDS'                # 更新区域为CDS
            # 获取CDS区域的分支数
            annotation['branch_count'] = branch_counts.get(chrom, {}).get(gene_id, {}).get('CDS', 0)
            # 预测CDS区域变异对氨基酸的影响
            return _predict_cds_effect(
                chrom, pos, ref, alt, strand, transcript_id, transcript_info, ref_genome, annotation
            )
    
    # 第三步：若不在exon/CDS，则为内含子区域
    annotation['region'] = 'intron'                  # 更新区域为内含子
    annotation['effect'] = 'intron_variant'          # 更新影响类型
    # 获取内含子区域的分支数
    annotation['branch_count'] = branch_counts.get(chrom, {}).get(gene_id, {}).get('intron', 0)
    return annotation

def _predict_cds_effect(chrom, pos, ref, alt, strand, transcript_id, transcript_info, ref_genome, annotation):
    """
    预测CDS区域变异对氨基酸的影响（如错义、同义、移码等）
    
    参数：
        同上
    
    返回：
        dict: 更新后的注释结果（含密码子和氨基酸变化）
    """
    # 获取转录本信息（若不存在，返回当前注释）
    transcript = transcript_info.get(transcript_id)
    if not transcript:
        return annotation
    
    cds_list = transcript['cds_list']       # CDS片段列表（已排序）
    gene_strand = transcript['strand']      # 基因链方向（+/-）
    
    # 拼接完整的CDS序列（用于分析密码子变化）
    full_cds = _get_full_cds_sequence(chrom, cds_list, ref_genome)
    
    # 计算变异在CDS序列中的位置（0-based）
    variant_in_cds_pos = _get_variant_cds_position(pos, cds_list)
    if variant_in_cds_pos is None:
        return annotation  # 变异不在CDS区域（理论上不应发生）
    
    # 处理负链：变异位置反转，参考和替代等位基因取反向互补
    if gene_strand == '-':
        # 负链上的位置 = CDS长度 - 1 - 正向位置
        variant_in_cds_pos = len(full_cds) - 1 - variant_in_cds_pos
        # 反向互补（A↔T，C↔G）
        ref = str(Seq(ref).reverse_complement())
        alt = str(Seq(alt).reverse_complement())
    
    # 根据变异类型（SNV/插入/缺失）预测影响
    if len(ref) == 1 and len(alt) == 1:
        # 单核苷酸变异（SNV）
        return _annotate_snv(full_cds, variant_in_cds_pos, ref, alt, annotation)
    elif len(ref) < len(alt):
        # 插入变异：判断是否移码（长度变化不是3的倍数）
        if (len(alt) - len(ref)) % 3 != 0:
            annotation['effect'] = 'frameshift_variant'  # 移码突变
        else:
            annotation['effect'] = 'inframe_insertion'   # 框内插入
    elif len(ref) > len(alt):
        # 缺失变异：判断是否移码
        if (len(ref) - len(alt)) % 3 != 0:
            annotation['effect'] = 'frameshift_variant'  # 移码突变
        else:
            annotation['effect'] = 'inframe_deletion'    # 框内缺失
    
    return annotation

def _get_full_cds_sequence(chrom, cds_list, ref_genome):
    """
    拼接完整的CDS序列（按位置排序的CDS片段拼接）
    
    参数：
        chrom (str): 染色体ID
        cds_list (list): CDS片段列表[(start, end, phase), ...]
        ref_genome (dict): 参考序列
    
    返回：
        str: 完整的CDS序列
    """
    cds_sequence = []
    for (start, end, _) in cds_list:
        # GFF是1-based，参考序列是0-based，因此取[start-1:end]
        cds_seq = ref_genome[chrom][start-1:end]
        cds_sequence.append(cds_seq)
    return ''.join(cds_sequence)

def _get_variant_cds_position(pos, cds_list):
    """
    计算变异在完整CDS序列中的位置（0-based）
    
    参数：
        pos (int): 变异在染色体上的位置（1-based）
        cds_list (list): CDS片段列表[(start, end, phase), ...]
    
    返回：
        int: 变异在CDS序列中的位置（0-based），若不在CDS中返回None
    """
    cds_offset = 0  # CDS序列的偏移量（累加前面CDS片段的长度）
    for (start, end, _) in cds_list:
        # 检查变异是否在当前CDS片段内
        if start <= pos <= end:
            # 位置 = 偏移量 + 变异在当前CDS中的相对位置
            return cds_offset + (pos - start)
        # 更新偏移量（当前CDS片段的长度：end - start + 1）
        cds_offset += (end - start + 1)
    return None  # 变异不在任何CDS片段内

def _annotate_snv(full_cds, pos, ref, alt, annotation):
    """
    注释单核苷酸变异（SNV）的密码子和氨基酸变化
    
    参数：
        full_cds (str): 完整CDS序列
        pos (int): 变异在CDS中的位置（0-based）
        ref (str): 参考等位基因
        alt (str): 替代等位基因
        annotation (dict): 当前注释结果
    
    返回：
        dict: 更新后的注释结果（含密码子和氨基酸变化）
    """
    # 检查参考等位基因是否与CDS序列一致（不一致则标记为参考不匹配）
    if full_cds[pos] != ref:
        annotation['effect'] = 'reference_mismatch'
        return annotation
    
    # 计算变异所在的密码子信息
    codon_index = pos // 3          # 密码子索引（0-based）
    codon_pos = pos % 3             # 变异在密码子中的位置（0-2）
    codon_start = codon_index * 3   # 密码子起始位置（0-based）
    codon_end = codon_start + 3     # 密码子结束位置（不包含）
    original_codon = full_cds[codon_start:codon_end]  # 原始密码子
    
    # 处理不完整密码子（CDS序列末尾可能不足3个碱基）
    if len(original_codon) < 3:
        annotation['effect'] = 'incomplete_codon_variant'
        return annotation
    
    # 计算变异后的密码子
    mutated_codon = list(original_codon)  # 转为列表便于修改
    mutated_codon[codon_pos] = alt        # 替换变异位置的碱基
    mutated_codon = ''.join(mutated_codon)# 转回字符串
    
    # 翻译密码子为氨基酸（*表示终止密码子）
    original_aa = str(Seq(original_codon).translate())   # 原始氨基酸
    mutated_aa = str(Seq(mutated_codon).translate())     # 变异后氨基酸
    
    # 根据氨基酸变化确定变异类型
    if mutated_aa == '*':
        # 变异后变为终止密码子（无义突变）
        annotation['effect'] = 'nonsense_mutation'
    elif original_aa == '*':
        # 原始为终止密码子，变异后变为氨基酸（无终止突变）
        annotation['effect'] = 'nonstop_mutation'
    elif original_aa != mutated_aa:
        # 氨基酸发生改变（错义突变）
        annotation['effect'] = 'missense_mutation'
    else:
        # 氨基酸不变（同义突变）
        annotation['effect'] = 'synonymous_variant'
    
    # 记录密码子和氨基酸变化
    annotation['codon_change'] = f"{original_codon}→{mutated_codon}"
    annotation['aa_change'] = f"{original_aa}{codon_index+1}{mutated_aa}"  # 1-based索引
    return annotation

# --------------------------
# 5. VCF文件处理模块
# 功能：读取输入VCF，添加注释，输出带注释的VCF
# --------------------------
def process_vcf(input_vcf, output_vcf, gff_data, transcript_info, ref_genome, branch_counts):
    """
    处理VCF文件，为每个变异添加注释并输出新VCF
    
    参数：
        input_vcf (str): 输入VCF文件路径
        output_vcf (str): 输出VCF文件路径
        其他参数：同annotate_variant
    """
    with open(input_vcf, 'r') as in_f, open(output_vcf, 'w') as out_f:
        # 处理VCF头部（添加新注释字段定义）
        _write_vcf_header(in_f, out_f)
        
        # 处理VCF数据行（添加注释）
        _process_vcf_data(in_f, out_f, gff_data, transcript_info, ref_genome, branch_counts)

def _write_vcf_header(in_f, out_f):
    """
    写入VCF头部，保留原始头部并添加新注释字段的定义
    
    参数：
        in_f: 输入VCF文件句柄（读）
        out_f: 输出VCF文件句柄（写）
    """
    for line in in_f:
        line = line.strip()
        # 保留原始##开头的头部行
        if line.startswith('##'):
            out_f.write(f"{line}\n")
        # 遇到#CHROM行（列名行），先添加新字段定义再写入列名
        elif line.startswith('#CHROM'):
            # 定义新增的INFO字段（符合VCF规范）
            out_f.write('##INFO=<ID=GENE_ID,Number=1,Type=String,Description="变异所在的基因ID">\n')
            out_f.write('##INFO=<ID=REGION,Number=1,Type=String,Description="变异所在的基因区域（intergenic/intron/exon/CDS）">\n')
            out_f.write('##INFO=<ID=EFFECT,Number=1,Type=String,Description="变异的生物学影响">\n')
            out_f.write('##INFO=<ID=CODON_CHANGE,Number=1,Type=String,Description="密码子变化（仅CDS区域）">\n')
            out_f.write('##INFO=<ID=AA_CHANGE,Number=1,Type=String,Description="氨基酸变化（仅CDS区域）">\n')
            out_f.write('##INFO=<ID=BRANCH_COUNT,Number=1,Type=Integer,Description="变异所在基因区域的分支数">\n')
            out_f.write('##INFO=<ID=IMPACT,Number=1,Type=String,Description="变异影响程度（HIGH/MODERATE/LOW/MODIFIER）">\n')
            # 写入列名行
            out_f.write(f"{line}\n")
            break  # 头部处理完毕，跳出循环

def _process_vcf_data(in_f, out_f, gff_data, transcript_info, ref_genome, branch_counts):
    """
    处理VCF数据行，为每个变异添加注释
    
    参数：
        同_write_vcf_header，增加注释所需的基因数据等
    """
    for line in in_f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue  # 跳过空行和剩余的注释行
        
        fields = line.split('\t')
        if len(fields) < 8:
            # 字段不足的行直接输出（格式错误）
            out_f.write(f"{line}\n")
            continue
        
        # 解析VCF数据行的关键字段
        chrom = fields[0]    # 染色体ID
        pos = fields[1]      # 位置（字符串，需转为整数）
        id_ = fields[2]      # 变异ID
        ref = fields[3]      # 参考等位基因
        alt = fields[4]      # 替代等位基因（可能多个，取第一个）
        qual = fields[5]     # 质量值
        filter_ = fields[6]  # 过滤标记
        info = fields[7]     # INFO字段（原始注释）
        other_fields = fields[8:]  # 其他字段（如样本基因型）
        
        # 注释当前变异（取第一个替代等位基因）
        ann = annotate_variant(
            chrom, int(pos), ref, alt.split(',')[0],
            gff_data, transcript_info, ref_genome, branch_counts
        )
        
        # 构建新的INFO字段（合并原始INFO和新注释）
        new_info = [
            f"GENE_ID={ann['gene_id']}",
            f"REGION={ann['region']}",
            f"EFFECT={ann['effect']}",
            f"CODON_CHANGE={ann['codon_change']}",
            f"AA_CHANGE={ann['aa_change']}",
            f"BRANCH_COUNT={ann['branch_count']}",
            f"IMPACT={ann['impact']}"
        ]
        # 若原始INFO非空，将其添加到新INFO的前面
        if info != '.':
            new_info = [info] + new_info
        
        # 拼接输出行并写入
        output_fields = fields[:7] + [';'.join(new_info)] + other_fields
        out_f.write('\t'.join(output_fields) + '\n')

# --------------------------
# 主函数
# 功能：解析命令行参数，调用各模块完成注释流程
# --------------------------
def main():
    # 检查命令行参数是否正确
    if len(sys.argv) != 5:
        print("用法：python complete_variant_annotator_with_details.py <输入vcf文件> <gff文件> <gfa文件> <fasta文件(可选，填none则不使用)>")
        print("示例：python complete_variant_annotator_with_details.py input.vcf genes.gff assembly.gfa ref_genome.fasta")
        sys.exit(1)
    
    # 解析命令行参数
    input_vcf = sys.argv[1]       # 输入VCF文件
    gff_file = sys.argv[2]        # GFF文件
    gfa_file = sys.argv[3]        # GFA文件
    # 处理FASTA参数（若为none则不使用）
    fasta_file = sys.argv[4] if sys.argv[4].lower() != 'none' else None
    # 输出文件名为：annotated_with_impact_branches_<输入文件名>
    output_vcf = f"annotated_with_impact_branches_{input_vcf}"
    
    try:
        # 1. 提取参考序列（优先GFA，次选FASTA）
        print("提取参考序列...")
        ref_genome, chrom_map = extract_reference(gfa_file, fasta_file)
        target_chroms = set(chrom_map.values())  # 需要处理的染色体ID
        
        # 2. 解析GFF文件，获取基因结构信息
        print("解析GFF文件...")
        gff_data, transcript_info = parse_gff(gff_file, target_chroms)
        
        # 3. 统计基因区域的分支数（基于GFA）
        print("统计基因区域分支数...")
        branch_counts = count_branches(gfa_file, gff_data, target_chroms)
        
        # 4. 处理VCF文件，添加注释并输出
        print("处理VCF并添加注释...")
        process_vcf(input_vcf, output_vcf, gff_data, transcript_info, ref_genome, branch_counts)
        
        print(f"完成！结果保存至：{output_vcf}")
    
    except Exception as e:
        # 捕获并显示错误信息
        print(f"运行出错：{str(e)}", file=sys.stderr)
        sys.exit(1)

# 程序入口
if __name__ == "__main__":
    main()
