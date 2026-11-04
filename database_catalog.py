"""Explicit access levels: discovery is not download, analysis, or trained prediction."""
DATABASES = [
 {"name":"PubMed", "role":"生物医学文献与摘要", "access":"complete_paginated_search", "url":"https://pubmed.ncbi.nlm.nih.gov/", "method":"四条疾病相关检索支路分页获取并按PMID去重；单查询超过10,000条时按发表日期分区；报告命中数、获取数和完整性。当前自动化只覆盖PubMed摘要，不等于全文检索或全球穷尽检索；中文库和需授权数据库需另行接入/导入"},
 {"name":"Europe PMC", "role":"生命科学文献、开放全文元数据与预印本", "access":"cursor_paginated_search", "url":"https://europepmc.org/", "method":"四条疾病支路使用REST游标分页取回摘要记录，与PubMed按PMID/DOI去重；能补充部分预印本/开放全文元数据，但不代表全文内容已自动审查；中文及付费库需另接入"},
 {"name":"GEO", "role":"疾病与干预表达数据", "access":"paginated_metadata_and_file_inventory", "url":"https://www.ncbi.nlm.nih.gov/geo/", "method":"按检索式分页获取研究、读取matrix和补充文件目录，提供样本和SRA入口。目录不等于下载或分析；通用分析器支持设计与注释明确的独立/配对log2矩阵及明确标注的原始整数RNA-seq计数探索分析；FASTQ/BAM仍需专用流程"},
 {"name":"BioStudies / ArrayExpress", "role":"补充表达队列", "access":"cursor_metadata_and_file_inventory", "url":"https://www.ebi.ac.uk/biostudies/arrayexpress", "method":"按同义词用游标分页，读取研究文件目录，E-GEOD链接GEO并提供ENA入口；未下载或自动分析全部文件"},
 {"name":"PRIDE / ProteomeXchange", "role":"蛋白组和磷酸化组", "access":"paginated_metadata_and_file_inventory", "url":"https://www.ebi.ac.uk/pride/", "method":"按疾病同义词分页查PXD并按项目去重，逐项目读取文件目录；未下载或解析蛋白定量表及原始质谱"},
 {"name":"Open Targets", "role":"疾病—靶点已知证据", "access":"association_api", "url":"https://platform.opentargets.org/", "method":"仅在疾病名称与规范名/同义词精确匹配后检索关联靶点；近似候选需人工确认本体范围，不向宽泛疾病继承靶点；总分不是因果概率"},
 {"name":"GDC / TCGA", "role":"肿瘤多组学", "access":"project_metadata_api", "url":"https://portal.gdc.cancer.gov/", "method":"仅肿瘤适用；检索项目再按患者ID选病例与对照，STAR计数用DESeq2/edgeR；受控数据需要授权；当前仅检索项目"},
 {"name":"LINCS / CMap / CLUE", "role":"药物及遗传扰动转录响应", "access":"bulk_or_credentials_required", "url":"https://clue.io/data/CMap2020", "method":"导入Level5签名与剂量、时间、细胞系信息；疾病上下调基因做反向连接评分和随机基因集校准；未下载全量签名，不输出真实CMap分数"},
 {"name":"DepMap / PRISM", "role":"CRISPR依赖性与药敏验证", "access":"bulk_required", "url":"https://depmap.org/portal/data_page/", "method":"固定版本下载Chronos、模型注释、PRISM反应和剂量信息；按ModelID及化学标识整合；主要肿瘤模型，不能外推所有疾病"},
 {"name":"ChEMBL", "role":"化合物—靶点机制与活性", "access":"mechanism_api", "url":"https://www.ebi.ac.uk/chembl/", "method":"人类单蛋白靶点精确映射后查mechanism；当前不将其当成疾病疗效；后续需活性、暴露与选择性筛选"},
 {"name":"HERB 2.0", "role":"中药—成分—方剂—疾病实验和临床线索", "access":"public_read_api_and_GEO_experiment_catalog", "url":"http://herb.ac.cn/v2/", "method":"平台只读调用HERB公开API解析精确疾病实体和药材/成分/方剂关联；共享靶点统计关联作为外部数据库先验独立展示，不当作本平台组学预测分数或疗效证据。HERB扰动实验目录可作为实测签名来源；平台自有签名连接分数需要取得完整匹配条件的扰动差异表后另行计算"},
 {"name":"CELLxGENE Census / Human Cell Atlas", "role":"细胞类型和组织背景", "access":"dedicated_worker_required", "url":"https://chanzuckerberg.github.io/cellxgene-census/", "method":"按疾病、组织、供体查询，原始计数供体级pseudobulk；细胞不是独立受试者；Census聚合记录需与原始GEO研究去重"},
 {"name":"GWAS Catalog / eQTL Catalogue", "role":"遗传定位和因果优先级", "access":"harmonized_summary_data_required", "url":"https://www.ebi.ac.uk/gwas/", "method":"下载匹配祖源的GWAS/eQTL/pQTL汇总数据，等位基因校正、共定位、必要时MR；关联P值或最近基因不能直接当因果靶点"},
 {"name":"GTEx", "role":"正常组织表达与eQTL背景", "access":"reference_download_required", "url":"https://gtexportal.org/", "method":"评估候选靶点在组织中是否表达；不能将GTEx和TCGA直接当作无批次差异的健康/疾病对照"},
 {"name":"Reactome / GO via g:Profiler", "role":"基因集功能富集", "access":"enrichment_api", "url":"https://biit.cs.ut.ee/gprofiler/", "method":"用实际被检测基因为背景，上下调分开并校正多重检验；富集不等同通路活性"}
]
