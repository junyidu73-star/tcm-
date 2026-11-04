(() => {
  const $ = (id) => document.getElementById(id);
  const state = { papers: [], decisions: new Map(), searchMeta: null, paperPage: 0, diagnosticPapers: [], lastMeta: null, lastDiagnostic: null };
  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  async function api(path, data) {
    const body=JSON.stringify(data);if(new TextEncoder().encode(body).length>950_000)throw new Error('CSV 经 JSON 编码后超过服务端请求上限；请拆分文件或减小文件体积后重试。');
    const response = await fetch(path, { method: 'POST', credentials: 'same-origin', headers: {'Content-Type':'application/json'}, body });
    const result = await response.json().catch(() => ({}));
    if (!response.ok) throw new Error(result.error || `请求失败 (${response.status})`);
    return result;
  }
  function setStatus(id, text, error=false) { const el=$(id); if(el){if(el.textContent!==text)el.textContent=text;el.classList.toggle('error',error);} }
  async function readCsvInput(id) {
    const file=$(id)?.files?.[0];
    if(!file) throw new Error('请先选择 CSV 文件。');
    if(file.size>700_000) throw new Error('文件请控制在 700 KB 以内，避免 JSON 编码后超过服务端请求上限。');
    return await file.text();
  }
  const promptTemplate = (need) => `你是循证医学文献检索与筛选方案助手。研究者的原始需求如下：\n${need}\n\n请严格按以下步骤生成可复核方案：\n1. 不擅自补充用户未提供的疾病亚型、人群、模型、干预、比较、结局或研究设计；未知项列为待确认，不要把它们设为排除条件。\n2. 拆分研究问题（PICO/PECO 或适合该问题的框架），给出同义词、缩写、拼写变体及主题词；中文需求同时考虑规范英文名和常见缩写。\n3. 输出可直接用于 PubMed 的英文检索式，并说明每个概念块；检索式优先保证召回率，不用过窄限制压掉相关研究。另列建议检索数据库，不得声称已经检索未实际调用的数据库。\n4. 提供题名/摘要初筛及全文筛选的纳入、排除标准，每条标准要可操作；模糊记录标“待核实”，不能自动排除。\n5. 说明需提取的研究设计、样本、物种/组织、干预剂量与时间、对照、结局、数据类型及偏倚风险字段。\n6. 给出至少 3 篇可作检索式哨兵校验的已知相关论文应如何选取（不能编造引用）；命中数及检索日期由真实检索后填写。\n7. 明确全文版权/数据库访问限制、筛选模型局限、人工双人复核和分歧裁决要求。\n\n回答中不要编造文献、数据库命中数、DOI/PMID 或已完成的筛选结果。`;
  function csvEscape(value) { return `"${String(value ?? '').replace(/"/g,'""')}"`; }
  function download(name, text, mime='text/csv;charset=utf-8') {
    const url=URL.createObjectURL(new Blob(['\uFEFF',text],{type:mime})); const a=document.createElement('a'); a.href=url; a.download=name; a.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  function renderPlan(plan, query) {
    $('literaturePlan').innerHTML=`<div class="research-prompt"><b>发给模型的精确任务提示词</b><pre>${esc(plan)}</pre><button type="button" class="button button-secondary" id="copyLiteraturePrompt">复制提示词</button></div><label class="field-label" for="editablePubmedQuery">可编辑 PubMed 检索式</label><textarea id="editablePubmedQuery" class="text-input text-area" rows="4" placeholder="模型生成后可在此修改；也可手动输入 PubMed 查询式。">${esc(query)}</textarea><div class="research-prompt-note">检索方案需要研究者确认；检索式通过真实 PubMed 接口执行后，才会展示记录。系统不把大模型生成的建议当作已检索结果。</div>`;
    $('copyLiteraturePrompt')?.addEventListener('click',async()=>{try{await navigator.clipboard.writeText(plan);setStatus('literatureStatus','提示词已复制，可用于外部模型或存档。')}catch{setStatus('literatureStatus','浏览器未授予剪贴板权限，请手动复制上方提示词。',true)}});$('editablePubmedQuery')?.addEventListener('input',()=>{$('runLiteratureSearch').disabled=!$('editablePubmedQuery').value.trim()});
  }
  async function makeLiteraturePlan() {
    const need=$('literatureNeed').value.trim(); if(!need){setStatus('literatureStatus','请先输入研究主题或文献类型。',true);return;}
    const prompt=promptTemplate(need); renderPlan(prompt,''); $('runLiteratureSearch').disabled=true; setStatus('literatureStatus','正在请求已配置模型生成检索式…');
    try {
      const response=await api('/api/chat',{messages:[{role:'user',content:`请依据以下任务生成 PubMed 检索方案。将完整检索式放入 pubmed_query；在 reply 中输出研究问题、待确认信息、纳入/排除标准、建议提取字段和筛选限制。不要声称已执行检索。\n\n${prompt}`}],literature:[]});
      renderPlan(prompt,response.pubmed_query||''); $('runLiteratureSearch').disabled=!(response.pubmed_query||'').trim();
      setStatus('literatureStatus',response.reply||'方案已生成。请先核对标准，再执行 PubMed 检索。');
    } catch(error) { setStatus('literatureStatus',`模型生成失败：${error.message}。仍可复制提示词，或手动填写 PubMed 检索式后执行。`,true); $('runLiteratureSearch').disabled=true; }
  }
  function renderPapers() {
    const meta=state.searchMeta||{};$('literatureCount').textContent=`PubMed 命中 ${meta.total??state.papers.length} 条；本次取回 ${state.papers.length} 条${meta.complete===false?'（未完整）':''}；导出包含全部当前标签`;
    if(!state.papers.length){$('literatureRecords').innerHTML='<div class="empty-state">没有结果。调整检索词后再试。</div>';return;}
    const perPage=50,start=state.paperPage*perPage,visible=state.papers.slice(start,start+perPage);
    $('literatureRecords').innerHTML=`${meta.complete===false?`<div class="research-warning">本次只取回 ${state.papers.length} 条，数据库报告总命中 ${meta.total??'未知'} 条；受单次上限或接口分页限制，当前列表不完整。可提高上限或按日期拆分检索。原始分页状态已记录。</div>`:''}<div class="research-table"><table><thead><tr><th>题名与来源</th><th>摘要</th><th>人工筛选</th><th>理由</th></tr></thead><tbody>${visible.map((p,offset)=>{const i=start+offset,key=String(p.pmid||i),decision=state.decisions.get(key)||'待核实';return `<tr><td><a href="${esc(p.url||'#')}" target="_blank" rel="noopener noreferrer">${esc(p.title||'无题名')} ↗</a><small>PMID ${esc(p.pmid||'—')} · ${esc(p.year||'年份未知')}</small></td><td>${esc(p.abstract||'无摘要')}</td><td><select class="text-input paper-decision" data-key="${esc(key)}"><option ${decision==='待核实'?'selected':''}>待核实</option><option ${decision==='纳入'?'selected':''}>纳入</option><option ${decision==='排除'?'selected':''}>排除</option></select></td><td><input class="text-input paper-reason" data-key="${esc(key+':reason')}" placeholder="简要理由（全文排除需记录理由）" value="${esc(state.decisions.get(key+':reason')||'')}"></td></tr>`}).join('')}</tbody></table></div><div class="research-page-controls"><button class="button button-secondary" id="paperPrev" ${state.paperPage===0?'disabled':''}>上一页</button><span>显示 ${start+1}–${Math.min(start+perPage,state.papers.length)} / ${state.papers.length}</span><button class="button button-secondary" id="paperNext" ${start+perPage>=state.papers.length?'disabled':''}>下一页</button></div>`;
    document.querySelectorAll('.paper-decision').forEach(el=>el.addEventListener('change',()=>state.decisions.set(el.dataset.key,el.value)));
    document.querySelectorAll('.paper-reason').forEach(el=>el.addEventListener('input',()=>state.decisions.set(el.dataset.key,el.value)));
    $('paperPrev')?.addEventListener('click',()=>{state.paperPage=Math.max(0,state.paperPage-1);renderPapers()});$('paperNext')?.addEventListener('click',()=>{if((state.paperPage+1)*50<state.papers.length){state.paperPage++;renderPapers()}});
  }
  async function runLiterature() {
    const query=$('editablePubmedQuery')?.value.trim(); if(!query){setStatus('literatureStatus','请填写 PubMed 检索式。',true);return;}
    setStatus('literatureStatus','正在通过 PubMed 检索…'); $('runLiteratureSearch').disabled=true;
    try { const result=await api('/api/literature',{query,limit:Number($('literatureLimit').value||500)}); state.papers=result.results||[];state.searchMeta=result;state.paperPage=0; state.decisions.clear(); renderPapers(); $('exportLiterature').disabled=!state.papers.length; setStatus('literatureStatus',`PubMed 命中 ${result.total??'未知'} 条，本次取回 ${state.papers.length} 条；${result.complete?'本次分页取回完整。':'当前取回受上限或接口限制，报告中保留完整性标记。'}列表每页显示 50 条，导出包含全部取回记录；跨库去重和全文筛选尚未在本页自动完成。`); }
    catch(error){setStatus('literatureStatus',`检索失败：${error.message}`,true);} finally{$('runLiteratureSearch').disabled=false;}
  }
  function exportPapers() {
    const headers=['pmid','title','year','abstract','url','screening_decision','screening_reason'];
    const lines=[headers.map(csvEscape).join(',')];
    state.papers.forEach((p,i)=>{const key=String(p.pmid||i);lines.push([p.pmid,p.title,p.year,p.abstract,p.url,state.decisions.get(key)||'待核实',state.decisions.get(key+':reason')||''].map(csvEscape).join(','));});
    download('literature-screening-asreview.csv',lines.join('\n'));
  }
  const metaTemplate=`study,yi,sei,subgroup,risk_of_bias\nStudy A,0.32,0.12,model_1,low\nStudy B,0.18,0.10,model_2,high\n`;
  function downloadMetaTemplate(){const kind=$('metaType').value;const csv=kind==='continuous'?'study,n_t,mean_t,sd_t,n_c,mean_c,sd_c,subgroup,risk_of_bias\nStudy A,48,5.2,1.4,45,4.7,1.2,model_1,low\nStudy B,32,6.1,1.8,30,5.4,1.5,model_2,high\n':kind==='binary'?'study,event_t,total_t,event_c,total_c,subgroup,risk_of_bias\nStudy A,12,80,18,82,model_1,low\nStudy B,9,64,16,61,model_2,high\n':metaTemplate;download('meta-analysis-template.csv',csv);}
  function table(headers,rows){return `<div class="research-table"><table><thead><tr>${headers.map(x=>`<th>${esc(x)}</th>`).join('')}</tr></thead><tbody>${rows.map(r=>`<tr>${r.map(x=>`<td>${esc(x)}</td>`).join('')}</tr>`).join('')}</tbody></table></div>`;}
  function renderMeta(data){
    state.lastMeta=data;
    const f=n=>Number(n).toFixed(4), o=data.overall;
    const eg=data.egger?.available?`Egger 截距 ${f(data.egger.intercept)}，正态近似 p=${f(data.egger.p_normal_approx)}。${esc(data.egger.caution)}`:esc(data.egger?.reason||'未计算发表偏倚检验。');
    let html=`<div class="meta-summary"><b>${esc(o.k)} 项研究</b><span>效应量：${esc(data.effect_type)}</span><span>随机效应合并：${f(o.estimate)} [${f(o.ci95[0])}, ${f(o.ci95[1])}]</span><span>I²=${f(o.i2_percent)}% · τ²=${f(o.tau2_dl)} · Q=${f(o.q)} (df=${o.q_df})</span><span>预测区间（正态近似）：[${f(o.prediction_interval_normal_approx[0])}, ${f(o.prediction_interval_normal_approx[1])}]</span></div>`;
    if(data.subgroups?.length) html+='<h3>预设亚组结果</h3>'+table(['亚组','研究数','合并效应','95% CI','I²'],data.subgroups.map(x=>[x.name,x.result.k,f(x.result.estimate),`${f(x.result.ci95[0])} – ${f(x.result.ci95[1])}`,`${f(x.result.i2_percent)}%`]))+`<p>亚组间 Q=${f(data.subgroup_test.q_between)}，df=${data.subgroup_test.df}，χ²近似 p=${f(data.subgroup_test.p_chi_square_approx)}。该检验是探索性结果，不能证明分组因果作用。</p>`;
    html+='<h3>逐项剔除敏感性分析</h3>'+table(['剔除研究','合并效应','I²'],data.leave_one_out.map(x=>[x.omitted_study,f(x.estimate),`${f(x.i2_percent)}%`]));
    if(data.high_risk_exclusion_sensitivity) html+='<h3>排除高偏倚风险研究</h3><p>敏感性合并效应 '+f(data.high_risk_exclusion_sensitivity.estimate)+'；I² '+f(data.high_risk_exclusion_sensitivity.i2_percent)+'%。</p>';
    const maxSe=Math.max(...data.funnel.map(x=>x.se),0.001);html+=`<h3>小研究效应/发表偏倚</h3><p>${eg}</p><p>${esc(data.heterogeneity_guidance)}</p><h3>漏斗图（效应量 × 标准误）</h3><div class="effect-points">${data.funnel.map(x=>{const px=Math.max(3,Math.min(97,50+((x.effect-o.estimate)/(1.96*maxSe))*45)),py=Math.max(5,Math.min(92,(x.se/maxSe)*88+4));return `<span title="${esc(x.study)}: 效应 ${f(x.effect)}，SE ${f(x.se)}" style="--point:${px}%;--depth:${py}%"></span>`}).join('')}</div><p class="hint">图形只作探索性显示；漏斗图不对称可能来自异质性、小样本效应或方法差异，不能单独确证发表偏倚。</p><h3>方法与限制</h3><ul>${data.warnings.map(x=>`<li>${esc(x)}</li>`).join('')}</ul><button class="button button-secondary" id="downloadMetaResult">下载完整 JSON 结果</button>`;
    $('metaResult').innerHTML=html; $('downloadMetaResult')?.addEventListener('click',()=>download('meta-analysis-result.json',JSON.stringify(data,null,2),'application/json;charset=utf-8'));
  }
  async function runMeta(){
    setStatus('metaStatus','正在读取并核对 CSV…'); $('runMeta').disabled=true;
    try {const csv=await readCsvInput('metaFile');const response=await api('/api/meta/analyze',{csv,outcome_type:$('metaType').value,group_column:$('metaGroup').value.trim(),risk_column:$('metaRisk').value.trim()});renderMeta(response.result);setStatus('metaStatus','分析已完成；请核对效应量定义、研究独立性和统计假设。');}
    catch(error){setStatus('metaStatus',`分析失败：${error.message}`,true);}finally{$('runMeta').disabled=false;}
  }
  async function runDiagnostics(){
    setStatus('diagnosticStatus','正在读取数据并检查格式、缺失值和样本结构…');$('runDiagnostics').disabled=true;
    try{const csv=await readCsvInput('diagnosticFile');const response=await api('/api/diagnostics/analyze',{csv,domain:$('diagnosticDomain').value,prediction:$('diagnosticPrediction').value.trim()});const data=response.result;state.lastDiagnostic=data;
      const summaries=Object.entries(data.numeric_summary).map(([name,x])=>[name,x.n,x.missing,Number(x.mean).toFixed(4),x.sd==null?'—':Number(x.sd).toFixed(4),`${x.min} – ${x.max}`]);
      $('diagnosticResult').innerHTML=`<div class="meta-summary"><b>${esc(data.domain)}</b><span>${data.row_count} 行</span><span>数值字段 ${summaries.length} 个</span></div><h3>字段概况</h3>${table(['字段','有效值','缺失','均值','SD','范围'],summaries)}<h3>需要核查的风险</h3><ul>${data.warnings.map(x=>`<li>${esc(x)}</li>`).join('')}</ul><h3>建议的改进顺序</h3><ol>${data.next_steps.map(x=>`<li>${esc(x)}</li>`).join('')}</ol><div class="research-warning">${esc(data.model_update)}</div><h3>待核实的科学约束草案</h3><pre class="constraint-preview">${esc(JSON.stringify({type:'unvalidated_constraint_draft',domain:data.domain,source_file:$('diagnosticFile').files[0]?.name||'',prediction:$('diagnosticPrediction').value.trim(),conditions:['species','tissue_or_cell','disease_stage','dose','time','control','assay'],rules:[],counterevidence:[],validation_plan:['independent_dataset_or_experiment','compare_to_previous_model','report_uncertainty'],status:'requires_human_review'},null,2))}</pre>`;
      $('downloadConstraints').disabled=false;setStatus('diagnosticStatus','检查完成：这是数据与流程审查，未执行算法训练。');
      $('draftDiagnosticSearch').disabled=false;
    }catch(error){setStatus('diagnosticStatus',`检查失败：${error.message}`,true);}finally{$('runDiagnostics').disabled=false;}
  }
  async function draftDiagnosticSearch(){
    if(!state.lastDiagnostic)return;setStatus('diagnosticEvidenceStatus','正在根据数据摘要拟定检索式…');$('draftDiagnosticSearch').disabled=true;
    try{const d=state.lastDiagnostic;const prompt=`请根据以下去标识化实验数据的结构摘要和研究者描述，为验证结果不一致的原因设计高召回 PubMed 检索。只使用提供的研究对象、实验条件和数据类型；未知信息列待确认。不要虚构论文或声称已搜索。请在 JSON 的 pubmed_query 字段提供可执行 PubMed 检索式，在 reply 中简要说明检索概念、需要人工核对的条件和主要偏倚风险。\n模块：${d.domain}\n研究者描述：${$('diagnosticPrediction').value.trim()||'未提供'}\n数据摘要：${JSON.stringify({columns:d.columns,row_count:d.row_count,numeric_summary:d.numeric_summary,warnings:d.warnings})}`;const response=await api('/api/chat',{messages:[{role:'user',content:prompt}],literature:[]});$('diagnosticEvidencePanel').hidden=false;$('diagnosticPubmedQuery').value=response.pubmed_query||'';setStatus('diagnosticEvidenceStatus',`${response.reply||'检索方案已生成。'}\n请检查并编辑检索式后，再点击“检索相关文献”。`);$('searchDiagnosticEvidence').disabled=!(response.pubmed_query||'').trim();}
    catch(error){setStatus('diagnosticEvidenceStatus',`检索式生成失败：${error.message}`,true);}finally{$('draftDiagnosticSearch').disabled=false;}
  }
  async function searchDiagnosticEvidence(){
    const query=$('diagnosticPubmedQuery').value.trim();if(!query){setStatus('diagnosticEvidenceStatus','请先输入并核对 PubMed 检索式。',true);return;}$('searchDiagnosticEvidence').disabled=true;setStatus('diagnosticEvidenceStatus','正在实际检索 PubMed…');
    try{const result=await api('/api/literature',{query,limit:100});state.diagnosticPapers=result.results||[];const shown=state.diagnosticPapers.slice(0,20);$('diagnosticEvidenceResults').innerHTML=`<p>PubMed 命中 ${esc(result.total??'未知')} 条，本次取回 ${state.diagnosticPapers.length} 条；以下显示前 ${shown.length} 条。${result.complete?'':'本次结果不完整。'}</p>${shown.length?`<div class="research-table"><table><thead><tr><th>题名 / PMID</th><th>摘要</th><th>年份</th></tr></thead><tbody>${shown.map(p=>`<tr><td><a href="${esc(p.url)}" target="_blank" rel="noopener noreferrer">${esc(p.title)} ↗</a><small>PMID ${esc(p.pmid)}</small></td><td>${esc(p.abstract||'无摘要')}</td><td>${esc(p.year||'—')}</td></tr>`).join('')}</tbody></table></div>`:'<p>没有检索到记录。请调整检索式。</p>'}<button class="button button-secondary" id="exportDiagnosticEvidence">导出本次全部记录 CSV</button>`;setStatus('diagnosticEvidenceStatus','检索完成。摘要与题名仅用于初筛线索，建议继续在文献筛选页记录纳入/排除理由。');$('exportDiagnosticEvidence')?.addEventListener('click',exportDiagnosticPapers);}
    catch(error){setStatus('diagnosticEvidenceStatus',`PubMed 检索失败：${error.message}`,true);}finally{$('searchDiagnosticEvidence').disabled=false;}
  }
  function exportDiagnosticPapers(){const lines=[['pmid','title','year','abstract','url'].map(csvEscape).join(',')];state.diagnosticPapers.forEach(p=>lines.push([p.pmid,p.title,p.year,p.abstract,p.url].map(csvEscape).join(',')));download('experiment-review-literature.csv',lines.join('\n'));}
  function exportConstraints(){if(!state.lastDiagnostic)return;download('unvalidated-science-constraint-draft.json',JSON.stringify({type:'unvalidated_constraint_draft',diagnostic:state.lastDiagnostic,prediction:$('diagnosticPrediction').value.trim(),source_file:$('diagnosticFile').files[0]?.name||'',status:'requires_human_review_and_external_validation'},null,2),'application/json;charset=utf-8');}
  function init(){
    $('planLiterature')?.addEventListener('click',makeLiteraturePlan);$('runLiteratureSearch')?.addEventListener('click',runLiterature);$('exportLiterature')?.addEventListener('click',exportPapers);
    $('downloadMetaTemplate')?.addEventListener('click',downloadMetaTemplate);$('metaType')?.addEventListener('change',()=>{const labels={effect:'下载效应量模板',continuous:'下载连续结局模板',binary:'下载二分类模板'};$('downloadMetaTemplate').textContent='↓ '+labels[$('metaType').value]});$('runMeta')?.addEventListener('click',runMeta);
    $('runDiagnostics')?.addEventListener('click',runDiagnostics);$('downloadConstraints')?.addEventListener('click',exportConstraints);$('draftDiagnosticSearch')?.addEventListener('click',draftDiagnosticSearch);$('searchDiagnosticEvidence')?.addEventListener('click',searchDiagnosticEvidence);$('diagnosticPubmedQuery')?.addEventListener('input',()=>{$('searchDiagnosticEvidence').disabled=!$('diagnosticPubmedQuery').value.trim()});
    document.querySelectorAll('[data-page="literature"]').forEach(el=>el.addEventListener('click',()=>{}));
    const syncResearchMirror=()=>{const source=$('researchResult'),mirror=$('assistantResearchResult');if(source&&mirror&&source.innerHTML.trim()&&mirror.innerHTML!==source.innerHTML)mirror.innerHTML=source.innerHTML;const status=$('researchStatus')?.textContent;if(status)setStatus('assistantResearchStatus',status);};
    new MutationObserver(syncResearchMirror).observe(document.body,{childList:true,subtree:true,characterData:true});
    document.addEventListener('click',e=>{if(e.target?.id==='downloadResearch'&&e.target.closest('#assistantResearchResult'))$('researchResult')?.querySelector('#downloadResearch')?.click();});
  }
  document.addEventListener('DOMContentLoaded',init);
})();
