// Optional science controls remain inside the single assistant workspace.
(()=>{
const panel=document.createElement('details');panel.className='panel research-tool-card';
panel.innerHTML=`<summary>科学约束、纠偏版本与模型训练</summary>
<p>把已核对的证据转为条件规则。规则影响后续任务的基因纳入与药物排序；AI提出的假设需要核对后才能采用。</p>
<label class="field-label" for="scienceRules">规则 JSON（所有适用条件须明确；any代表明确不限条件）</label>
<textarea class="text-input text-area" id="scienceRules" rows="10" spellcheck="false"></textarea>
<label class="field-label" for="scienceReason">本次修订原因或反证</label><input class="text-input" id="scienceReason" placeholder="原要求偏差、实验结果、适用条件修正等">
<label><input type="checkbox" id="scienceConfirmed">我已核对原文、关系方向及适用条件，确认这些规则可用于后续研究筛选。</label>
<div class="research-tool-actions"><button class="button button-primary" id="scienceSave">保存规则版本</button><button class="button button-secondary" id="scienceRefresh">刷新版本</button></div>
<label class="field-label" for="scienceVersions">已保存版本（可选择旧版本回滚）</label><select class="text-input" id="scienceVersions"></select>
<div class="research-tool-actions"><button class="button button-primary" id="scienceActivate">后续任务采用所选版本</button><button class="button button-secondary" id="sciencePreview">用所选版本复核本会话最近结果</button><button class="button button-secondary" id="scienceDisable">停用自定义规则</button></div>
<p class="hint">预览保留旧报告，仅复核已获得的记录；条件未知不自动套用规则。修订版本继承所选父版本的标识，需在JSON中保留仍然适用的规则。</p>
<details><summary>用核验标签训练约束神经网络原型（可选）</summary><p>上传训练JSON：task为gene或drug，rows中每项包含label、source_group、source_id、adjudicated=true、features。至少40条及8个独立研究组。特征：direct_replication、genetic_functional_support、perturbation_rescue、direction_match、exposure_match、safety_risk，数值0—1，未知null。标签对应预先定义的实验验证结局，不能用模型预测充当标签。</p>
<input class="text-input" id="scienceTrainingFile" type="file" accept=".json,application/json"><button class="button button-secondary" id="scienceTrain">训练并与逻辑回归基线比较</button><p class="hint">产物为离线研究模型，尚须独立验证；不会自动部署到当前筛选。</p></details>
<pre id="scienceOutput" style="white-space:pre-wrap;overflow-wrap:anywhere" role="status">尚未加载规则版本。</pre>`;
document.getElementById('page-assistant').append(panel);
const $=id=>document.getElementById(id),output=x=>{$('scienceOutput').textContent=typeof x==='string'?x:JSON.stringify(x,null,2);};
const context=document.createElement('details');context.innerHTML='<summary>规则适用场景（可选；未知留空）</summary><label class="field-label">物种</label><input class="text-input" id="scienceContextSpecies" placeholder="例如 Homo sapiens"><label class="field-label">组织或细胞</label><input class="text-input" id="scienceContextTissue"><label class="field-label">病程阶段</label><input class="text-input" id="scienceContextStage">';document.getElementById('diseaseResearchForm').append(context);
window.scientificContext=()=>({species:$('scienceContextSpecies').value.trim(),tissue:$('scienceContextTissue').value.trim(),stage:$('scienceContextStage').value.trim()});
$('scienceRules').value=JSON.stringify([{kind:'target_direction',value:'inhibit',subject:'',disease:'',species:'',tissue:'',stage:'',source_id:'',quote:'',source_text:''}],null,2);
const call=async(action,body={})=>researchApi('/api/science/'+action,body);
async function refresh(){const value=await call('state');$('scienceVersions').replaceChildren();for(const row of value.versions){const op=document.createElement('option');op.value=row.id;op.textContent=`${row.id} · ${row.rule_count}条 · ${row.reason}`;$('scienceVersions').append(op);}output({active:value.active||'无自定义规则',versions:value.versions});}
function bind(id,fn){$(id).onclick=async()=>{$(id).disabled=true;try{await fn();}catch(e){output(e.message);}finally{$(id).disabled=false;}};}
bind('scienceRefresh',refresh);
bind('scienceSave',async()=>{const value=await call('save',{rules:JSON.parse($('scienceRules').value),confirmed:$('scienceConfirmed').checked,parent:$('scienceVersions').value||null,reason:$('scienceReason').value});await refresh();$('scienceVersions').value=value.id;output({saved:value.id,message:'版本已保存，选择采用后仅作用于后续任务。'});});
bind('scienceActivate',async()=>{if(!$('scienceVersions').value)throw Error('请先保存并选择版本。');output(await call('activate',{id:$('scienceVersions').value}));});
bind('sciencePreview',async()=>{output(await call('preview',{id:$('scienceVersions').value}));});
bind('scienceDisable',async()=>output(await call('activate',{id:null})));
bind('scienceTrain',async()=>{const file=$('scienceTrainingFile').files[0];if(!file)throw Error('请先选择已核验标签JSON。');if(file.size>500000)throw Error('训练文件须小于500KB。');output('正在训练并做固定研究组留出验证…');output(await call('train',JSON.parse(await file.text())));});
})();
