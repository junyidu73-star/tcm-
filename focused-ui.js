// Two workspaces; optional data review lives inside the assistant.
document.querySelectorAll('.main-nav [data-page]').forEach(el=>{
  if(!['assistant','discovery'].includes(el.dataset.page))el.remove();
});
const assistantPage=document.getElementById('page-assistant');
const modes=document.createElement('div');
modes.className='panel research-tool-card';
modes.innerHTML='<label class="field-label" for="assistantFeedbackMode">这次希望怎样处理？</label><select class="text-input" id="assistantFeedbackMode"><option value="research">开展研究或继续当前任务</option><option value="requirements">结果不符合要求：按原要求重新检查</option><option value="experiment">有实验数据：复核实验与预测差异</option></select><p id="feedbackModeHelp" class="hint">直接在对话中输入疾病或研究要求。</p>';
assistantPage.querySelector('.assistant-mode-banner').after(modes);
const dataReview=document.createElement('details');
dataReview.className='panel research-tool-card';dataReview.hidden=true;
dataReview.innerHTML='<summary>实验数据复核（仅在有数据时使用）</summary>';
document.querySelectorAll('#page-diagnostics > .research-tool-card').forEach(panel=>dataReview.append(panel));
modes.after(dataReview);
const selector=document.getElementById('assistantFeedbackMode');
selector.addEventListener('change',()=>{
  const mode=selector.value;
  dataReview.hidden=mode!=='experiment';dataReview.open=mode==='experiment';
  document.getElementById('feedbackModeHelp').textContent=mode==='requirements'?'在对话中写原要求、当前结果和不符合之处；无需上传数据。助手会逐项核对、指出证据缺口并提出修订。':mode==='experiment'?'可展开下方上传 CSV 检查数据；也可先在对话中描述实验差异。':'直接在对话中输入疾病或研究要求。';
});
document.querySelector('.assistant-side .workflow-panel').hidden=true;
document.querySelector('.brand').href='#assistant';
dataReview.querySelectorAll('[data-page="literature"]').forEach(el=>el.remove());
runTask=async function(id){
  if(id==='literature'){await searchPapers(latestQuery||queryFor(activeBrief));return;}
  if(id==='omics'||id==='candidates'){navigate('discovery');return;}
  navigate('assistant');document.getElementById('chatInput').value='请根据当前研究要求和已有证据，设计能区分竞争机制的验证实验，并标明未知条件。';
};
document.querySelectorAll('#page-assistant .assistant-results-mirror').forEach(el=>{el.hidden=true;});
// Results are available on the gene-to-drug page; tools stay within the assistant.
navigate(location.hash.slice(1));
