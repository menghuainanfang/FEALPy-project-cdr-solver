/* Complete accepted-state playback of nested mesh optimization. */
(function () {
  'use strict';
  function timeline(caseData) {
    return caseData.rounds.flatMap((round,k)=>{
      const last=round.inner.accepted_step_trace.length-1;
      return [{k,j:0,phase:'setup'},
        ...Array.from({length:last},(_,i)=>({k,j:i+1,phase:'inner'})),
        {k,j:last,phase:'mapping'},{k,j:last,phase:'check'}];
    });
  }
  function position(caseData, event) {
    const round=caseData.rounds[event.k];
    return {round, mesh:round.inner.accepted_step_trace[event.j],
      diagnostic:round.inner.convergence_trace[event.j],
      complete:event.phase==='check'&&round.outer_converged};
  }
  const api={timeline,position};
  if(typeof module!=='undefined'&&module.exports){module.exports=api;return;}
  const data=JSON.parse(document.getElementById('flow-data').textContent);
  const el=id=>document.getElementById(id);
  let events=[],at=0,timer=null,tracked=0;
  const currentCase=()=>data.cases[el('case').value];
  const scientific=x=>x===0?'0':x.toExponential(2);
  const labels={setup:'建立计算映射',inner:'优化计算网格',mapping:'映射回物理网格',check:'检查外层收敛'};
  function pause(){clearTimeout(timer);timer=null;el('play').textContent='播放全流程';}
  function mesh(id,nodes,ghost){
    const c=el(id),rect=c.getBoundingClientRect(),ratio=globalThis.devicePixelRatio||1;
    const w=rect.width,h=rect.width;c.width=w*ratio;c.height=h*ratio;
    const ctx=c.getContext('2d');ctx.scale(ratio,ratio);
    const pad=24,span=w-2*pad,point=p=>[pad+p[0]*span,h-pad-p[1]*span];
    const model=currentCase();
    function lines(points,color,width,dash=[]){ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);ctx.beginPath();
      for(const [a,b] of model.edges){ctx.moveTo(...point(points[a]));ctx.lineTo(...point(points[b]));}ctx.stroke();ctx.setLineDash([]);}
    if(ghost)lines(ghost,'#d6dfdc',.65,[3,4]);
    lines(nodes,id==='computational'?'#57718a':'#087b69',.85);
    for(let i=0;i<nodes.length;i++)if(model.boundary[i]){ctx.fillStyle='#bb7a27';ctx.beginPath();ctx.arc(...point(nodes[i]),2.5,0,2*Math.PI);ctx.fill();}
    const p=point(nodes[tracked]);ctx.fillStyle='#cb6438';ctx.beginPath();ctx.arc(...p,5,0,2*Math.PI);ctx.fill();
    c.onclick=event=>{const x=(event.clientX-rect.left-pad)/span,y=1-(event.clientY-rect.top-pad)/span;
      tracked=nodes.reduce((best,p,i)=>Math.hypot(p[0]-x,p[1]-y)<Math.hypot(nodes[best][0]-x,nodes[best][1]-y)?i:best,0);render();};
  }
  function render(){
    const model=currentCase(),event=events[at],view=position(model,event),r=view.round,d=view.diagnostic;
    const cfg=model.outer_convergence,last=r.inner.accepted_step_trace.length-1;
    el('timeline').value=at;el('previous').disabled=at===0;el('next').disabled=at===events.length-1;
    el('progress').textContent=`外迭代 ${event.k+1} / ${model.rounds.length} · ${labels[event.phase]} · 内层接受步 ${event.j} / ${last}`;
    el('result').textContent=view.complete?'外层收敛，优化结束':labels[event.phase];
    el('result').className=view.complete?'result complete':'result';
    for(const node of document.querySelectorAll('[data-phase]'))node.setAttribute('aria-current',String(node.dataset.phase===event.phase));
    el('rounds').value=event.k;
    const descriptions={
      setup:`以第 ${event.k+1} 轮物理网格建立映射；计算坐标重新从固定参考网格开始。本轮物理网格和度量保持冻结。`,
      inner:`只优化计算坐标 ξ，右侧同步展示映射得到的物理网格预览。内层${d.criterion_met?'已满足完整收敛判据':'继续优化'}。`,
      mapping:'内层优化完成。采用右侧映射结果作为本轮输出物理网格，并在新单元中心重新采样度量。',
      check:view.complete?'网格变化与度量变化均在阈值内，且已连续两轮达标。全流程在此结束。':
        '本轮尚未满足外层停止条件：'+r.unmet_reasons.join('；')+'。以当前物理网格进入下一轮。'};
    el('explanation').textContent=descriptions[event.phase];
    el('physical-caption').textContent=event.phase==='setup'?'本轮输入物理网格':event.phase==='inner'?'当前计算坐标映射得到的物理网格预览':'本轮输出物理网格；灰色虚线为本轮输入';
    el('computational-caption').textContent=`ξ(s) · 网格伪时间 s=${scientific(view.mesh.pseudo_time)}；橙色节点沿边移动，四角固定`;
    mesh('computational',view.mesh.computational,model.reference);
    mesh('physical',view.mesh.physical,r.input_physical);
    el('residual').textContent=`内层残量 ${scientific(d.mesh_residual_rms)} / 阈值 ${scientific(d.residual_threshold)}`;
    el('inner-check').textContent=d.criterion_met?'内层已收敛':`内层连续达标 ${d.consecutive_success} / ${d.required_successes}`;
    el('outer-check').textContent=event.phase==='check'?`外层网格变化 ${scientific(r.mesh_change)} · 度量变化 ${scientific(r.metric_change)} · 连续 ${r.consecutive_success} / ${cfg.consecutive_rounds}`:'外层判定在本轮映射完成后进行';
    const q=view.mesh.computational[tracked],x=view.mesh.physical[tracked];
    el('mapping-point').textContent=`同一顶点 ${tracked}：计算坐标 (${q.map(v=>v.toFixed(4)).join(', ')}) → 物理坐标 (${x.map(v=>v.toFixed(4)).join(', ')})。点击网格可更换跟踪点。`;
    const windowValue=d.window_ready?`位移 ${scientific(d.displacement)} ≤ ${scientific(d.displacement_threshold)}；能量变化率 ${scientific(d.energy_change_rate)} ≤ ${scientific(d.energy_rate_threshold)}`:'固定诊断窗口尚未就绪';
    el('criteria').textContent=`内层：残量低于阈值，固定窗口位移和能量变化率达标，满足最小伪时间及几何条件，连续 ${d.required_successes} 个接受步通过。当前${windowValue}。\n外层：max|Δx|/h_ref ≤ ${cfg.mesh_tolerance}，度量相对变化 ≤ ${cfg.metric_tolerance}，内层收敛且几何合法，连续 ${cfg.consecutive_rounds} 轮通过。度量采用固定 ${model.config.metric_relaxation} 松弛更新；变化判据使用松弛前差值。`;
  }
  function select(index){pause();at=Math.max(0,Math.min(index,events.length-1));render();}
  function choose(){pause();at=0;const model=currentCase();events=timeline(model);el('timeline').max=events.length-1;
    tracked=model.boundary.findIndex((b,i)=>b&&model.reference[i][0]>0&&model.reference[i][0]<1);
    el('rounds').replaceChildren();model.rounds.forEach((r,k)=>{const o=document.createElement('option');o.value=k;
      o.textContent=`第 ${k+1} 轮 · ${r.inner.solver_stats.accepted_bdf_steps} 个内层接受步`;el('rounds').append(o);});
    el('case-summary').textContent=`已记录至外层收敛 · ${model.rounds.length} 轮外迭代 · ${model.rounds.reduce((sum,r)=>sum+r.inner.solver_stats.accepted_bdf_steps,0)} 个内层接受步`;
    render();}
  function schedule(){const event=events[at],delay=event.phase==='inner'?Number(el('speed').value):850;
    timer=setTimeout(()=>{if(at===events.length-1){pause();return;}at++;render();schedule();},delay);}
  el('play').onclick=()=>{if(timer!==null){pause();return;}if(at===events.length-1){at=0;render();}el('play').textContent='暂停';schedule();};
  el('previous').onclick=()=>select(at-1);el('next').onclick=()=>select(at+1);el('reset').onclick=()=>select(0);
  el('final').onclick=()=>select(events.length-1);el('timeline').oninput=()=>select(Number(el('timeline').value));
  el('rounds').onchange=()=>select(events.findIndex(e=>e.k===Number(el('rounds').value)));
  for(const [key,model] of Object.entries(data.cases)){const o=document.createElement('option');o.value=key;o.textContent=model.label;el('case').append(o);}
  el('case').value='band';el('case').onchange=choose;
  window.addEventListener('resize',render);document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});
  choose();
})();
