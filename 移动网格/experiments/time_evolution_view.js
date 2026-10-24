/* Physical time, coupling rounds and accepted mesh states remain separate. */
(function(){
  'use strict';
  function trace(inner){
    if(inner.convergence_trace)return inner.convergence_trace;
    return inner.trace_rows.map(values=>({...inner.trace_constants,...Object.fromEntries(inner.trace_columns.map((key,i)=>[key,values[i]]))}));
  }
  function coordinates(inner,j){
    const saved=inner.mesh_states.find(r=>r.accepted_step_index===j);
    if(saved)return saved;
    if(j===0)return inner.endpoint_states[0];
    if(j===trace(inner).length-1)return inner.endpoint_states.at(-1);
    return null;
  }
  function displacement(state,n){return Math.max(...state.nodes.map((p,i)=>Math.hypot(p[0]-state.input_physical[i][0],p[1]-state.input_physical[i][1])))*n;}
  const api={trace,coordinates,displacement};
  if(typeof module!=='undefined'&&module.exports){module.exports=api;return;}
  const library=JSON.parse(document.getElementById('time-data').textContent),el=id=>document.getElementById(id);
  const terms=Object.fromEntries(library.field_dictionary.fields.map(f=>[f.key,f]));
  el('terms').textContent=['physical_time','mesh_pseudo_time','wall_clock_seconds','physical_dt','outer_residual','outer_status','pde_status','coordinate_coverage'].map(key=>terms[key].label+'：'+terms[key].definition).join('\n\n');
  for(const [id,key] of [['time','physical_time'],['outer','outer_iteration'],['inner','mesh_pseudo_time'],['outer-status','outer_status'],['pde-status','pde_status']])el(id).title=terms[key].definition;
  let ti=0,ki=0,ji=0,timer=null,states=[],tracked=0,detailed=null;
  const model=()=>library.models[el('scenario').value],state=()=>states[ti],round=()=>(detailed||state()).outer_iterations[ki];
  const fmt=v=>v==null?'不适用':Number(v).toExponential(3);
  const names={converged:'已收敛',not_checked:'尚未检查',not_applicable:'不适用',not_converged:'继续迭代',cap_reached:'到达上限',failed:'失败',initial_interpolation:'初值插值'};
  function pause(){clearInterval(timer);timer=null;el('play-time').textContent=model().scenario==='static'?'回放静态优化':'播放物理时间';}
  function badge(node,label,status){node.textContent=label+'：'+(names[status]||status);node.dataset.status=status;}
  function surface(id){const c=el(id),w=c.getBoundingClientRect().width,h=w,dpr=globalThis.devicePixelRatio||1;c.width=w*dpr;c.height=h*dpr;const ctx=c.getContext('2d');ctx.scale(dpr,dpr);return {c,ctx,w,h,p:p=>[20+p[0]*(w-40),h-20-p[1]*(w-40)]};}
  function draw(id,nodes,{values=null,cellsValues=null,ghost=null,arrows=false}={}){
    const {c,ctx,w,h,p}=surface(id),m=model();
    const valuesUsed=values||cellsValues;
    if(valuesUsed){const lo=values?m.display_range[0]:Math.min(...valuesUsed),hi=values?m.display_range[1]:Math.max(...valuesUsed);
      m.cells.forEach((cell,i)=>{const v=cellsValues?cellsValues[i]:cell.reduce((sum,j)=>sum+values[j],0)/3,z=(v-lo)/(hi-lo||1);ctx.fillStyle=`hsl(${220-200*z} 58% ${88-35*z}%)`;ctx.beginPath();cell.forEach((j,k)=>k?ctx.lineTo(...p(nodes[j])):ctx.moveTo(...p(nodes[j])));ctx.closePath();ctx.fill();});
      ctx.font='10px system-ui';ctx.fillStyle='#31564f';ctx.fillText('范围 '+lo.toFixed(3)+' — '+hi.toFixed(3),20,12);}
    function lines(points,color,dash=[]){ctx.strokeStyle=color;ctx.lineWidth=.65;ctx.setLineDash(dash);ctx.beginPath();m.edges.forEach(([a,b],i)=>{if(el('density').value==='sparse'&&i%2)return;ctx.moveTo(...p(points[a]));ctx.lineTo(...p(points[b]));});ctx.stroke();ctx.setLineDash([]);}
    if(ghost)lines(ghost,'#c3cfca',[3,3]);lines(nodes,'#237265');
    if(arrows&&ghost){const a=p(ghost[tracked]),b=p(nodes[tracked]);ctx.strokeStyle='#b16828';ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(...a);ctx.lineTo(...b);ctx.stroke();}
    const point=p(nodes[tracked]);ctx.fillStyle='#ae5c30';ctx.beginPath();ctx.arc(...point,3,0,Math.PI*2);ctx.fill();
    c.onclick=e=>{const rect=c.getBoundingClientRect(),x=(e.clientX-rect.left-20)/(w-40),y=1-(e.clientY-rect.top-20)/(h-40);tracked=nodes.reduce((best,p,i)=>Math.hypot(p[0]-x,p[1]-y)<Math.hypot(nodes[best][0]-x,nodes[best][1]-y)?i:best,0);renderStory();};
  }
  function blank(id){const {ctx,w}=surface(id);ctx.font='13px system-ui';ctx.fillStyle='#647670';ctx.fillText('此接受步未保存完整坐标',20,60);ctx.fillText('诊断数据完整；请选关键物理步展开网格',20,90);}
  function physicalPlot(){const c=el('motion-plot'),w=c.getBoundingClientRect().width,h=150,dpr=globalThis.devicePixelRatio||1;c.width=w*dpr;c.height=h*dpr;const ctx=c.getContext('2d');ctx.scale(dpr,dpr);
    const staticMode=model().scenario==='static',values=states.map(s=>displacement(s,model().config.n)),max=Math.max(...values,1e-12),end=staticMode?states.length-1:states.at(-1).physical_time,x=t=>50+(w-70)*t/(end||1),y=v=>15+105*(1-v/max),xx=(s,i)=>staticMode?i:s.physical_time;
    ctx.font='10px system-ui';ctx.fillStyle='#647670';ctx.fillText(max.toExponential(2),0,20);ctx.fillText('0',30,125);ctx.fillText(staticMode?'优化前（t=0）':'物理时间 t=0',50,145);ctx.fillText(staticMode?'优化后（t=0）':'t='+end.toFixed(2),w-100,145);
    ctx.strokeStyle='#087b69';ctx.beginPath();states.forEach((s,i)=>i?ctx.lineTo(x(xx(s,i)),y(values[i])):ctx.moveTo(x(xx(s,i)),y(values[i])));ctx.stroke();ctx.fillStyle='#a76329';ctx.beginPath();ctx.arc(x(xx(state(),ti)),y(values[ti]),4,0,Math.PI*2);ctx.fill();
    c.onclick=e=>{const target=(e.clientX-c.getBoundingClientRect().left-50)/(w-70)*end;let best=0;states.forEach((s,i)=>{if(Math.abs(xx(s,i)-target)<Math.abs(xx(states[best],best)-target))best=i;});chooseTime(best);};}
  const diagnostics=MeshDiagnostics.mount(el('mesh-convergence-panel'),{onSelectAccepted:(j)=>{ji=j;el('inner').value=j;renderInner();}});
  async function loadDetailed(){
    if(detailed)return detailed;
    const selectedModel=model(),index=ti;
    const bytes=Uint8Array.from(atob(selectedModel.algorithm_blobs[index]),c=>c.charCodeAt(0));
    const text=await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).text();
    if(selectedModel!==model()||index!==ti)return null;
    detailed=JSON.parse(text);return detailed;
  }
  async function renderInner(){
    const s=state();let r=round();el('algorithm-content').hidden=!r;el('initial-note').hidden=!!r;
    if(!r)return;
    el('coverage').textContent='正在按需解压此物理步的实际算法记录…';
    if(!await loadDetailed())return;
    r=round();
    const inner=r.mesh_inner,rows=trace(inner);inner.convergence_trace=rows;const d=rows[ji],coords=coordinates(inner,ji);
    el('inner').max=rows.length-1;el('inner').value=ji;
    el('inner-label').textContent=`BDF接受步 j=${ji}/${rows.length-1} · 网格伪时间 s=${d.mesh_pseudo_time.toFixed(6)} · 残量 ${fmt(d.mesh_residual_rms)} / ${fmt(d.residual_threshold)}`;
    el('coupling-label').textContent=`物理步 n=${s.physical_step_index}，t=${s.physical_time.toFixed(3)}；选中 k=${r.outer_iteration}，本步实际 ${s.k_actual} 轮，预算最多 ${s.k_max} 轮；M版本=${r.metric_version}。`;
    el('outer-reason').textContent=`外层 ${names[r.outer_status]||r.outer_status}；残量 ${r.outer_residual==null?'未检查':fmt(r.outer_residual)} / ${r.outer_tol}。${r.outer_unmet_reasons.length?r.outer_unmet_reasons.join('；'):'本轮耦合条件已满足。'}${r.pde?' 选中轮PDE残量 '+fmt(r.pde.relative_nonlinear_residual)+' / 1e-9；状态 '+r.pde.status:' 解析目标案例，PDE残量不适用。'}`;
    el('coverage').textContent=inner.coordinate_coverage==='all_accepted'?'此物理步保留全部实际接受网格。':'此物理步保留初末网格与全部接受步诊断；中间坐标未保存，不作插值播放。关键物理步可查看完整网格。';
    if(coords){draw('inner-computational',coords.computational);draw('inner-physical',coords.physical,{ghost:r.frozen_physical});}else{blank('inner-computational');blank('inner-physical');}
    diagnostics.update(inner,{pseudo_time:d.mesh_pseudo_time,physical_time:s.physical_time,display_frame_index:null});
  }
  function renderStory(){
    const s=state(),m=model(),last=s.outer_iterations.at(-1),isStatic=m.scenario==='static';
    el('time').value=ti;el('previous-time').disabled=ti===0;el('next-time').disabled=ti===states.length-1;
    el('time-label').textContent=isStatic?`物理推进0步 · t=0 · ${ti?'静态网格优化后':'优化前'}`:`第 ${s.physical_step_index}/${states.at(-1).physical_step_index} 个物理时间步 · t=${s.physical_time.toFixed(3)} · Δt=${s.physical_dt.toFixed(3)}`;
    el('motion-caption').textContent=isStatic?'静态优化前后的网格变化；两个状态均为t=0，点击定位。':'实际每步网格重分布强度 max|Δx|/h；点击曲线定位物理时间节点。';
    const reference=el('reference').value==='fixed'?(m.scenario==='burgers'?m.fixed_grid_baseline[ti].nodes:m.reference):s.input_physical;
    draw('field',s.nodes,{values:s.values});draw('physical',s.nodes,{ghost:el('overlay').checked?reference:null,arrows:true});draw('metric',s.nodes,{cellsValues:s.metric.map(a=>a[0][0]+a[1][1])});
    el('field-caption').textContent=m.scenario==='burgers'?(ti?'本步最终数值 PDE 解 U_h':'制造解初值的有限元插值'):'解析目标在本步最终网格上的采样（未求解PDE）';
    el('mesh-caption').textContent=ti?'高亮：本步最终网格；灰虚线：'+(el('reference').value==='fixed'?'固定基线网格':'上一物理时刻')+'。点击选点，线段表示与参照的实际位移。':'初始物理网格。';
    badge(el('pde-status'),'本步最终PDE',m.scenario==='burgers'?(s.pde?.status||'initial_interpolation'):'not_applicable');
    badge(el('outer-status'),'本步最终外层',s.outer_status);
    badge(el('mesh-status'),'本步全部内层',last?(s.outer_iterations.every(r=>r.mesh_inner.criterion_met)?'converged':'cap_reached'):'not_applicable');
    const error=m.scenario==='burgers'?s.pde?.l2_error:s.target_interpolation_rmse;
    const base=m.scenario==='burgers'?m.fixed_grid_baseline?.[ti]?.pde?.l2_error:s.fixed_grid_target_rmse;
    el('numbers').textContent=`本步 max|Δx|/h=${fmt(displacement(s,m.config.n))} · 最小有向面积=${fmt(s.min_signed_area)} · ${m.scenario==='burgers'?'PDE L²误差':'目标插值 RMSE'}：移动网格 ${fmt(error)} / 固定网格 ${fmt(base)}${m.scenario==='burgers'?' · 本步最终PDE残量 '+fmt(s.pde?.relative_nonlinear_residual)+' / 1e-9':''}`;
    el('layer-reasons').textContent=ti?`物理层：${isStatic?'静态输入，物理时间未推进':'实际推进至 t='+s.physical_time.toFixed(3)}；外层：${names[s.outer_status]||s.outer_status}（实际 ${s.k_actual} 轮）；内层：${last?last.mesh_inner.stop_reason:'不适用'}。`:'初始状态，尚未进入物理时间与耦合循环。';
    physicalPlot();if(el('algorithm').open)renderInner();
  }
  function chooseTime(index,keepPlaying=false){if(!keepPlaying)pause();ti=Math.max(0,Math.min(index,states.length-1));detailed=null;ki=Math.max(0,state().outer_iterations.length-1);ji=round()?round().mesh_inner.solver_stats.accepted_bdf_steps:0;
    el('outer').replaceChildren();state().outer_iterations.forEach((r,k)=>{const o=document.createElement('option');o.value=k;o.textContent=`k=${r.outer_iteration} · M版本${r.metric_version}`;el('outer').append(o);});el('outer').value=ki;el('outer').disabled=!round();renderStory();}
  function rebuild(){pause();const m=model(),max=Math.max(...m.states.map(s=>s.physical_step_index));el('horizon').replaceChildren();
    for(const count of (max===0?[0]:[10,30,60,100].filter(n=>n<=max))){const o=document.createElement('option');o.value=count;o.textContent=count===0?'静态输入（0个物理步）':count+' 个物理步';el('horizon').append(o);}el('horizon').value=max>=60?60:max;setHorizon();
    const type=m.scenario==='burgers'?'已数值求解Burgers':m.scenario==='static'?'静态解析输入，不推进物理时间':'教学用时变解析目标网格演示（不求解非线性PDE）';el('type').textContent=type;
    el('title').textContent=m.scenario==='burgers'?'二维黏性 Burgers：解与网格一起演化':m.scenario==='static'?'静态输入：只优化网格':'移动界面：观察网格连续跟随';
    el('problem').textContent=m.scenario==='burgers'?'u_t + u(u_x+u_y) − 0.05Δu = f；P1有限元、后向Euler、Picard迭代。新时间步继承上一时刻网格。':'解析陡变界面沿水平方向周期移动，度量来自已知目标；每个时间节点重新计算目标并更新网格。';
    if(m.scenario==='static')el('problem').textContent='解析目标与物理时间均固定在t=0，仅记录网格优化前后状态。';
    [...el('guide').children].forEach((b,i)=>b.textContent=(i+1)+' '+(i===1&&m.scenario==='static'?'检查静态时间输入':descriptions[i]));
    el('sources').textContent=`网格 ${m.config.n}×${m.config.n}，${m.reference.length} 个P1自由度；数据版本 ${m.schema_version}。项目 ${m.provenance.project.sha}；FEALPy ${m.provenance.fealpy.sha}。完整 ${m.states.at(-1).physical_step_index} 个物理步数据集实测 ${m.performance.wall_seconds.toFixed(2)} 墙钟秒；选择较短回放不代表重新计时。\n复现：python -m experiments.three_loop_time_demo_v2 --scenario ${m.scenario} --physical-steps ${m.config.steps||10} --dt-physical ${m.config.physical_dt||.02} --nx ${m.config.n} --save-inner-traces key`;
  }
  function setHorizon(){const m=model(),count=Number(el('horizon').value);states=m.scenario==='static'?m.states:m.states.filter(s=>s.physical_step_index<=count);el('time').max=states.length-1;
    const candidates=states.map((s,i)=>s.outer_iterations.some(r=>r.mesh_inner.coordinate_coverage==='all_accepted')||i===0?i:null).filter(i=>i!==null);
    const chosen=[...new Set([0,.25,.5,.75,1].map(f=>candidates.reduce((best,i)=>Math.abs(i-f*(states.length-1))<Math.abs(best-f*(states.length-1))?i:best,candidates[0])))];
    el('bookmarks').replaceChildren();chosen.forEach(i=>{const b=document.createElement('button');b.textContent=m.scenario==='static'?(i===0?'静态输入 · 优化前':'静态输入 · 优化后'):`n=${states[i].physical_step_index} · t=${states[i].physical_time.toFixed(2)}`;b.onclick=()=>chooseTime(i);el('bookmarks').append(b);});
    const maxAt=states.reduce((best,s,i)=>displacement(s,m.config.n)>displacement(states[best],m.config.n)?i:best,0);el('strongest').textContent=`重分布最强：n=${states[maxAt].physical_step_index}`;el('strongest').onclick=()=>chooseTime(maxAt);chooseTime(0);}
  function playTime(){if(ti===states.length-1)chooseTime(0);el('play-time').textContent='暂停';timer=setInterval(()=>{chooseTime(ti+1,true);if(ti===states.length-1)pause();},Number(el('speed').value));}
  el('play-time').onclick=()=>{if(timer){pause();return;}playTime();};
  // chooseTime stops user-driven playback; the timer advances through the same selection path.
  el('time').oninput=()=>chooseTime(Number(el('time').value));el('previous-time').onclick=()=>chooseTime(ti-1);el('next-time').onclick=()=>chooseTime(ti+1);el('reset-time').onclick=()=>chooseTime(0);
  el('outer').onchange=()=>{pause();ki=Number(el('outer').value);ji=0;renderInner();};el('inner').oninput=()=>{pause();ji=Number(el('inner').value);renderInner();};el('algorithm').ontoggle=()=>{if(el('algorithm').open){pause();renderInner();}};
  el('speed').onchange=()=>{if(timer){pause();playTime();}};
  el('horizon').onchange=setHorizon;el('scenario').onchange=rebuild;el('overlay').onchange=renderStory;el('density').onchange=renderStory;el('reference').onchange=renderStory;
  const descriptions=['初始网格与场','物理时间推进，场发生变化','由新场建立度量 M','冻结 M，展开网格内迭代','重建网格并检查外层'];
    descriptions.forEach((label,g)=>{const b=document.createElement('button');b.textContent=(g+1)+' '+label;b.onclick=async()=>{
    chooseTime(g===0?0:1);let message='';const s=state(),r=s.outer_iterations[0];
    if(g===0)message=`初始 t=${s.physical_time}，${s.nodes.length} 个节点，物理推进0步。`;
    if(g===1)message=model().scenario==='static'?'物理时间保持t=0，解析目标未随时间变化；仅执行网格优化。':`实际 t=${s.physical_time.toFixed(3)}；${model().scenario==='burgers'?'数值PDE状态 '+s.pde.status:'解析目标更新，未求解PDE'}。`;
    if(g===2)message=`M版本 ${r.metric_version}，来自${model().scenario==='burgers'?'当前数值P1解梯度':'当前时刻解析目标'}；本轮冻结。`;
    if(g===3){el('algorithm').open=true;ki=0;ji=Math.floor(r.mesh_inner.solver_stats.accepted_bdf_steps/2);el('outer').value=0;await renderInner();const d=trace(round().mesh_inner)[ji];message=`n=${s.physical_step_index},k=1,j=${ji}；s=${d.mesh_pseudo_time.toFixed(4)}，残量 ${fmt(d.mesh_residual_rms)} / ${fmt(d.residual_threshold)}。`;}
    if(g===4)message=`本步外层 ${names[s.outer_status]}，实际 ${s.k_actual} 轮。${model().scenario==='burgers'?'非守恒P1转移引起质量变化 '+fmt(s.transfer_mass_change):'在新物理网格重新采样解析场与度量'}。`;
    el('guide-text').textContent=message;};el('guide').append(b);});
  for(const [key,m] of Object.entries(library.models)){m.display_range=m.states.flatMap(s=>s.values).reduce(([lo,hi],v)=>[Math.min(lo,v),Math.max(hi,v)],[Infinity,-Infinity]);const o=document.createElement('option');o.value=key;o.textContent=({analytic:'解析移动界面',burgers:'数值Burgers',static:'静态输入'})[m.scenario]+` · nx=${m.config.n} · Δt=${m.config.physical_dt}`;el('scenario').append(o);}
  const preferred=Object.keys(library.models).find(key=>key.startsWith('analytic_n8_dt0.02'));if(preferred)el('scenario').value=preferred;
  window.addEventListener('resize',renderStory);document.addEventListener('visibilitychange',()=>{if(document.hidden)pause();});rebuild();
})();
