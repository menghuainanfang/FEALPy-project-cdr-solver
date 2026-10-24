/* Shared accepted-state explanations for frozen mesh epochs. */
(function () {
  "use strict";
  const dictionary = globalThis.MESH_FIELD_DICTIONARY;
  const definitions = Object.fromEntries(dictionary.fields.map(f => [f.key, f]));
  const number = x => x === null || x === undefined ? "未提供" : Number(x).toExponential(3);
  const count = x => x === null || x === undefined ? "未提供" : String(x);

  function normalizeScene(scene) {
    const modern = Array.isArray(scene.convergence_trace);
    const accepted = scene.accepted_step_trace;
    const raw = modern ? scene.convergence_trace : accepted || scene.frames || [];
    const trace = raw.map((r, i) => modern ? r : ({
      mesh_pseudo_time: r.pseudo_time,
      mesh_energy: r.energy,
      mesh_residual_rms: r.scaled_residual_rms ?? null,
      accepted_step_index: accepted ? r.accepted_step_index : null,
      pseudo_step_size: accepted ? r.accepted_dt : null,
      residual_threshold: scene.residual_threshold ?? null,
      nfev: r.nfev ?? null, njev: r.njev ?? null, nlu: r.nlu ?? null,
      epoch_id: "历史冻结问题", criterion_met: null, geometry_valid: null,
      consecutive_success: null, required_successes: null,
      unmet_reasons: ["旧数据未提供逐步监控判据；不补写当前收敛结论。"],
      stop_reason: i === raw.length-1 ? scene.stop_reason : null
    }));
    return {scene, modern, trace, savedFrames: scene.saved_frames ?? scene.frames?.length ?? null,
      acceptedSteps: scene.solver_stats?.accepted_bdf_steps ?? null};
  }

  function acceptedAtOrBefore(trace, s) {
    let index = 0;
    for (let i=1; i<trace.length && trace[i].mesh_pseudo_time<=s; i++) index=i;
    return index;
  }

  function counterRows(model, row, frameIndex) {
    const j = row.accepted_step_index;
    return [
      ["saved_frames", model.savedFrames === null ? "未提供" : frameIndex == null ? `输出采样共${model.savedFrames}帧；此处按BDF接受状态查看` : `可视化 第${frameIndex+1}/${model.savedFrames}帧`],
      ["accepted_bdf_steps", j === null || model.acceptedSteps === null ? "未提供" : j===0 ? `BDF初始状态 j=0 / 总接受步${model.acceptedSteps}` : `BDF接受步 第${j}/${model.acceptedSteps}步`],
      ["nfev", `当前检查 ${count(row.nfev)}；终值 ${count(model.scene.solver_stats?.nfev)}`],
      ["njev", `当前检查 ${count(row.njev)}；终值 ${count(model.scene.solver_stats?.njev)}`],
      ["nlu", `当前检查 ${count(row.nlu)}；终值 ${count(model.scene.solver_stats?.nlu)}`],
      ["rejected_bdf_attempts", model.scene.solver_stats?.rejected_bdf_attempts == null ? "未提供（SciPy公开API未暴露）" : model.scene.solver_stats.rejected_bdf_attempts],
      ["pseudo_step_size", number(row.pseudo_step_size)],
      ["wall_clock_seconds", number(model.scene.wall_clock_seconds ?? model.scene.wall_seconds)]
    ];
  }

  function guideIndices(model) {
    const rows=model.trace, last=rows.length-1;
    const first=predicate=>{const i=rows.findIndex(predicate); return i<0?last:i;};
    return [0, Math.floor(last*.4), first(r=>r.conditions_met===true),
      first(r=>r.criterion_met===true), last];
  }

  function mount(root, {onSelectAccepted=()=>{}}={}) {
    if (!document.getElementById("mesh-diag-style")) {
      const css=document.createElement("style"); css.id="mesh-diag-style";
      css.textContent=`.mesh-diag{margin:20px 0;padding:20px;border:1px solid #cbdad5;border-radius:12px;background:white;color:#17343c}.mesh-diag h2{font-size:22px;margin:0 0 10px}.mesh-diag h3{font-size:17px}.mesh-diag .checks{display:grid;grid-template-columns:repeat(3,1fr);gap:10px}.mesh-diag .check{border:1px solid #d5e0df;padding:12px;border-radius:8px;overflow-wrap:anywhere}.mesh-diag .check strong{display:block}.mesh-diag .help{font:inherit;border:1px solid #a6bfba;border-radius:50%;padding:0 6px;margin-left:5px;background:white;color:#356775;cursor:help}.mesh-diag .charts{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}.mesh-diag canvas{width:100%;height:210px;aspect-ratio:auto;cursor:crosshair}.mesh-diag table{width:100%;border-collapse:collapse;font-size:14px}.mesh-diag th,.mesh-diag td{text-align:left;border-bottom:1px solid #d5e0df;padding:8px;overflow-wrap:anywhere}.mesh-diag .guide{display:flex;gap:8px;flex-wrap:wrap}.mesh-diag .guide button{font:inherit;background:#edf5f2;border:1px solid #b8cec8;padding:8px;cursor:pointer;border-radius:5px}.mesh-diag button:focus-visible{outline:3px solid #d28c32;outline-offset:3px}.mesh-diag .position{font-variant-numeric:tabular-nums}.mesh-diag .outcome{padding:12px;background:#edf5f2;white-space:pre-wrap}.mesh-diag .subtle{font-size:13px;color:#58717b}@media(max-width:850px){.mesh-diag .checks,.mesh-diag .charts{grid-template-columns:1fr}.mesh-diag{padding:14px}}`;
      document.head.append(css);
    }
    root.classList.add("mesh-diag");
    root.innerHTML=`<h2>这个网格算收敛了吗？</h2><p class="position"></p><div class="checks"></div><p class="outcome"></p><div class="charts"><section><h3>r_RMS 对网格伪时间 s</h3><canvas data-chart="residual" role="img" aria-label="网格残量与阈值曲线，可点击定位"></canvas></section><section><h3>冻结度量能量 I_h 对 s</h3><canvas data-chart="energy" role="img" aria-label="单个冻结epoch的能量曲线，可点击定位"></canvas></section><section><h3>Δs_BDF 对接受步 j</h3><canvas data-chart="step" role="img" aria-label="BDF真实接受步长曲线，可点击定位"></canvas></section></div><p class="subtle">蓝线为当前轨迹；残量图虚线为阈值，竖线为首次完整收敛检查。点击曲线可选择实际接受状态。每张能量曲线仅属于当前冻结epoch。</p><h3>数字是什么意思？</h3><table><thead><tr><th>字段</th><th>当前检查 / 本次总计</th></tr></thead><tbody></tbody></table><p class="subtle">帧数、接受步数与RHS次数没有固定大小顺序。网格伪时间不是墙钟秒；网格收敛与积分器成功、PDE状态分别判断。</p><h3>3分钟看懂收敛</h3><div class="guide"></div><p class="guide-text"></p><p class="subtle">满足这些网格平衡判据不证明局部或全局能量极小值。</p>`;
    const find=s=>root.querySelector(s);
    let model=null, position=null, selected=0;
    function help(key) {
      const b=document.createElement("button"); b.className="help"; b.type="button"; b.textContent="?";
      b.title=definitions[key].definition; b.setAttribute("aria-label",definitions[key].label+"："+b.title); return b;
    }
    function check(label, value, description, key) {
      const div=document.createElement("div"); div.className="check";
      const title=document.createElement("span"); title.textContent=label;
      if(key) title.append(help(key)); else {const b=document.createElement("button"); b.className="help"; b.type="button"; b.textContent="?"; b.title=description; b.setAttribute("aria-label",description); title.append(b);}
      const v=document.createElement("strong"); v.textContent=value; div.append(title,v); return div;
    }
    function plot(canvas, kind) {
      const rect=canvas.getBoundingClientRect(), ratio=globalThis.devicePixelRatio||1;
      canvas.width=Math.max(1,rect.width*ratio); canvas.height=210*ratio;
      const ctx=canvas.getContext("2d"); ctx.scale(ratio,ratio);
      const rows=model.trace, w=rect.width, h=210, pad={l:60,r:15,t:20,b:30};
      const xValue=r=>kind==="step"?r.accepted_step_index:r.mesh_pseudo_time;
      const yValue=r=>kind==="step"?r.pseudo_step_size:kind==="energy"?r.mesh_energy:r.mesh_residual_rms;
      const data=rows.filter(r=>xValue(r)!=null&&yValue(r)!=null&&(kind!=="step"||yValue(r)>0));
      if(!data.length){ctx.fillText("本记录未提供此曲线",15,50);return;}
      const log=kind!=="energy", transform=v=>log?Math.log10(Math.max(v,1e-16)):v;
      const threshold=kind==="residual"?rows[0].residual_threshold:null;
      const values=data.map(r=>transform(yValue(r))); if(threshold!=null)values.push(transform(threshold));
      let lo=Math.min(...values),hi=Math.max(...values); if(hi-lo<1e-12){lo-=.5;hi+=.5;}
      const xmax=Math.max(...data.map(xValue),1e-12), xx=v=>pad.l+(w-pad.l-pad.r)*v/xmax, yy=v=>pad.t+(h-pad.t-pad.b)*(hi-transform(v))/(hi-lo);
      ctx.font="11px system-ui";ctx.fillStyle="#58717b";
      for(let i=0;i<3;i++){const transformed=lo+(hi-lo)*i/2, y=pad.t+(h-pad.t-pad.b)*(hi-transformed)/(hi-lo);ctx.fillText(log?"10^"+transformed.toFixed(1):transformed.toPrecision(3),2,y+4);ctx.strokeStyle="#e3ece7";ctx.beginPath();ctx.moveTo(pad.l,y);ctx.lineTo(w-pad.r,y);ctx.stroke();}
      ctx.fillText("0",pad.l,h-8);ctx.fillText(xmax.toPrecision(4),Math.max(pad.l,w-65),h-8);
      ctx.strokeStyle="#087f82";ctx.lineWidth=1.7;ctx.beginPath(); let epoch=null;
      data.forEach(r=>{if(epoch!==r.epoch_id){ctx.moveTo(xx(xValue(r)),yy(yValue(r)));epoch=r.epoch_id;}else ctx.lineTo(xx(xValue(r)),yy(yValue(r)));});ctx.stroke();
      if(threshold!=null){ctx.strokeStyle="#b07a22";ctx.setLineDash([5,4]);ctx.beginPath();ctx.moveTo(pad.l,yy(threshold));ctx.lineTo(w-pad.r,yy(threshold));ctx.stroke();ctx.setLineDash([]);}
      const first=rows.find(r=>r.criterion_met===true);if(first){ctx.strokeStyle="#b07a22";ctx.setLineDash([3,3]);ctx.beginPath();ctx.moveTo(xx(xValue(first)),pad.t);ctx.lineTo(xx(xValue(first)),h-pad.b);ctx.stroke();ctx.setLineDash([]);}
      const current=rows[selected]; if(xValue(current)!=null){ctx.strokeStyle="#647d97";ctx.beginPath();ctx.moveTo(xx(xValue(current)),pad.t);ctx.lineTo(xx(xValue(current)),h-pad.b);ctx.stroke();}
      canvas.onclick=event=>{const target=Math.max(0,Math.min(1,(event.clientX-rect.left-pad.l)/(w-pad.l-pad.r)))*xmax;
        let best=0;rows.forEach((r,i)=>{if(xValue(r)!=null&&Math.abs(xValue(r)-target)<Math.abs(xValue(rows[best])-target))best=i;});select(best);};
    }
    function render() {
      const r=model.trace[selected];
      const physical=position.physical_time==null?"":`物理时间 t=${position.physical_time.toFixed(3)}；`;
      find(".position").textContent=`${physical}可视化网格伪时间 s=${position.pseudo_time.toFixed(6)}；判据检查取实际接受状态 s=${r.mesh_pseudo_time.toFixed(6)}；epoch=${r.epoch_id}。`;
      const geom=r.geometry_valid==null?"未提供":r.geometry_valid?"是":"否";
      find(".checks").replaceChildren(
        check("主判据 · 网格残量",`${number(r.mesh_residual_rms)} / ${number(r.residual_threshold)}`,"", "mesh_residual_rms"),
        check("辅助 · 窗口相对位移",`${number(r.displacement)} / ${number(r.displacement_threshold)}`,`max |Δxi|/h_ref；固定窗口Δs=${number(r.window_pseudo_time)}，窗口就绪=${r.window_ready===true?"是":"否"}。`),
        check("辅助 · 冻结能量变化率",`${number(r.energy_change_rate)} / ${number(r.energy_rate_threshold)}`,`|ΔI|/max(|I_window|,I_scale)/(window/tau)；本次参考能量规模=${number(r.energy_reference)}。`),
        check("硬约束 · 几何合法",geom,"检查两个网格的最小正有向面积、边界法向与角点；固定边界另检查切向位移。"),
        check("连续达标次数",`${r.consecutive_success??"未提供"} / ${r.required_successes??"未提供"}`,"仅在实际接受步检查；所有单次条件满足才累加，否则归零。"),
        check("当前接受检查",r.criterion_met===true?"网格判据已满足":r.criterion_met===false?"尚未满足网格判据":"旧记录：逐步判据未提供","数值成功不等于网格收敛；最终停止与当前检查分别显示。")
      );
      const terminal=model.trace.at(-1), scene=model.scene;
      find(".outcome").textContent=`当前判据说明：${r.unmet_reasons?.length?r.unmet_reasons.join("；"):"无未达条件，当前完整判据已满足"}\n${model.modern?"本次最终停止":"历史输出的最终停止"}：${scene.stop_reason}；积分器=${scene.solver_status??"历史记录未分层"}；网格=${scene.mesh_status??"历史记录未分层"}；criterion_met=${scene.criterion_met??"未提供"}；cap_reached=${scene.cap_reached??"未提供"}。\n首次完整收敛检查 s=${terminal.first_converged_s==null?(model.modern?"本次未达到":"旧记录未提供"):terminal.first_converged_s.toFixed(6)}。`;
      const body=find("tbody");body.replaceChildren();
      for(const [key,value] of counterRows(model,r,position.display_frame_index)){const tr=document.createElement("tr"),name=document.createElement("th"),v=document.createElement("td");name.textContent=definitions[key].label;name.append(help(key));v.textContent=value;tr.append(name,v);body.append(tr);}
      for(const canvas of root.querySelectorAll("canvas"))plot(canvas,canvas.dataset.chart);
    }
    function select(i) {const target=model.trace[i]; onSelectAccepted(i,target);selected=i;render();}
    const stages=["1 初始残量","2 逐步下降","3 首次单次全条件","4 连续达标","5 最终停止"];
    stages.forEach((label,i)=>{const b=document.createElement("button");b.type="button";b.textContent=label;b.onclick=()=>{const index=guideIndices(model)[i];select(index);const r=model.trace[index];find(".guide-text").textContent=`${label}：实际接受状态 s=${r.mesh_pseudo_time.toFixed(6)}，残量=${number(r.mesh_residual_rms)}；${r.criterion_met===true?"已达到本例完整判据":r.unmet_reasons.join("；")}。`;};find(".guide").append(b);});
    return {update(scene,context){model=normalizeScene(scene);position=context;selected=acceptedAtOrBefore(model.trace,context.pseudo_time);render();}};
  }
  const api={normalizeScene,acceptedAtOrBefore,counterRows,guideIndices,mount};
  if(typeof module!=="undefined"&&module.exports)module.exports=api;
  else globalThis.MeshDiagnostics=api;
})();
