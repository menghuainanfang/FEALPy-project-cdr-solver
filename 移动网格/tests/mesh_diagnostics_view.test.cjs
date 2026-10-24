const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const root=path.resolve(__dirname,'..');
globalThis.MESH_FIELD_DICTIONARY=JSON.parse(fs.readFileSync(path.join(root,'src/mesh_field_dictionary.json'),'utf8'));
const api=require('../experiments/mesh_diagnostics_view.js');
const read=relative=>JSON.parse(fs.readFileSync(path.join(root,relative),'utf8'));

test('archived E velocity is not relabelled as a residual or accepted step count',()=>{
  const data=read('results/e_variational_mesh.json');
  const scene=data.scenes.smoothing,model=api.normalizeScene(scene);
  assert.equal(model.savedFrames,61);
  assert.equal(model.acceptedSteps,null);
  assert.equal(model.trace[0].mesh_residual_rms,null);
  assert.equal(model.trace[0].criterion_met,null);
  const rows=api.counterRows(model,model.trace.at(-1),60);
  assert.match(rows[0][1],/第61\/61帧/);
  assert.equal(rows[1][1],'未提供');
  assert.ok(rows.every(row=>!row[1].includes('undefined')));
});

test('archived H 111 accepted steps remain distinct from eight output samples',()=>{
  const scene=read('results/mentor_2026-10-23/burgers/h_burgers.json').states[1].outer_iterations[0].mesh_inner;
  const model=api.normalizeScene(scene),last=model.trace.at(-1);
  assert.equal(model.acceptedSteps,111);
  assert.equal(model.savedFrames,8);
  assert.equal(last.criterion_met,null);
  const rows=api.counterRows(model,last,null);
  assert.match(rows[0][1],/共8帧/);
  assert.match(rows[1][1],/第111\/111步/);
  assert.match(rows[5][1],/未提供/);
});

test('modern trace uses recorded facts and guide never creates a convergence event',()=>{
  const trace=[0,1,2].map(j=>({mesh_pseudo_time:j*.1,accepted_step_index:j,mesh_residual_rms:1/(j+1),
    mesh_energy:3-j,criterion_met:false,conditions_met:false,epoch_id:'actual',unmet_reasons:['残量未达标']}));
  const model=api.normalizeScene({convergence_trace:trace,frames:Array(7),solver_stats:{accepted_bdf_steps:2,nfev:9},stop_reason:'pseudo_time_cap'});
  assert.equal(api.acceptedAtOrBefore(trace,.15),1);
  assert.deepEqual(api.guideIndices(model),[0,0,2,2,2]);
  assert.match(api.guideLabels(model)[2],/未达到/);
  assert.match(api.guideLabels(model)[3],/未达到/);
  assert.equal(model.trace.at(-1).criterion_met,false);
  assert.match(api.counterRows(model,trace[2],6)[2][1],/终值 9$/);
});

test('complete library replays the exact accepted meshes and contains distinct real outcomes',()=>{
  const data=read('results/mentor_2026-10-23/p1_complete/e_variational_mesh.json');
  assert.equal(Object.keys(data.scenes).length,8);
  const numeric=['mesh_residual_rms','residual_threshold','displacement_threshold',
    'energy_rate_threshold','consecutive_success','required_successes','nfev','njev','nlu'];
  for(const scene of Object.values(data.scenes)){
    const model=api.normalizeScene(scene);
    assert.equal(model.modern,true);
    assert.equal(scene.boundary_mode,'tangent');
    assert.equal(scene.frames.length,scene.convergence_trace.length);
    assert.equal(scene.frames.length,scene.solver_stats.accepted_bdf_steps+1);
    scene.frames.forEach((f,j)=>{
      const r=model.trace[j];
      assert.equal(f.pseudo_time,r.mesh_pseudo_time);
      assert.equal(f.accepted_step_index,r.accepted_step_index);
      assert.equal(f.energy,r.mesh_energy);
      assert.equal(f.physical.length,data.reference.length);
      assert.equal(f.computational.length,data.reference.length);
      assert.equal(f.accepted_step_index,scene.accepted_step_trace[j].accepted_step_index);
      for(const [coordinates,areaKey] of [[f.physical,'physical_min_area'],[f.computational,'computational_min_area']]){
        const minArea=Math.min(...data.cells.map(([a,b,c])=>{
          const p=coordinates[a],q=coordinates[b],t=coordinates[c];
          return ((q[0]-p[0])*(t[1]-p[1])-(q[1]-p[1])*(t[0]-p[0]))/2;
        }));
        assert.ok(minArea>0);
        assert.ok(Math.abs(minArea-r[areaKey])<1e-14);
      }
      for(const key of numeric)assert.equal(typeof r[key],'number',`${scene.case_key}: ${key}`);
      assert.equal(r.geometry_valid,true);
      assert.equal(typeof f.interpolation_rmse,'number');
      if(r.window_ready){assert.equal(typeof r.displacement,'number');assert.equal(typeof r.energy_change_rate,'number');}
    });
    for(const i of api.guideIndices(model))assert.equal(scene.frames[i].pseudo_time,model.trace[i].mesh_pseudo_time);
  }
  for(const key of ['adaptive','smoothing']){
    assert.equal(data.scenes[key].stop_reason,'converged');
    assert.equal(data.scenes[key+'_short'].criterion_met,false);
    assert.equal(data.scenes[key+'_short'].cap_reached,true);
    assert.deepEqual(data.scenes[key].fixed_physical,data.scenes[key+'_short'].fixed_physical);
  }
  assert.equal(data.scenes.equilibrium.convergence_trace[0].criterion_met,false);
  assert.equal(data.scenes.equilibrium.criterion_met,true);
  assert.equal(data.scenes.smoothing_continue.criterion_met,true);
  assert.equal(data.scenes.smoothing_continue.stop_reason,'pseudo_time_cap');
  assert.ok(data.scenes.smoothing_continue.first_converged_s<data.scenes.smoothing_continue.terminal_pseudo_time);
  assert.equal(data.scenes.adaptive_budget.stop_reason,'max_accepted_steps');
  assert.equal(data.scenes.adaptive_budget.cap_reached,false);
});

test('all current E and Poisson pages expose only tangential moving cases',()=>{
  for(const file of ['results/e_variational_mesh.html',
    'results/mentor_2026-10-23/e_convergence/e_variational_mesh.html',
    'results/mentor_2026-10-23/p1_v2/e/e_variational_mesh.html',
    'results/mentor_2026-10-23/p1_complete/e_variational_mesh.html']){
    const html=fs.readFileSync(path.join(root,file),'utf8');
    const data=JSON.parse(html.match(/<script id="data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
    assert.ok(Object.values(data.scenes).every(s=>s.boundary_mode==='tangent'));
    assert.ok(!html.includes('value="adaptive_fixed"'));
  }
  const html=fs.readFileSync(path.join(root,'results/f_poisson_comparison.html'),'utf8');
  const data=JSON.parse(html.match(/<script id="data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
  assert.ok(data.rows.every(r=>r.method!=='adaptive_fixed'));
  assert.ok(Object.values(data.datasets).every(s=>!('adaptive_fixed' in s.variants)));
  assert.ok(!html.includes('自适应·固定边界'));
});
