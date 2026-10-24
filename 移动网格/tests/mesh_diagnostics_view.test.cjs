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
  assert.equal(model.trace.at(-1).criterion_met,false);
  assert.match(api.counterRows(model,trace[2],6)[2][1],/终值 9$/);
});
