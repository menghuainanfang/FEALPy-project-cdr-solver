const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const api=require('../experiments/time_evolution_view.js');
const root=path.resolve(__dirname,'../results/mentor_2026-10-23/p2_v2');

test('packed diagnostics preserve accepted states and missing coordinates remain explicit',()=>{
  const data=JSON.parse(fs.readFileSync(path.join(root,'burgers_n8_steps30_dt0.02/dataset.json'),'utf8'));
  const key=data.states[1].outer_iterations[0].mesh_inner;
  const rows=api.trace(key);
  assert.equal(rows.length,key.solver_stats.accepted_bdf_steps+1);
  assert.equal(rows[0].accepted_step_index,0);
  assert.equal(rows.at(-1).criterion_met,true);
  assert.equal(api.coordinates(key,1).accepted_step_index,1);
  const thin=data.states[2].outer_iterations[0].mesh_inner;
  assert.equal(thin.coordinate_coverage,'endpoints_only');
  assert.equal(api.coordinates(thin,1),null);
  assert.ok(api.coordinates(thin,0).physical);
  assert.ok(api.coordinates(thin,api.trace(thin).length-1).computational);
  assert.ok(api.displacement(data.states[1],data.config.n)>0);
});
