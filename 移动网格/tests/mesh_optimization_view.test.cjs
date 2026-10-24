const {test}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const api=require('../experiments/mesh_optimization_view.js');
const root=path.resolve(__dirname,'..');
const data=JSON.parse(fs.readFileSync(path.join(root,'results/mentor_2026-10-23/optimization_flow/optimization_flow.json'),'utf8'));

test('playback visits every accepted state and every mapping before outer convergence',()=>{
  for(const model of Object.values(data.cases)){
    const events=api.timeline(model);
    assert.equal(events[0].phase,'setup');
    assert.equal(events.at(-1).phase,'check');
    assert.equal(api.position(model,events.at(-1)).complete,true);
    model.rounds.forEach((round,k)=>{
      const slice=events.filter(e=>e.k===k);
      const steps=round.inner.accepted_step_trace.length-1;
      assert.equal(slice.length,steps+3);
      assert.deepEqual(slice.slice(0,-2).map(e=>e.j),Array.from({length:steps+1},(_,i)=>i));
      assert.equal(slice.at(-2).phase,'mapping');
      assert.equal(slice.at(-1).phase,'check');
      assert.equal(slice.at(-1).j,steps);
      for(const e of slice){const p=api.position(model,e);assert.equal(p.mesh.pseudo_time,p.diagnostic.mesh_pseudo_time);}
    });
    assert.ok(events.slice(0,-1).every(e=>!api.position(model,e).complete));
  }
});
