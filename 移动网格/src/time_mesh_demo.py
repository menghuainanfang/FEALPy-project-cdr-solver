"""Continuous analytic-target mesh motion using the existing nested optimizer."""

import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.mesh import TriangleMesh
from matplotlib.tri import LinearTriInterpolator, Triangulation

from src.mesh_optimization import run as optimize
from src.variational_mesh import band_metric
from src.geometry import signed_triangle_areas


KEY_STEPS = {1,5,10,15,20,30,45,60,75,100}


def target_shift(time):
    return -.18*np.cos(2*np.pi*time/1.2)


def target_field(nodes, time):
    return np.tanh(25*(nodes[:,0]-.5-target_shift(time)-.15*np.sin(2*np.pi*nodes[:,1])))


def target_metric(centers, time):
    return band_metric(centers-np.array([target_shift(time),0.]))


def sample_error(nodes, cells, time):
    q=(np.arange(40)+.5)/40
    x,y=np.meshgrid(q,q);points=np.column_stack((x.ravel(),y.ravel()))
    interpolated=LinearTriInterpolator(Triangulation(nodes[:,0],nodes[:,1],cells),target_field(nodes,time))(points[:,0],points[:,1])
    return float(np.sqrt(np.mean((np.asarray(interpolated)-target_field(points,time))**2)))


def run(*, n=8, steps=60, dt=.02, inner_cap=40., max_outer=30,
        static=False, save_inner="key", on_step=None):
    bm.set_backend("numpy")
    mesh=TriangleMesh.from_box(nx=n,ny=n)
    reference,cells=np.array(mesh.node),np.array(mesh.cell)
    boundary=np.array(mesh.boundary_node_flag())
    nodes=reference.copy()
    nodes[~boundary]+=np.random.default_rng(42).uniform(-.2/n,.2/n,(sum(~boundary),2))
    def state(index,time,nodes,iterations,previous):
        return {"physical_step_index":index,"physical_time":time,"physical_dt":0. if static or index==0 else dt,
                "t_start":0. if static or index==0 else (index-1)*dt,"t_end":time,
                "nodes":nodes.tolist(),"values":target_field(nodes,time).tolist(),
                "metric":target_metric(nodes[cells].mean(axis=1),time).tolist(),
                "input_physical":previous.tolist(),"outer_iterations":iterations,
                "pde":{"status":"not_applicable","relative_nonlinear_residual":None,"l2_error":None},
                "target_origin":"prescribed time-dependent analytic target; no PDE solve",
                "target_position":float(.5+target_shift(time)),
                "target_interpolation_rmse":sample_error(nodes,cells,time),
                "fixed_grid_target_rmse":sample_error(reference,cells,time),
                "min_signed_area":float(signed_triangle_areas(nodes,cells).min()),
                "k_actual":len(iterations),"k_max":max_outer,"outer_tol":1.,
                "outer_status":"not_applicable" if not iterations else iterations[-1]["outer_status"],
                "outer_residual":None if not iterations else iterations[-1]["outer_residual"],
                "outer_converged":bool(iterations and iterations[-1]["outer_converged"]),
                "outer_stop_reason":"not_started" if not iterations else "converged" if iterations[-1]["outer_converged"] else "outer_iteration_cap"}
    states=[state(0,0.,nodes,[],nodes)]
    if on_step:on_step(states[0])
    metric_version=0
    for index in range(1,(1 if static else steps)+1):
        time=0. if static else index*dt
        previous=nodes.copy()
        result=optimize("band",n=n,initial_physical=nodes,perturbation=0.,
                        prescribed_metric=lambda centers:target_metric(centers,time),
                        epoch_prefix=f"{'static' if static else 'analytic'}_physical_step_{0 if static else index}",
                        inner_cap=inner_cap,max_outer=max_outer,
                        save_inner_meshes=save_inner=="all" or save_inner=="key" and index in KEY_STEPS,
                        verbose=False)
        iterations=[]
        for row in result["rounds"]:
            metric_version+=1
            iterations.append({"outer_iteration":row["outer_iteration"],"physical_time":time,
                               "metric_version":metric_version,"frozen_physical":row["input_physical"],
                               "output_physical":row["output_physical"],"frozen_metric":row["frozen_metric"],
                               "provisional_values":target_field(np.array(row["input_physical"]),time).tolist(),
                               "mesh_inner":row["inner"],"outer_converged":row["outer_converged"],
                               "outer_status":"converged" if row["outer_converged"] else "not_converged",
                               "outer_residual":max(row["mesh_change"]/1e-3,row["metric_change"]/1e-3),"outer_tol":1.,
                               "relative_mesh_change":row["mesh_change"],"relative_metric_change":row["metric_change"],
                               "outer_consecutive_success":row["consecutive_success"],"outer_consecutive_required":2,
                               "outer_unmet_reasons":row["unmet_reasons"]})
        nodes=np.array(result["final_physical"])
        states.append(state(0 if static else index,time,nodes,iterations,previous))
        states[-1]["overall_success"]=result["outer_converged"]
        if not result["outer_converged"]:states[-1]["outer_status"]="cap_reached"
        if on_step:on_step(states[-1])
        if not result["outer_converged"]:break
    return {"config":{"n":n,"steps":0 if static else steps,"physical_dt":0. if static else dt,
                      "inner_pseudo_time_cap":inner_cap,"max_outer":max_outer,"static":static,"seed":42,
                      "save_inner_traces":save_inner,"metric_relaxation":.5},
            "reference":reference.tolist(),"cells":cells.tolist(),"edges":np.array(mesh.entity("edge")).tolist(),
            "states":states,"method":"analytic moving band / nested FEALPy mesh optimization; no PDE solve"}
