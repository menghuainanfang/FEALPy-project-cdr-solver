"""Reproduce multidimensional CDR verification, figures and report."""
from pathlib import Path
import csv
import json
import platform
import importlib.metadata as metadata
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fealpy
from fealpy.functionspace import LagrangeFESpace
from pde import CDRData
from cdr_lfem_solver import CDRLFEMSolver, box_mesh

ROOT = Path(__file__).parent
OUT = ROOT/'results_extension'
OUT.mkdir(exist_ok=True)


def model(data, n, p):
    return CDRLFEMSolver(LagrangeFESpace(box_mesh(data, n), p=p), data)


def errors(s, uh, t=0):
    return (float(s.mesh.error(s.pde.at_time(s.pde.solution, t), uh, q=7)),
            float(s.mesh.error(s.pde.at_time(s.pde.gradient, t), uh.grad_value, q=7)))


def check_history(s):
    assert all(r['residual'] < 1e-10 and r['boundary'] < 1e-11 for r in s.history), s.history


def run():
    spatial = []
    for dim in (1, 2, 3):
        x = sp.symbols(f'x0:{dim}')
        exact = 1 + sum(x) + sp.prod(sp.sin(sp.pi*v) for v in x)
        data = CDRData(dim, a=1+sum(v*v for v in x), b=[1+v for v in x],
                       c=2+sum(x), exact=exact)
        for p in (1, 2):
            previous = None
            grids = (8, 16, 32, 64) if dim == 1 else ((4, 8, 16, 32) if dim == 2 else (2, 4, 8, 16))
            for n in grids:
                s = model(data, n, p)
                uh = s.solve()
                e = errors(s, uh)
                check_history(s)
                orders = np.log2(np.array(previous)/e) if previous is not None else [0, 0]
                row = dict(dim=dim, p=p, n=n, dofs=s.space.number_of_global_dofs(),
                           L2=e[0], H1=e[1], L2_order=float(orders[0]), H1_order=float(orders[1]), **s.history[-1])
                spatial.append(row)
                previous = e
                print('space', row, flush=True)
            assert orders[0] > p+.7 and orders[1] > p-.2, row
    temporal = []
    independent = []
    t = sp.Symbol('t')
    for dim in (1, 2, 3):
        x = sp.symbols(f'x0:{dim}')
        w = 1+sum(v*v for v in x)
        a = 1+sum(v*v for v in x)
        b = [1+v for v in x]
        c = 2+sum(x)
        d = 1+t+sum(x)
        # Hand-derived spatial operator on w; source bypasses automatic differentiation.
        Lw = -2*dim*a - 4*sum(v*v for v in x) + 2*sum((1+v)*v for v in x) + c*w
        linear = CDRData(dim, a=a, b=b, c=c, d=d, exact=(1+t)*w,
                         source=d*w+(1+t)*Lw)
        s = model(linear, 4, 2)
        uh = s.solve_time(.2, 4)
        e = errors(s, uh, .2)
        check_history(s)
        assert e[0] < 1e-10 and e[1] < 1e-9, e
        independent.append(dict(dim=dim, L2=e[0], H1=e[1], **s.history[-1]))
        data = CDRData(dim, a=a, b=b, c=c, d=d, exact=sp.exp(t)*w,
                       source=sp.exp(t)*(d*w+Lw))
        previous = None
        for steps in (5, 10, 20, 40):
            s = model(data, 4, 2)
            uh = s.solve_time(.5, steps)
            e = errors(s, uh, .5)
            check_history(s)
            order = float(np.log2(previous/e[0])) if previous is not None else 0.
            temporal.append(dict(dim=dim, steps=steps, dt=.5/steps, L2=e[0], H1=e[1], order=order,
                                 max_residual=max(r['residual'] for r in s.history),
                                 max_boundary=max(r['boundary'] for r in s.history)))
            previous = e[0]
            print('time', temporal[-1], flush=True)
        assert .85 < order < 1.15, temporal[-1]
    # A changed d with fixed f must change the transient solution.
    base = CDRData(1, d=1, exact='exp(t)*x0*(1-x0)')
    changed = CDRData(1, d=2, exact='exp(t)*x0*(1-x0)', source=base.expressions['source'])
    s1, s2 = model(base, 16, 2), model(changed, 16, 2)
    u1, u2 = s1.solve_time(.2, 20), s2.solve_time(.2, 20)
    capacity_effect = float(np.max(np.abs(u1[:]-u2[:])))
    assert capacity_effect > 1e-3
    # The d=1 steady limit reproduces the existing one-dimensional interface.
    from cdr_solver import solve_cdr
    old = solve_cdr(base.diffusion, base.convection, base.reaction,
                    base.source, base.dirichlet, n=16, p=2)[2]
    new = model(base, 16, 2).solve()
    regression = float(np.max(np.abs(old[:]-new[:])))
    assert regression < 1e-11
    # Export a changing 2D field and verify the actual VTK reader can read it.
    demo = CDRData(2, d='1+x0+x1', exact='exp(-t)*sin(pi*x0)*sin(pi*x1)')
    s = model(demo, 20, 1)
    uh = s.solve_time(1., 20, OUT/'vtu')
    check_history(s)
    import vtk
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(OUT/'vtu'/'solution_0020.vtu'))
    reader.Update()
    grid = reader.GetOutput()
    assert grid.GetNumberOfPoints() == s.mesh.number_of_nodes()
    assert grid.GetPointData().GetArray('u').GetNumberOfTuples() == s.mesh.number_of_nodes()
    for name, rows in [('spatial', spatial), ('temporal', temporal)]:
        with (OUT/f'{name}.csv').open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader(); writer.writerows(rows)
    results = dict(spatial=spatial, temporal=temporal, independent=independent,
                   capacity_effect=capacity_effect, regression=regression,
                   environment=dict(python=platform.python_version(), fealpy_path=fealpy.__file__,
                                    versions={k:metadata.version(k) for k in ('fealpy','numpy','scipy','sympy','vtk')}),
                   vtk_points=grid.GetNumberOfPoints(), vtk_cells=grid.GetNumberOfCells())
    (OUT/'verification.json').write_text(json.dumps(results, indent=2), encoding='utf-8')
    figures(spatial, temporal, s, uh)
    report(results)
    print('PASS: spatial/time convergence, independent polynomial, capacity effect, 1D regression, VTU readback', flush=True)


def figures(spatial, temporal, demo, uh):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for dim, ax in enumerate(axes, 1):
        for p in (1, 2):
            rows = [r for r in spatial if r['dim']==dim and r['p']==p]
            for key, marker in [('L2','o'), ('H1','s')]:
                ax.loglog([1/r['n'] for r in rows], [r[key] for r in rows], marker+'-', label=f'P{p} {key}')
        ax.set(title=f'{dim}D spatial convergence', xlabel='h', ylabel='Error')
        ax.grid(alpha=.3); ax.legend()
    fig.tight_layout(); fig.savefig(OUT/'spatial.png', dpi=170); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4))
    for dim in (1,2,3):
        rows = [r for r in temporal if r['dim']==dim]
        ax.loglog([r['dt'] for r in rows], [r['L2'] for r in rows], 'o-', label=f'{dim}D')
    ax.set(xlabel='Time step', ylabel='L2 error', title='Backward Euler: temporal convergence')
    ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(OUT/'temporal.png',dpi=170); plt.close(fig)
    import matplotlib.tri as mtri
    points = demo.mesh.node
    tri = mtri.Triangulation(points[:,0],points[:,1],demo.mesh.cell)
    values = np.empty(demo.mesh.number_of_nodes())
    values[np.asarray(demo.mesh.cell)] = np.asarray(uh(np.eye(3)))
    exact = demo.pde.solution(points,1.)
    fig, axes = plt.subplots(1,3,figsize=(13,4))
    for ax, val, title in zip(axes,[demo.pde.solution(points,0), values, values-exact],
                              ['Initial field','Numerical field, t=1','Vertex error, t=1']):
        m = ax.tripcolor(tri,val,shading='gouraud'); fig.colorbar(m,ax=ax)
        ax.set(title=title,xlabel='x',ylabel='y',aspect='equal')
    fig.tight_layout(); fig.savefig(OUT/'field.png',dpi=170); plt.close(fig)


def report(r):
    last = [v for v in r['spatial'] if v['n']==(64 if v['dim']==1 else 32 if v['dim']==2 else 16)]
    lines = ['# 多维与含时 CDR 求解器任务报告', '',
        '模型采用 $d(\\boldsymbol{x},t)u_t-\\nabla\\cdot(a\\nabla u)+\\boldsymbol b\\cdot\\nabla u+cu=f$。完成 1、2、3 维连续 P1/P2 空间离散、稳态求解、后向欧拉时间推进及 VTU/PVD 输出。数值结果见下表。', '',
        '## 模型与离散', '',
        '时间项中的 d 按独立容量系数解释，允许常数、空间或时间的符号表达式。这是依据沟通转述采用的模型约定；d=1 包含标准时间导数。采用 d u_t，不是 ∂t(d u)。扩散 a 为正标量，对流 b 为维数对应的向量，c 为标量。区域为区间、矩形或长方体，全 Dirichlet 边界。', '',
        '弱形式为 $\\int d u_t v+\\int a\\nabla u\\cdot\\nabla v+\\int(\\boldsymbol b\\cdot\\nabla u)v+\\int cuv=\\int fv$。稳态不组装时间质量矩阵。含时格式为 $(M_d^{n+1}+\\Delta t A^{n+1})U^{n+1}=M_d^{n+1}U^n+\\Delta t F^{n+1}$，边界值在新时刻施加。时间变化的 d 在等式两边均使用新时刻值。', '',
        '## 程序结构', '',
        '- `pde/cdr_data.py`：PDE 数据类，封装 a/b/c/d/f、边界、初值、精确解及梯度；SymPy 求导构造制造源项，NumPy 函数供积分器求值。',
        '- `cdr_lfem_solver.py`：网格工厂、积分器初始化、系数检查、空间组装、质量矩阵、边界处理、稳态与时间求解、结果导出。线性求解器在初始化时确定。',
        '- `verify_extension.py`：全部数值验证、CSV/JSON、图片和本报告生成。', '',
        '## 空间算例与验证', '',
        '取 $u=1+\\sum_i x_i+\\prod_i\\sin(\\pi x_i)$，$a=1+\\sum_i x_i^2$，$b_i=1+x_i$，$c=2+\\sum_i x_i$。由强形式生成源项，边界取精确解。1D 网格 n=8/16/32/64，2D 为 4/8/16/32，3D 为 2/4/8/16，每维均验证 P1/P2。h=1/n，误差阶为 log2(Eh/Eh/2)。H1 表示半范误差。', '',
        '|维数|阶次|n|自由度|L2误差|H1半范误差|L2阶|H1阶|', '|---|---|---|---|---|---|---|---|']
    lines += [f"|{v['dim']}|P{v['p']}|{v['n']}|{v['dofs']}|{v['L2']:.5e}|{v['H1']:.5e}|{v['L2_order']:.3f}|{v['H1_order']:.3f}|" for v in last]
    lines += ['', '![空间收敛](../results_extension/spatial.png)', '', '## 时间算例与独立验证', '',
        '令 $w=1+\\sum_i x_i^2$，保留上述 a/b/c，取 $d=1+t+\\sum_i x_i$。手工展开 $Lw=-2ma-4\\sum_i x_i^2+2\\sum_i(1+x_i)x_i+cw$，其中 m 为维数。', '',
        '独立多项式算例取 $u=(1+t)w$，显式输入 $f=dw+(1+t)Lw$。空间 P2 能精确表示 w，后向欧拉能精确表示线性时间变化；它同时检验变系数 d、非零时变边界及空间算子。', '',
        '|维数|L2误差|H1半范误差|', '|---|---|---|']
    lines += [f"|{v['dim']}|{v['L2']:.5e}|{v['H1']:.5e}|" for v in r['independent']]
    lines += ['', '时间收敛取 $u=e^t w$、显式源项 $f=e^t(dw+Lw)$，固定 n=4、P2，终止时间 0.5，步数 5/10/20/40。空间能表示该精确解，用于分离时间误差。', '', '|维数|步数|dt|L2误差|时间阶|', '|---|---|---|---|---|']
    lines += [f"|{v['dim']}|{v['steps']}|{v['dt']:.4f}|{v['L2']:.5e}|{v['order']:.3f}|" for v in r['temporal'] if v['steps']==40]
    lines += ['', '![时间收敛](../results_extension/temporal.png)', '',
        f"固定源项，仅将 d 从 1 改为 2，终态自由度最大差为 {r['capacity_effect']:.5e}，说明容量系数参与求解。一维稳态新旧接口最大差为 {r['regression']:.5e}。",
        '', '每次线性求解检查施加边界后的相对残差小于 1e-10，边界误差小于 1e-11。空间最终阶要求 L2>p+0.7、H1>p−0.2；时间阶要求 0.85 到 1.15。完整数据位于 CSV 和 verification.json。', '',
        '## 直观结果与 ParaView', '',
        '二维展示取 $u=e^{-t}\\sin(\\pi x)\\sin(\\pi y)$，d=1+x+y，a=1，b=0，c=0，n=20、P1、dt=0.05。左图是初值，中图是 t=1 的数值解，右图是顶点误差。', '',
        '![二维解与误差](../results_extension/field.png)', '',
        f"VTU 读取验证得到 {r['vtk_points']} 个顶点、{r['vtk_cells']} 个单元，并核对 u 数组长度。输出共 21 个时刻，PVD 记录实际物理时间。",
        '', '在 ParaView 中打开 `results_extension/vtu/solution.pvd`，点击 Apply，颜色字段选择 u，使用播放按钮查看变化。VTU 由 mesh.to_vtk 写入。导出的是网格顶点采样；P2 的单元内部高阶变化不包含在这种线性网格展示中。', '',
        '## 复现', '', '```powershell', 'python verify_extension.py', '```', '',
        f"Python {r['environment']['python']}；依赖：" + '，'.join(f'{k} {v}' for k,v in r['environment']['versions'].items()) + '。',
        f"FEALPy 加载位置：`{r['environment']['fealpy_path']}`。", '',
        '## 限制', '',
        '实现与验证范围为 NumPy/CPU、1/2/3 维盒形区域、连续 P1/P2、全 Dirichlet、正标量扩散和正容量。一般维数的数学形式不等于已支持高维网格。未实现张量扩散、混合边界、自适应时间步或对流稳定化；强对流工况仍可能振荡。积分点正性检查不是整个区域上的正性证明。后向欧拉为一阶时间格式。符号数据支持可由 NumPy 求值的 SymPy 表达式。', '']
    (ROOT/'docs'/'多维含时求解器_任务报告.md').write_text('\n'.join(lines), encoding='utf-8')


if __name__ == '__main__':
    run()
