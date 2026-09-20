# 从数组到有限元求解器：代码与 FEALPy 源码学习讲义

这份文档对应仓库中的 `cdr/` 包与 `example_vector.py`。读者可以没有 FEALPy 经验，但最好能读懂 Python 的函数、列表和类。数学式子的详细来源见[数学讲义](CDR数学推导与学习讲义.md)。本篇把“这一句在做什么”与“为什么要这样写”连在一起，文末附完整核心源码的带行号阅读本及逐段解释。

## 1. 先运行一个最小闭环

在仓库根目录打开终端，创建并激活虚拟环境，再安装依赖：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements_extension.txt
python example_vector.py --dim 2 --n 8 --p 2 --steps 20
```

程序计算两个耦合分量，到时间 0.5 停止，打印 $L^2$ 误差、梯度误差和最后一步代数残差，写入 `output/vector/solution.pvd`。其中 `n=8` 是每个坐标方向分 8 段，不是全网格只有 8 个单元。

模块关系如下：

```text
example_vector.py  配置问题并运行
    ├── cdr/cases.py  具体的数学数据
    │       └── cdr/pde.py  数据接口、时间绑定、可选符号制造解
    └── cdr/solver.py  网格、空间、装配、边界、时间推进、误差、VTU
            └── cdr/stabilization.py  SUPG 单元矩阵和单元载荷
```

`verify_vector.py` 跑数值验证并保存数据；`build_final_report.py` 用这些数据画图并写报告。画图脚本不重新求解方程。旧版根目录的 `cdr_lfem_solver.py` 等文件保留历史标量接口；新增功能的入口是 `from cdr import ...`。

## 2. 读代码前先认清这些 Python 写法

`class` 是把相关数据和函数放在一起的方式。`self.dim` 是当前这个问题对象的空间维数。两个对象可以使用不同维数和系数，彼此独立。`super().__init__` 调用父类初始化公共属性。

`@cartesian` 是装饰器。在这里它主要给函数标记“输入是笛卡尔坐标”，不是自动求导工具。`@staticmethod` 表示方法不使用实例的 `self`。

`p[...,0]` 中的 `...` 保留前面的所有轴，最后的 0 选择 x 坐标。例如 p 的形状为 `(8,6,2)`，表示 8 个单元、每单元 6 个点、每点 2 个坐标，则 x 的形状为 `(8,6)`。同一方法也可接收 `(100,2)` 的普通点集。

`a[...,None]` 在末尾增加长度为 1 的轴。若 a 是 `(NC,NQ)`，它变为 `(NC,NQ,1)`，就能与 `(NC,NQ,LDOF)` 的基函数数组逐点相乘。

`bm.stack([u1,u2],axis=-1)` 把两个标量场沿末轴合成向量场。`axis=-2` 用在梯度中，让结果排列为“点、分量、坐标方向”。`bm.broadcast_to` 把常数或共享数组扩展为需要的形状；它常返回视图，本项目只读取这些广播结果。

`@` 是矩阵乘法；`*` 是逐元素乘法。`uh[:] = values` 把系数填入 FEALPy 的有限元函数对象，保留该对象的函数空间信息。`uh(bcs)` 则在重心坐标处计算函数值，和取系数切片不是同一个操作。

## 3. 数据接口：把方程交给求解器

### 3.1 每个方法的约定

令 `S=p.shape[:-1]` 表示所有采样点的前置轴形状。

| 方法 | 数学意义 | 返回形状 |
| --- | --- | --- |
| `diffusion_coef(p,t)` | a | S |
| `diffusion_gradient(p,t)` | $\nabla a$ | S+(D,) |
| `convection_coef(p,t)` | b | S+(D,) |
| `reaction_coef(p,t)` | C | S+(m,m) |
| `capacity_coef(p,t)` | d | S |
| `source(p,t)` | f | S+(m,) |
| `dirichlet(p,t)` | 边界 g | S+(m,) |
| `initial(p)` | 初值 | S+(m,) |
| `solution(p,t)` | 精确解，仅误差验证需要 | S+(m,) |
| `gradient(p,t)` | 精确解梯度，仅误差验证需要 | S+(m,D) |

即使只有一个未知量，也保留末尾长度为 1 的分量轴。这样标量与向量使用同一条求解路径。直接编写数据类时，常系数也返回与采样点匹配的数组，例如 `bm.ones_like(p[...,0])`；不应把任意 Python 标量直接当成满足上述约定的场。

### 3.2 为什么普通类方法比动态改 self.source 更适合阅读

FEALPy 的 PDE 示例通常直接定义 `def source(self,p,...)`。编辑器能找到方法定义，子类可以自然重写，读者也能把源项公式与代码逐行对应。`CoupledData` 就采用这种方式。

`SymbolicCDR` 是可选的制造解工具：它在初始化时用 SymPy 求导和编译表达式，但 `source`、`solution`、`diffusion_coef` 等仍是显式方法。若题目已经给出了 f，不必知道精确解，也不必使用符号微分。

### 3.3 一个不用精确解的可运行例子

把下面内容保存为仓库根目录的 `my_case.py`：

```python
from fealpy.backend import backend_manager as bm
from cdr import CDRSolver, SymbolicCDR, box_mesh

bm.set_backend('numpy')
pde = SymbolicCDR(
    2, components=2,
    a='0.01*(1+t)', b=['1+t', '0.2'], d='1+t',
    reaction=[[2, -1], [-1, 2]],
    source=['exp(-20*((x-.3)**2+(y-.5)**2))', '0'],
    boundary=[0, 0], initial=[0, 0],
)
solver = CDRSolver(pde, box_mesh(pde, 12), p=1)
solver.solve_time(0.2, 20, output='output/my_case')
print(solver.history[-1])
```

这里第一分量有局部源，第二分量可以通过耦合受到影响。没有提供精确解，因此不要调用 `errors()`。网格和时间步仍须根据问题做收敛研究；这个例子展示接口，不替代该物理问题的精度验收。

### 3.4 一步步理解 CoupledData.source

代码先由 `_spatial` 得到 $w,\nabla w,\Delta w$，由 `_time` 得到 $s,s'$。精确解设为 $s(w,2w+1)^T$。

`spatial=-a*lap` 给出 $-a\Delta w$；再加 `(b-grad_a)*grad` 对末轴求和，得到 $(b-\nabla a)\cdot\nabla w$。合起来就是 $Lw=-\nabla\cdot(a\nabla w)+b\cdot\nabla w$。

`ut=ds*stack([w,2*w+1])` 是时间导数。第二分量的空间微分是第一分量的两倍，常数 1 的空间导数为零。最后加上 `einsum('...ij,...j->...i',C,u)`，就是在每个点做一次矩阵向量乘法。三个部分相加成为 f。

## 4. FEALPy 在底层做了什么

以下源码定位使用本机检出的 FEALPy 提交 `e033533bd8e24f796f9c5284685b20a98b4d11cb`。安装版 3.4.0 和该源码版都通过了本项目接口测试。链接固定提交，避免随着主分支更新而改变含义。只摘取理解本项目所需的路径。

### 4.1 PDE 示例与坐标标记

阅读 [parabolic2d.py](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/pde/parabolic2d.py) 的 `source(p,t)`，可看到先取 x、y 再计算公式的形式；[CDR model 的 Exp0001](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/model/diffusion_convection_reaction/exp0001.py) 提供了扩散、对流、反应和源项的组织参考。

坐标标记的实际作用在 [utils/utils.py 的 process_coef_func](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/utils/utils.py#L21)：如果是 Cartesian 函数，先 `mesh.bc_to_point(bcs)` 得到实际坐标，然后 `coef(ps)`；如果是重心坐标函数，则沿另一条调用路径。

标准积分器只需要一个关于 p 的函数，因此 `at_time` 先固定 t：

```python
@cartesian
def evaluate(p):
    value = function(p, t)
    return value if component is None else value[(...,) + component]
```

这里的闭包只是把已确定的时间传给数据对象。它不把 PDE 的 `self.source` 改写成匿名属性。`component=(i,j)` 表示选择 C 的一项；`component=(i,)` 选择 f 的一个分量。

### 4.2 从网格到基函数

[LagrangeFESpace](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/functionspace/lagrange_fe_space.py#L146) 的 `basis` 返回形函数值；`grad_basis` 转给网格的梯度形函数；`hess_basis` 转给网格的 Hessian 形函数。

重心坐标 bcs 与实际坐标 p 不同。三角形的重心坐标有三个数且和为 1；实际平面点只有 x、y 两个数。同一套 bcs 可供所有单元使用，`bc_to_point` 把它们映到各自单元中。

| 数组 | 本项目仿射单纯形中的形状 | 含义 |
| --- | --- | --- |
| bcs | (NQ,D+1) | 参考求积点 |
| points | (NC,NQ,D) | 实际空间点 |
| measure | (NC,) | 单元长度、面积或体积 |
| phi | (NC,NQ,LDOF) | 广播后的基函数值 |
| grad | (NC,NQ,LDOF,D) | 物理坐标梯度 |
| Hessian | (NC,NQ,LDOF,D,D) | 二阶导数矩阵 |
| cell_to_dof | (NC,LDOF) | 局部编号到全局编号 |
| 单元矩阵 | (NC,LDOF,LDOF) | 每单元一块小矩阵 |

原始标量 `basis` 常为 `(1,NQ,LDOF)`，因为所有仿射单元共享参考基函数值；代码把它广播到每个单元。梯度与单元几何有关，不能同样忽略单元轴。

### 4.3 向量空间的排列约定

[TensorFunctionSpace 的初始化](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/functionspace/tensor_space.py#L13) 使用 `-1` 标记自由度轴位置。`shape=(m,-1)` 对应本项目的分量块排列。

若节点函数值矩阵为

```text
[[10, 100],
 [20, 200],
 [30, 300]]
```

则 `.T.reshape(-1)` 得到 `[10,20,30,100,200,300]`，与 `BlockForm` 的块顺序一致。如果直接 `.reshape(-1)` 会得到另一种排列，使耦合项和边界位置错位。

`boundary_interpolate`、`value`、`grad_value` 都依赖这个约定。特别地，向量函数的梯度返回 `(NC,NQ,m,D)`，所以误差积分能同时对分量和方向求和。

### 4.4 标准积分器怎样表示一项积分

[ScalarConvectionIntegrator](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/fem/scalar_convection_integrator.py#L33) 的 `fetch` 取得积分点、权重、基函数、梯度和单元度量；`assembly` 求出系数，把 b 与试探函数梯度内积，然后用 `bilinear_integral` 与测试函数相乘积分。

本项目不重写这些已有项。对应关系为：

| 弱形式 | FEALPy 类 | 构造关键参数 |
| --- | --- | --- |
| $\int a\nabla u\cdot\nabla v$ | ScalarDiffusionIntegrator | `coef=...` |
| $\int(b\cdot\nabla u)v$ | ScalarConvectionIntegrator | `coef=...` |
| $\int cuv$、$\int duv$ | ScalarMassIntegrator | `coef=...` |
| $\int fv$ | ScalarSourceIntegrator | `source=...` |

质量积分器不只用于时间项，因为反应项与质量项的积分结构相同。

### 4.5 自定义 SUPG 为什么能进入同一个 Form

[ConstIntegrator](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/fem/integrator.py#L264) 接收已经算好的局部数组 `value` 和局部到全局映射 `to_gdof`。它的 `assembly` 返回该数组，`to_global_dof` 返回映射。

这里的“Const”表示交给它的数组在本次装配中已确定，不表示整个 PDE 的系数恒定。每个新时刻重新计算 SUPG 数组，再新建对应积分器。

`BilinearForm.add_integrator(...)` 收集积分项，`assembly()` 求和并散装成稀疏矩阵。`LinearForm` 对应载荷向量。库内 COO 表示非零元素的行、列与值，`coalesce()` 合并共享自由度贡献，CSR 更适于后续矩阵运算。我们没有把稀疏矩阵转成稠密矩阵。

### 4.6 块矩阵如何拼接

[BlockForm.assembly](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/fem/block_form.py#L67) 先计算每个块的行列偏移，然后逐块装配。例如第二分量的所有行号加上标量自由度数 $N_s$，第二分量的列号也加 $N_s$。最终合并 COO 项并转为 CSR。

因此 `forms[i][j]` 必须表示第 i 个方程对第 j 个未知分量的作用。非对角块只有反应耦合及其 SUPG 修正；本模型没有交叉扩散。

### 4.7 边界消元与稀疏求解

[DirichletBC.apply](https://github.com/weihuayi/fealpy/blob/e033533bd8e24f796f9c5284685b20a98b4d11cb/fealpy/fem/dirichlet_bc.py#L101) 先处理向量，再处理矩阵。`apply_vector` 插入边界数据并扣除其对内部方程的影响；`apply_matrix` 清理边界行列，设置单位对角项。这个操作不消除对流带来的非对称性。

`spsolve(A,F,solver='scipy')` 使用库提供的稀疏求解入口。`linear_solver` 回调可以替换它；输入是 FEALPy 稀疏张量和后端数组，返回必须是同形状有限向量。向回调传副本，使残差检查仍使用未被回调修改的 A 和 F。

## 5. SUPG 的每一类数组怎么计算

打开 `cdr/stabilization.py` 对照以下步骤。

1. 取得 bcs、权重、物理点和度量。这四项定义了“在哪些点、用什么权重积分”。
2. 取得 phi 和 grad，并调用 PDE 类取得新时刻 a、b、C、d。所有系数具有显式点轴。
3. `speed=norm(b,axis=-1)` 是每个求积点的速度模长。`mesh.grad_lambda()` 是 P1 重心坐标梯度，用于计算沿流向的单元尺寸，即使当前空间是 P2 也仍用它。
4. `directional=einsum('cqd,cid->cqi',b,linear_grad)` 对空间方向 d 求和，保留单元 c、点 q、重心基函数 i。绝对值求和得到尺度分母。
5. 计算 h、Pe 和 tau。小 Pe 的级数避免两个大数相减；无流动时 tau 为零。分母占位是为避免 `where` 两分支求值引入零除。
6. `streamline=einsum('cqd,cqid->cqi',b,grad)` 得到每个基函数的方向导数。`test=tau[...,None]*streamline` 是新增的测试权。
7. P2 的 `einsum('cqidd->cqi',hessian)` 取 Hessian 对角线并求和，得到 Laplacian。P1 的单元内 Laplacian 为零。再加入 $-a\Delta\phi-\nabla a\cdot\nabla\phi$。
8. 对每一对分量组装反应耦合；只有同一分量才加扩散和对流残差。
9. 容量和源项分别产生时间质量修正和载荷修正。返回的仍是单元数组，之后交给 FEALPy 全局装配。

最重要的一句是

```python
bm.einsum('q,c,cqi,cqj->cij', w, measure, test, trial)
```

它等价于

$$S[c,i,j]=|K_c|\sum_qw[q]T[c,q,i]R[c,q,j].$$

箭头右侧没有 q，所以对 q 求和；c,i,j 被保留。这个公式正好对应数学讲义中的单元积分，不是另一种数值方法。

## 6. 时间循环怎样避免混淆历史解

进入第 k 步前，`old` 表示 $U^{k-1}$，`previous` 在第二步开始表示 $U^{k-2}$。`assemble(t)` 在新时刻计算 A、M、F。第一步用 BE，后面使用

```python
left = 1.5*M + dt*A
right = M @ (2*old - 0.5*previous) + dt*F
```

解完后执行 `previous,old=old,uh[:].copy()`。Python 先计算右侧再赋值，所以旧的 old 能正确保存为 previous；copy 保证下一次写入有限元函数时不会把历史数组一起改变。

每次调用 `solve_time` 都从 `initial` 重新开始并重置 history，不是从上次终态续算。第一步会用 t=0 的边界覆盖初值边界位置；实际建模应尽量提供相容的初边值。

这里没有把 A、M 缓存为时间常数。代价是每步重新组装；收益是随时间变化的系数与边界不会被错误复用。未来只有在明确判断哪些项不变后，才适合增加装配或分解复用。

## 7. 结果和文件如何产生

`errors` 调用网格积分器，在求积点比较精确函数与 `uh`，再比较精确梯度与 `uh.grad_value`。它不是只比较节点，所以能发现“节点吻合但单元内部有误差”的情况。

`export_vtu` 传入 `eye(D+1)` 作为重心坐标，这些点恰好是单元顶点。`uh(...)` 得到每个单元顶点处的分量值，再按网格 cell 的节点编号放回全局节点数组。连续有限元在共享顶点上的值一致。

输出先浅复制网格，再单独复制 `nodedata` 字典；新写入 `u_0`、`u_1` 等数组不会替换调用者的原字典条目。几何和拓扑共享但不修改。`u_magnitude` 是欧氏模长；最多三分量时额外写三分量 `u`，缺失位置补零。

每步 VTU 文件名含补零序号，例如 `solution_0007.vtu`；PVD 用真实时间记录文件。打开 PVD 可以直接播放。详见[ParaView 指南](向量含时结果_ParaView指南.md)。

## 8. 验证脚本不只是“运行没报错”

| 实验 | 对应文件位置 | 能发现的错误 |
| --- | --- | --- |
| 多维 P1/P2 空间阶 | verify_vector.py 的 spatial 循环 | 形函数、维数、符号、耦合装配错误 |
| 空间二次、时间一次再现 | polynomial 循环 | 时间质量、边界及 SUPG 缺项 |
| BE/BDF2 时间阶 | temporal 循环 | 历史系数或新旧时刻混用 |
| 手工 f 与独立符号微分 | source crosscheck | 手工求导或坐标方向错误 |
| 固定 f 去掉耦合 | coupling_effect | 非对角块被遗漏 |
| Galerkin/SUPG 边界层 | strong 循环 | 稳定项失效或符号错误 |
| 层的网格加密 | resolved_layer | 把节点拟合误当成场精确 |
| VTK 读回 | vtu_vector_values | 导出分量排列或数据错误 |
| 共享网格、重复运行 | verify_vector_interfaces.py | 输入数据污染或历史状态串用 |
| 四分量无解析解 | 同上 | 把空间维数与分量数混淆、误依赖精确解 |
| 导入时保持其他后端 | 同上 | 全局后端或模块命名空间污染 |

完整运行 `python verify_vector.py` 后会生成 JSON 和 CSV；断言失败会终止，不写出完整通过记录。接口脚本的额外 PyTorch 依赖只用于模拟调用者先选择了另一个后端，求解本身仍限定 NumPy/CPU。

## 9. 使用限制与排查顺序

本实现支持 1—3 维全维仿射单纯形、P1/P2、盒形全 Dirichlet 边界、共享正标量 a/d、共享 b 和反应矩阵 C。它不支持曲面、非线性、交叉/张量扩散或任意边界类型。导出是顶点采样；SUPG 不保证任意网格的正性。详细数值限制见任务报告。

如果出现数组维数错误，先检查 PDE 输出形状与上面的表，而不是修改 FEALPy 源码。若边界错误，先检查 `shape=(m,-1)` 与 `.T.reshape(-1)` 是否一致，再看边界掩码。残差小但误差大时，检查 f 的符号、对流是否写成了散度形式，以及误差比较的时间是否相同。非有限数先看系数和源项的表达式是否溢出。

## 10. 自己修改一次模型

可以从 `CoupledData` 复制一个数据类，保留形状约定，修改 a、b、C 和精确解，手工推导 source。先用能被 P2 表示的二次函数验证，再用三角函数做网格加密。若只改 C 而仍想保留同一个制造解，也必须随之修改 f。

学习顺序建议为：先跑示例，再读 PDE 数据类，手算一维单元矩阵，阅读 assemble，最后阅读时间循环和 SUPG。下面的源码阅读本提供与文件一致的行号，方便逐段定位。


## 11. 完整核心源码阅读本

下面覆盖核心包和示例的全部源码。每段标注原文件行号，代码块保持原文，段后解释该组语句的作用。连续的参数列表、括号和空行随所属语句一起解释。

### cdr/__init__.py

#### 包入口（第 1—6 行）

```python
"""Vector convection-diffusion-reaction finite elements."""

from .solver import CDRSolver, box_mesh
from .pde import CDRModel, SymbolicCDR

__all__ = ["CDRSolver", "box_mesh", "CDRModel", "SymbolicCDR"]
```

模块说明用于说明职责；两个相对导入把求解器与数据类放到包的公共入口。__all__ 列出公开名称，因此使用者可以写 from cdr import CDRSolver。这里没有设置后端，也没有启动计算。

### example_vector.py

#### 导入（第 1—7 行）

```python
"""Run a coupled time-dependent CDR example and export a ParaView series."""

import argparse
from pathlib import Path
from fealpy.backend import backend_manager as bm
from cdr import CDRSolver, box_mesh
from cdr.cases import CoupledData
```

argparse 读取命令行参数；Path 表示文件路径；bm 是 FEALPy 的后端管理器；最后两行导入当前求解器和具体算例。

#### 定义命令行接口（第 8—18 行）

```python


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dim", type=int, choices=(1, 2, 3), default=2)
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--p", type=int, choices=(1, 2), default=2)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--scheme", choices=("bdf2", "backward_euler"), default="bdf2")
    parser.add_argument("--output", type=Path, default=Path("output/vector"))
    args = parser.parse_args()
```

main 把运行行为收在函数内。dim 和 p 用 choices 限定已验证范围；n 控制每方向网格数；steps 控制时间步数；scheme 选择差分格式；output 转为 Path。parse_args 将命令行文字转为有类型的参数对象。

#### 运行一遍求解流程（第 19—25 行）

```python
    bm.set_backend("numpy")
    pde = CoupledData(args.dim)
    solver = CDRSolver(pde, box_mesh(pde, args.n), p=args.p)
    solver.solve_time(0.5, args.steps, scheme=args.scheme, output=args.output)
    print("L2 / H1 seminorm:", solver.errors(0.5))
    print("Last step:", solver.history[-1])
    print("ParaView:", args.output / "solution.pvd")
```

选择 NumPy 是应用程序入口的决定，不在库导入时进行。CoupledData 创建方程数据；box_mesh 创建匹配维数的网格；CDRSolver 创建离散空间。solve_time 推进到 0.5，errors 用同一个终止时间比较。history[-1] 取最后一步记录，最后打印 PVD 路径。

#### 脚本入口（第 26—29 行）

```python


if __name__ == "__main__":
    main()
```

只有直接执行此文件时才调用 main。其他 Python 程序 import 它，不会自动计算或写文件。

### cdr/pde.py

#### 公共数据契约（第 1—14 行）

```python
"""PDE interfaces: d*u_t - div(a*grad(u)) + b.grad(u) + C*u = f."""

import sympy as sp
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian


class CDRModel:
    """Shared positive scalar a,d; spatial velocity b; component reaction C.

    solution/source/dirichlet return (..., components), even for one component.
    gradient returns (..., components, dim). Coefficients depend on (p,t).
    """

```

SymPy 用于符号表达式；bm 用于数组；cartesian 标记坐标类型。CDRModel 的说明约定了共享系数、分量轴和梯度轴，这是后续所有数组操作的前提。

#### 维数、分量数和区域（第 15—23 行）

```python
    def __init__(self, dim, components=1, box=None):
        if dim not in (1, 2, 3) or components < 1:
            raise ValueError("Use dimensions 1/2/3 and a positive component count")
        self.dim, self.components = dim, components
        self.box = list(box) if box is not None else [0.0, 1.0] * dim

    def domain(self):
        return list(self.box)

```

初始化保存 dim、components；默认 box 是 [0,1] 重复 dim 次。list(box) 保存独立列表；domain 再返回列表副本，使使用者修改返回值时不会改变数据对象中的区域。

#### 边界、初值和默认精确解（第 24—39 行）

```python
    @cartesian
    def is_dirichlet_boundary(self, p):
        mask = bm.zeros(p.shape[:-1], dtype=bm.bool)
        for axis in range(self.dim):
            mask |= bm.abs(p[..., axis] - self.box[2 * axis]) < 1e-12
            mask |= bm.abs(p[..., axis] - self.box[2 * axis + 1]) < 1e-12
        return mask

    @cartesian
    def dirichlet(self, p, t=0.0):
        return self.solution(p, t)

    @cartesian
    def initial(self, p):
        return self.solution(p, 0.0)

```

mask 的每个位置对应一个采样点；逐坐标检测盒子两侧并用逻辑或累积。容差用于浮点坐标判断。默认 dirichlet 使用精确解，initial 使用 t=0 的精确解；真实无解析解问题应重写这两个方法。

#### 固定时间的系数函数（第 40—50 行）

```python

def at_time(function, t, component=None):
    """Bind time while preserving the Cartesian coefficient convention."""

    @cartesian
    def evaluate(p):
        value = function(p, t)
        return value if component is None else value[(...,) + component]

    return evaluate

```

外层 at_time 接收数据方法及时间，内层 evaluate 只接收 p。component 为 None 就保留完整结果，否则用 (...,)+component 选择最后一个或两个轴。返回函数而非函数值，积分器之后在自己的求积点调用它。

#### 符号问题的入口与符号统一（第 51—87 行）

```python

class SymbolicCDR(CDRModel):
    """Optional NumPy symbolic data; public PDE functions remain class methods.

    Use x,y,z,t symbols. reaction is an m-by-m matrix; scalar reaction means c*I.
    Supply source/boundary/initial for physical data without an exact solution.
    """

    def __init__(
        self,
        dim,
        *,
        exact=None,
        source=None,
        boundary=None,
        initial=None,
        a=1,
        b=None,
        reaction=0,
        d=1,
        components=None,
        box=None,
    ):
        count = len(exact) if exact is not None else components
        if count is None:
            raise ValueError("Specify components when no exact solution is supplied")
        super().__init__(dim, count, box)
        x, y, z, t = sp.symbols("x y z t")
        coordinates = (x, y, z)[:dim]
        symbols = {str(v): v for v in (x, y, z, t)}

        def expr(value):
            value = sp.sympify(value)
            return value.xreplace(
                {v: symbols[str(v)] for v in value.free_symbols if str(v) in symbols}
            )

```

关键字参数减少多个系数位置写错的机会。若有 exact，就由列表长度确定分量数，否则必须给 components。symbols 建立 x,y,z,t 标准符号；expr 把字符串或表达式转为 SymPy 对象，并把同名但带不同假设的外来符号替换成标准符号，保证求导针对同一个变量。

#### 标量和矩阵系数（第 88—98 行）

```python
        a, d = expr(a), expr(d)
        b = [expr(v) for v in (b if b is not None else [0] * dim)]
        C = (
            sp.Matrix(reaction)
            if isinstance(reaction, (list, tuple, sp.MatrixBase))
            else sp.eye(count) * expr(reaction)
        )
        C = C.applyfunc(expr)
        if C.shape != (count, count) or len(b) != dim:
            raise ValueError("Coefficient dimensions do not match the PDE")
        u = [expr(v) for v in exact] if exact is not None else None
```

a、d 逐个规范化，b 默认为零向量。reaction 若是矩阵数据就保留其结构，否则按 cI 解释标量反应。applyfunc 规范化每个矩阵元素。维数检查防止 b 的长度或 C 的行列与问题不符。u 是解析分量列表，也允许不存在。

#### 按强形式制造源项（第 99—110 行）

```python
        if source is None:
            if u is None:
                raise ValueError("Provide source or exact solution")
            source = [
                d * sp.diff(u[i], t)
                + sum(
                    -sp.diff(a * sp.diff(u[i], v), v) + b[j] * sp.diff(u[i], v)
                    for j, v in enumerate(coordinates)
                )
                + sum(C[i, j] * u[j] for j in range(count))
                for i in range(count)
            ]
```

没有 source 时必须有 u。对每个分量先做 d*u_t；每个空间方向上做 -∂(a*∂u)/∂x 加 b*∂u；再把 C 的一行与所有 u 分量相乘求和。这里直接对乘积求导，因此包含变扩散系数的导数。

#### 边界、初值和编译器（第 111—123 行）

```python
        boundary = u if boundary is None else boundary
        initial = [v.subs(t, 0) for v in u] if initial is None and u is not None else initial
        if boundary is None:
            raise ValueError("Dirichlet data are required")
        self._functions = {}

        def compile_field(name, values, shape):
            expressions = [expr(v) for v in values]
            self._functions[name] = (
                [sp.lambdify((x, y, z, t), v, "numpy") for v in expressions],
                shape,
            )

```

默认边界取解析解；默认初值把 t 替换为零。没有解析解时必须明确提供边界。_functions 是私有字典。compile_field 把每个标量表达式 lambdify 成 NumPy 函数，同时记录该场末尾应有的形状。

#### 注册全部场（第 124—135 行）

```python
        compile_field("a", [a], ())
        compile_field("d", [d], ())
        compile_field("b", b, (dim,))
        compile_field("C", list(C), (count, count))
        compile_field("grad_a", [sp.diff(a, v) for v in coordinates], (dim,))
        for name, values in [("f", source), ("g", boundary), ("initial", initial), ("u", u)]:
            if values is not None:
                if len(values) != count:
                    raise ValueError(f"{name} needs {count} components")
                compile_field(name, values, (count,))
        if u is not None:
            compile_field("gradient", [sp.diff(v, c) for v in u for c in coordinates], (count, dim))
```

a、d 的末尾形状为空，b 为 dim，C 为 count×count，grad_a 为 dim。f、g、initial、u 都要求 count 个分量。解析梯度按“先分量、后坐标”展开，后续再 reshape 成 (count,dim)。

#### 将符号函数变为数组值（第 136—145 行）

```python

    def _evaluate(self, name, p, t):
        x = p[..., 0]
        y = p[..., 1] if self.dim > 1 else bm.zeros_like(x)
        z = p[..., 2] if self.dim > 2 else bm.zeros_like(x)
        functions, shape = self._functions[name]
        values = [
            bm.broadcast_to(bm.asarray(f(x, y, z, t), dtype=p.dtype), x.shape) for f in functions
        ]
        return bm.stack(values, axis=-1).reshape(x.shape + shape) if shape else values[0]
```

取 x，必要时取 y、z；缺少的坐标传零，因为 lambdify 统一接收四个参数。常数函数可能返回标量，broadcast_to 把它扩展到所有点。stack 将标量表达式组装，再 reshape 成约定的末尾轴。p 应是浮点坐标数组，保留其 dtype。

#### 显式公开方法（第 146—187 行）

```python

    @cartesian
    def diffusion_coef(self, p, t=0.0):
        return self._evaluate("a", p, t)

    @cartesian
    def diffusion_gradient(self, p, t=0.0):
        return self._evaluate("grad_a", p, t)

    @cartesian
    def convection_coef(self, p, t=0.0):
        return self._evaluate("b", p, t)

    @cartesian
    def reaction_coef(self, p, t=0.0):
        return self._evaluate("C", p, t)

    @cartesian
    def capacity_coef(self, p, t=0.0):
        return self._evaluate("d", p, t)

    @cartesian
    def source(self, p, t=0.0):
        return self._evaluate("f", p, t)

    @cartesian
    def solution(self, p, t=0.0):
        return self._evaluate("u", p, t)

    @cartesian
    def gradient(self, p, t=0.0):
        return self._evaluate("gradient", p, t)

    @cartesian
    def dirichlet(self, p, t=0.0):
        return self._evaluate("g", p, t)

    @cartesian
    def initial(self, p):
        if "initial" not in self._functions:
            raise ValueError("Transient problems require initial data")
        return self._evaluate("initial", p, 0.0)
```

每个 @cartesian 方法有固定数学职责，分别从私有字典中读取 a、grad_a、b、C、d、f、u、gradient、g。它们的公共形式和手工数据类相同。initial 单独检查是否提供初值；没有解析解的模型也能求解，但不能调用依赖 u 或 gradient 的误差评估。

### cdr/cases.py

#### 边界层数据类（第 1—14 行）

```python
"""Explicit coordinate-based FEALPy data for a coupled manufactured system."""

from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from .pde import CDRModel


class LayerData(CDRModel):
    """Stable evaluation of -eps*u_xx+u_x=0 on the unit box."""

    def __init__(self, dim=1, eps=0.001):
        super().__init__(dim, 1)
        self.eps = eps

```

继承 CDRModel 并设 components=1，保存正扩散参数 eps。虽然测试允许 1—3 维，精确解只沿 x 方向变化。

#### 稳定的指数值与梯度（第 15—27 行）

```python
    @cartesian
    def solution(self, p, t=0.0):
        x = p[..., 0]
        tail = bm.exp(bm.asarray(-1 / self.eps))
        return ((bm.exp((x - 1) / self.eps) - tail) / (1 - tail))[..., None]

    @cartesian
    def gradient(self, p, t=0.0):
        value = bm.zeros(p.shape[:-1] + (1, self.dim), dtype=p.dtype)
        value[..., 0, 0] = bm.exp((p[..., 0] - 1) / self.eps) / (
            self.eps * (-bm.expm1(bm.asarray(-1 / self.eps)))
        )
        return value
```

solution 的指数 (x-1)/eps 在单位盒内不为正，避免计算巨大指数；末尾 None 保留单分量轴。gradient 先创建 (点轴,1,dim) 的零数组，只给分量0、方向0赋值。expm1(v) 精确计算 exp(v)-1，分母负号恢复 1-exp(-1/eps)。

#### 边界层各项系数（第 28—54 行）

```python

    @cartesian
    def diffusion_coef(self, p, t=0.0):
        return bm.full_like(p[..., 0], self.eps)

    @cartesian
    def diffusion_gradient(self, p, t=0.0):
        return bm.zeros_like(p)

    @cartesian
    def convection_coef(self, p, t=0.0):
        value = bm.zeros_like(p)
        value[..., 0] = 1.0
        return value

    @cartesian
    def reaction_coef(self, p, t=0.0):
        return bm.zeros(p.shape[:-1] + (1, 1), dtype=p.dtype)

    @cartesian
    def capacity_coef(self, p, t=0.0):
        return bm.ones_like(p[..., 0])

    @cartesian
    def source(self, p, t=0.0):
        return bm.zeros(p.shape[:-1] + (1,), dtype=p.dtype)

```

full_like 生成常扩散；其空间梯度为零。b 只有 x 分量等于1。C 的末尾形状即使为单分量也保留 (1,1)。容量为1，源项为0，数组形状都与点集一致。

#### 耦合算例及显式坐标（第 55—71 行）

```python

class CoupledData(CDRModel):
    """Two components with time-dependent a,b,C,d in dimensions 1/2/3."""

    def __init__(self, dim=2, profile="sine", temporal="exp"):
        super().__init__(dim, 2)
        self.profile, self.temporal = profile, temporal

    def _spatial(self, p):
        x = p[..., 0]
        coordinates = [x]
        if self.dim > 1:
            y = p[..., 1]
            coordinates.append(y)
        if self.dim > 2:
            z = p[..., 2]
            coordinates.append(z)
```

CoupledData 固定两个分量，profile 选择空间形状，temporal 选择时间因子。_spatial 先拿 x，再根据 dim 增加 y、z，所以一维不会访问不存在的坐标轴。

#### 两种空间形状（第 72—80 行）

```python
        if self.profile == "quadratic":
            w = 1 + sum(v * v for v in coordinates)
            grad = bm.stack([2 * v for v in coordinates], axis=-1)
            lap = bm.full_like(x, 2 * self.dim)
        else:
            product = bm.ones_like(x)
            for v in coordinates:
                product *= bm.sin(bm.pi * v)
            w = 1 + sum(coordinates) + product
```

quadratic 分支设置 w=1+Σx²、grad=2x、lap=2dim；这能被 P2 表示。sine 分支逐方向累乘 sin(pi*x)，再加入常数及坐标和，为空间收敛提供非多项式光滑函数。

#### 三角函数的梯度与 Laplacian（第 81—90 行）

```python
            derivatives = []
            for i, v in enumerate(coordinates):
                term = bm.pi * bm.cos(bm.pi * v)
                for j, other in enumerate(coordinates):
                    if j != i:
                        term *= bm.sin(bm.pi * other)
                derivatives.append(1 + term)
            grad = bm.stack(derivatives, axis=-1)
            lap = -self.dim * bm.pi**2 * product
        return w, grad, lap
```

对每个方向 i，把该方向的 sin 求导为 pi*cos，其他方向仍乘 sin；线性部分的导数加1。每次二阶方向导数都给出 -pi²*product，dim 个方向相加得到 lap。返回三个数组供精确解和源项共同使用。

#### 时间因子和导数（第 91—100 行）

```python

    def _time(self, t):
        if self.temporal == "steady":
            return 1.0, 0.0
        return (
            (1 + t, 1.0)
            if self.temporal == "linear"
            else (bm.exp(bm.asarray(t)), bm.exp(bm.asarray(t)))
        )

```

steady 返回 (1,0)，linear 返回 (1+t,1)，exp 返回 (exp(t),exp(t))。第二项始终是第一项对 t 的导数，因此 source 不需要猜测时间变化。

#### 随时间变化的所有系数（第 101—120 行）

```python
    @cartesian
    def diffusion_coef(self, p, t=0.0):
        return (1 + t) * (1 + bm.sum(p * p, axis=-1))

    @cartesian
    def diffusion_gradient(self, p, t=0.0):
        return 2 * (1 + t) * p

    @cartesian
    def convection_coef(self, p, t=0.0):
        return (1 + t) * (1 + p)

    @cartesian
    def reaction_coef(self, p, t=0.0):
        C = bm.asarray([[2.0, -0.5], [-0.25, 3.0]])
        return bm.broadcast_to((1 + t) * C, p.shape[:-1] + (2, 2))

    @cartesian
    def capacity_coef(self, p, t=0.0):
        return 1 + t + bm.sum(p, axis=-1)
```

扩散同时随空间和时间变化，其梯度是2(1+t)p。对流对每个坐标取 (1+t)(1+x_k)。反应矩阵包含两个非零非对角项，广播到所有点。容量 1+t+Σx 在本实验时空区域内为正。

#### 解析解和解析梯度（第 121—133 行）

```python

    @cartesian
    def solution(self, p, t=0.0):
        w, _, _ = self._spatial(p)
        s, _ = self._time(t)
        return s * bm.stack([w, 2 * w + 1], axis=-1)

    @cartesian
    def gradient(self, p, t=0.0):
        _, grad, _ = self._spatial(p)
        s, _ = self._time(t)
        return s * bm.stack([grad, 2 * grad], axis=-2)

```

solution 把 w 与2w+1堆叠到末轴后乘s。gradient 把 grad 与2grad堆叠到倒数第二轴，末轴仍留给空间方向。下划线表示这一次不需要的返回值。

#### 手工源项的完整表达式（第 134—147 行）

```python
    @cartesian
    def source(self, p, t=0.0):
        w, grad, lap = self._spatial(p)
        s, ds = self._time(t)
        spatial = -self.diffusion_coef(p, t) * lap
        spatial += bm.sum(
            (self.convection_coef(p, t) - self.diffusion_gradient(p, t)) * grad, axis=-1
        )
        ut = ds * bm.stack([w, 2 * w + 1], axis=-1)
        return (
            self.capacity_coef(p, t)[..., None] * ut
            + s * bm.stack([spatial, 2 * spatial], axis=-1)
            + bm.einsum("...ij,...j->...i", self.reaction_coef(p, t), self.solution(p, t))
        )
```

spatial 依次组装 -aΔw 和 (b-∇a)·∇w；ut 使用时间因子导数。返回值三部分分别是容量乘时间导数、两个分量的空间算子、反应矩阵乘向量解。einsum 的 ... 保留任意前置点轴。

### cdr/stabilization.py

#### 几何、求积与基函数（第 1—13 行）

```python
"""Consistent cell residuals for streamline-upwind Petrov-Galerkin FEM."""

from fealpy.backend import backend_manager as bm


def supg_terms(space, pde, t, q):
    """Return local (spatial blocks, time mass, loads) using tau*b.grad(v)."""
    mesh = space.mesh
    bcs, w = mesh.quadrature_formula(q, "cell").get_quadrature_points_and_weights()
    points = mesh.bc_to_point(bcs)
    measure = mesh.entity_measure("cell")
    grad = space.grad_basis(bcs)
    phi = bm.broadcast_to(space.basis(bcs), grad.shape[:-1])
```

space.mesh 取得网格，quadrature_formula 提供参考点和权重；bc_to_point 转换为实际坐标。measure 分别是长度、面积或体积。grad 有单元轴，phi 广播到相同的点与局部自由度轴。

#### 系数与流线尺寸（第 14—24 行）

```python
    a, b, C, d = (
        pde.diffusion_coef(points, t),
        pde.convection_coef(points, t),
        pde.reaction_coef(points, t),
        pde.capacity_coef(points, t),
    )
    speed = bm.linalg.norm(b, axis=-1)
    # Streamline length uses linear barycentric gradients, independent of p.
    linear_grad = mesh.grad_lambda()
    directional = bm.einsum("cqd,cid->cqi", b, linear_grad)
    denominator = bm.sum(bm.abs(directional), axis=-1)
```

同一组 points,t 求出全部系数。speed 对 b 的方向轴求范数。grad_lambda 是线性重心坐标梯度；einsum 对空间方向求和；sum(abs(...)) 形成沿流向尺度的分母。

#### tau 的数值计算（第 25—31 行）

```python
    active = speed > 1e-14
    h = 2 * speed / bm.where(denominator > 0, denominator, 1.0) / space.p
    pe = speed * h / (2 * a)
    small = pe < 1e-3
    safe_pe = bm.where(small, 1.0, pe)
    langevin = bm.where(small, pe / 3 - pe**3 / 45, 1 / bm.tanh(safe_pe) - 1 / safe_pe)
    tau = bm.where(active, h * langevin / (2 * bm.where(active, speed, 1.0)), 0.0)
```

active 区分近零速度；where 中的1仅用作安全分母。h 除以空间阶次 p。Pe 小于阈值时用级数，其他位置用 coth(Pe)-1/Pe；safe_pe 避免未选择分支中的零除。最后给静止点 tau=0。

#### 完整的试探函数残差（第 32—36 行）

```python
    streamline = bm.einsum("cqd,cqid->cqi", b, grad)
    test = tau[..., None] * streamline
    laplace = bm.einsum("cqidd->cqi", space.hess_basis(bcs)) if space.p == 2 else bm.zeros_like(phi)
    residual = streamline - a[..., None] * laplace
    residual -= bm.einsum("cqd,cqid->cqi", pde.diffusion_gradient(points, t), grad)
```

streamline 是 b·∇φ。test 乘 tau。P2 用 Hessian 迹得到 Δφ，P1 取零；residual 再扣除 aΔφ 和 ∇a·∇φ。这里还没加入反应项，因为反应需要按分量对 i,j 分开。

#### 单元积分函数（第 37—40 行）

```python

    def integrate(trial):
        return bm.einsum("q,c,cqi,cqj->cij", w, measure, test, trial)

```

integrate 接收任意试探项 trial，与 test、权重、单元度量相乘，对求积点求和。结果保留 c,i,j 三轴，表示每个单元的一块局部矩阵。

#### 空间块、时间质量与载荷（第 41—49 行）

```python
    m = pde.components
    blocks = [
        [integrate(C[..., i, j, None] * phi + (residual if i == j else 0.0)) for j in range(m)]
        for i in range(m)
    ]
    mass = integrate(d[..., None] * phi)
    source = pde.source(points, t)
    loads = [bm.einsum("q,c,cqi,cq->ci", w, measure, test, source[..., i]) for i in range(m)]
    return blocks, mass, loads
```

每个分量对都包含 Cij*phi；仅 i==j 的块加空间微分残差。mass 用 d*phi 形成时间修正；loads 每个分量分别乘源项。三个返回值在 solver 中交给 ConstIntegrator。

### cdr/solver.py

#### 导入已有库能力（第 1—22 行）

```python
"""Block Lagrange FEM, SUPG and BDF2 for vector CDR systems."""

from pathlib import Path
from copy import copy
from fealpy.backend import backend_manager as bm
from fealpy.mesh import IntervalMesh, TriangleMesh, TetrahedronMesh
from fealpy.functionspace import LagrangeFESpace, TensorFunctionSpace
from fealpy.fem import (
    BilinearForm,
    LinearForm,
    BlockForm,
    DirichletBC,
    ScalarDiffusionIntegrator,
    ScalarConvectionIntegrator,
    ScalarMassIntegrator,
    ScalarSourceIntegrator,
)
from fealpy.fem.integrator import ConstIntegrator
from fealpy.solver import spsolve
from .pde import at_time
from .stabilization import supg_terms

```

Path 管理文件，copy 用于输出网格副本。三种网格对应空间维数；标量和张量空间管理自由度；Form 与积分器管理装配；DirichletBC 管理边界。ConstIntegrator 来自实际存在的 integrator 模块，spsolve 使用库自己的求解入口。

#### 盒形网格工厂（第 23—30 行）

```python

def box_mesh(pde, n):
    if pde.dim == 1:
        return IntervalMesh.from_interval_domain(pde.domain(), nx=n)
    if pde.dim == 2:
        return TriangleMesh.from_box(pde.domain(), nx=n, ny=n)
    return TetrahedronMesh.from_box(pde.domain(), nx=n, ny=n, nz=n)

```

依 dim 选择线段、三角形或四面体网格，domain 给出边界坐标，每个方向都用相同分段数 n。维数已经在数据类中限制为1、2、3。

#### 求解器初始化（第 31—49 行）

```python

class CDRSolver:
    """Scalar (m=1) or coupled vector solver; NumPy/CPU, affine simplices P1/P2."""

    def __init__(self, pde, mesh, p=1, *, stabilization="supg", q=None, linear_solver=None):
        self._check_backend()
        if p not in (1, 2) or stabilization not in ("supg", "galerkin"):
            raise ValueError("Use P1/P2 and supg/galerkin")
        if mesh.geo_dimension() != pde.dim or mesh.top_dimension() != pde.dim:
            raise ValueError("PDE and full-dimensional mesh must agree")
        self.pde, self.mesh = pde, mesh
        self.scalar_space = LagrangeFESpace(mesh, p=p)
        self.space = TensorFunctionSpace(self.scalar_space, shape=(pde.components, -1))
        self.stabilization, self.q = stabilization, p + 4 if q is None else q
        self.linear_solver = linear_solver
        scalar_boundary = pde.is_dirichlet_boundary(self.scalar_space.interpolation_points())
        self.boundary = bm.tile(scalar_boundary, pde.components)
        self.history = []

```

先核对 NumPy 及已验证的阶次、稳定化方式和网格维数；保存问题与网格；创建标量空间及分量优先的张量空间。q 未给出时取 p+4。边界掩码在全部插值点上求值，因此包括 P2 的边界高阶自由度；tile 按分量复制。history 初始为空。

#### 把环境与系数检查独立出来（第 50—65 行）

```python
    @staticmethod
    def _check_backend():
        if bm.get_current_backend().backend_name != "numpy":
            raise RuntimeError("The validated solver and sparse direct engine require NumPy/CPU")

    def check_coefficients(self, t):
        self._check_backend()
        bcs, _ = self.mesh.quadrature_formula(self.q, "cell").get_quadrature_points_and_weights()
        p = self.mesh.bc_to_point(bcs)
        for field in (self.pde.diffusion_coef, self.pde.capacity_coef):
            value = field(p, t)
            if not bool(bm.all(bm.isfinite(value))) or bool(bm.any(value <= 0)):
                raise ValueError(
                    "Diffusion and capacity must be finite and positive at quadrature points"
                )

```

_check_backend 不修改全局状态，只拒绝不匹配的后端。check_coefficients 在本次时刻的求积点检查 a、d 有限且为正，避免违反当前模型的基本假设；它不是对整个连续区域正性的数学证明。

#### 一次装配的准备（第 66—75 行）

```python
    def assemble(self, t):
        """Assemble all coefficients at the requested time; no stale time cache."""
        self.check_coefficients(t)
        space, m, q = self.scalar_space, self.pde.components, self.q
        dofs = space.cell_to_dof()
        if self.stabilization == "supg":
            blocks, mass, loads = supg_terms(space, self.pde, t, q)
        forms = []
        rhs = []
        massforms = [[None for _ in range(m)] for _ in range(m)]
```

每次调用先检查系数，再取空间、分量数、求积阶和局部映射。SUPG 只计算一次本时刻的所有局部修正；forms 收集空间块，rhs 收集各分量载荷，massforms 先用 None 表示零块。

#### 逐块添加空间积分项（第 76—93 行）

```python
        for i in range(m):
            row = []
            for j in range(m):
                form = BilinearForm(space)
                form.add_integrator(
                    ScalarMassIntegrator(coef=at_time(self.pde.reaction_coef, t, (i, j)), q=q)
                )
                if i == j:
                    form.add_integrator(
                        ScalarDiffusionIntegrator(coef=at_time(self.pde.diffusion_coef, t), q=q)
                    )
                    form.add_integrator(
                        ScalarConvectionIntegrator(coef=at_time(self.pde.convection_coef, t), q=q)
                    )
                if self.stabilization == "supg":
                    form.add_integrator(ConstIntegrator(blocks[i][j], dofs))
                row.append(form)
            forms.append(row)
```

外循环 i 对应方程行，内循环 j 对应未知分量列。所有块加入反应 Cij 的质量积分；对角块再加入扩散和对流。若启用 SUPG，把对应局部修正加入同一个 form。row 与 forms 的嵌套结构正是块矩阵。

#### 载荷、容量和最终装配（第 94—105 行）

```python
            load = LinearForm(space)
            load.add_integrator(
                ScalarSourceIntegrator(source=at_time(self.pde.source, t, (i,)), q=q)
            )
            mf = BilinearForm(space)
            mf.add_integrator(ScalarMassIntegrator(coef=at_time(self.pde.capacity_coef, t), q=q))
            if self.stabilization == "supg":
                load.add_integrator(ConstIntegrator(loads[i], dofs))
                mf.add_integrator(ConstIntegrator(mass, dofs))
            rhs.append(load.assembly())
            massforms[i][i] = mf
        return BlockForm(forms).assembly(), BlockForm(massforms).assembly(), bm.concatenate(rhs)
```

每个分量建立一个 LinearForm，选 f_i。mf 用容量 d 构造质量积分。SUPG 给载荷和质量都加修正。rhs 收集装配向量，massforms 只放对角块；最终返回空间 CSR、质量 CSR 和按分量连接的载荷。

#### 边界与线性求解（第 106—118 行）

```python

    def _solve(self, A, F, t):
        self.A, self.F = DirichletBC(
            self.space, gd=at_time(self.pde.dirichlet, t), threshold=self.boundary
        ).apply(A, F)
        if self.linear_solver is None:
            values = spsolve(self.A.copy(), self.F.copy(), solver="scipy")
        else:
            values = self.linear_solver(self.A.copy(), self.F.copy())
        if values.shape != self.F.shape or not bool(bm.all(bm.isfinite(values))):
            raise ValueError("Linear solver must return a finite vector with the system shape")
        self.uh = self.space.function()
        self.uh[:] = values
```

DirichletBC 接收向量空间、新时间边界函数和全局布尔掩码。默认走 SciPy 稀疏法，也可用回调。复制输入保护用于检查的矩阵与向量。结果形状必须匹配右端且无 NaN/Inf，然后写入有限元函数。

#### 两项独立诊断（第 119—129 行）

```python
        residual = float(bm.linalg.norm(self.A @ values - self.F)) / max(
            1.0, float(bm.linalg.norm(self.F))
        )
        exact_boundary = self.pde.dirichlet(self.scalar_space.interpolation_points(), t).T.reshape(
            -1
        )
        boundary_error = float(
            bm.max(bm.abs(values[self.boundary] - exact_boundary[self.boundary]))
        )
        self.history.append(dict(time=float(t), residual=residual, boundary_error=boundary_error))
        return self.uh
```

残差使用边界处理后的方程组，并除以 max(1,右端范数)。边界数据按分量优先排列，取掩码处的最大绝对差。history 保存真实时间、残差、边界误差；返回仍携带函数空间的 uh。

#### 稳态入口（第 130—135 行）

```python

    def solve_steady(self, t=0.0):
        self.history = []
        A, _, F = self.assemble(t)
        return self._solve(A, F, t)

```

重置记录，在给定 t 装配，忽略质量矩阵，求 A U=F。若制造解原本是含时的，调用者必须使用对应稳态源项，不能期待忽略 u_t 后答案仍相同。

#### 含时初始化（第 136—153 行）

```python
    def solve_time(self, end_time, steps, *, scheme="bdf2", output=None):
        """Uniform-step BDF2 with one backward-Euler start; optional BE comparison."""
        self._check_backend()
        if (
            end_time <= 0
            or not isinstance(steps, int)
            or steps < 1
            or scheme not in ("bdf2", "backward_euler")
        ):
            raise ValueError("Use positive end_time, integer steps and bdf2/backward_euler")
        self.history = []
        points = self.scalar_space.interpolation_points()
        old = self.pde.initial(points).T.reshape(-1).copy()
        initial_boundary = self.pde.dirichlet(points, 0.0).T.reshape(-1)
        old[self.boundary] = initial_boundary[self.boundary]
        previous = None
        dt = end_time / steps
        frames = []
```

检查终止时间、步数及格式。初值在全部插值点上求值，转为分量块排列并复制；边界位置使用 t=0 的 Dirichlet 数据。previous 起步时不存在，dt 为均匀步长；frames 保存将来的 PVD 条目。

#### 保存初始帧（第 154—160 行）

```python
        if output is not None:
            output = Path(output)
            output.mkdir(parents=True, exist_ok=True)
            first = self.space.function()
            first[:] = old
            self.export_vtu(first, output / "solution_0000.vtu")
            frames.append((0.0, "solution_0000.vtu"))
```

只有请求 output 才建目录和输出文件。构造有限元函数 first 承载初值，写 solution_0000.vtu 并记录时间0。输出行为不是求解的必要条件。

#### 主时间循环（第 161—173 行）

```python
        for k in range(1, steps + 1):
            t = k * dt
            A, M, F = self.assemble(t)
            if scheme == "bdf2" and k > 1:
                left, right = 1.5 * M + dt * A, M @ (2 * old - 0.5 * previous) + dt * F
            else:
                left, right = M + dt * A, M @ old + dt * F
            uh = self._solve(left, right, t)
            previous, old = old, uh[:].copy()
            if output is not None:
                name = f"solution_{k:04d}.vtu"
                self.export_vtu(uh, output / name)
                frames.append((t, name))
```

每步 t=k*dt，重新装配全部系数。第二步起 BDF2 使用两层历史，其他情形用 BE。_solve 应用当前边界后求解。历史更新保留旧数组，并复制新系数。若要求输出，每步保存真实帧名和时间。

#### 时间索引文件（第 174—185 行）

```python
        if output is not None:
            lines = [
                '<?xml version="1.0"?>',
                '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"><Collection>',
            ]
            lines += [
                f'<DataSet timestep="{t:.16g}" group="" part="0" file="{name}"/>'
                for t, name in frames
            ]
            lines += ["</Collection></VTKFile>"]
            (output / "solution.pvd").write_text("\n".join(lines), encoding="utf-8")
        return uh
```

XML 声明之后建立 VTK Collection，每个 DataSet 写一个 timestep 和相对文件名；末尾闭合标签。相对文件名使整个目录可搬动。UTF-8 保存文本，返回最终有限元函数。

#### 误差积分（第 186—191 行）

```python

    def errors(self, t=0.0):
        l2 = self.mesh.error(at_time(self.pde.solution, t), self.uh, q=self.q + 2)
        h1 = self.mesh.error(at_time(self.pde.gradient, t), self.uh.grad_value, q=self.q + 2)
        return float(l2), float(h1)

```

mesh.error 分别比较函数值与梯度，积分阶比装配高2级。返回 Python float 方便打印和 JSON 序列化。H1 这一项是梯度半范数。

#### VTU 输出（第 192—206 行）

```python
    def export_vtu(self, uh, filename):
        """Vertex samples; individual components and a VTK vector for m<=3."""
        self._check_backend()
        values = bm.zeros((self.mesh.number_of_nodes(), self.pde.components), dtype=self.mesh.ftype)
        values[self.mesh.cell] = uh(bm.eye(self.pde.dim + 1))
        mesh = copy(self.mesh)
        mesh.nodedata = dict(self.mesh.nodedata)
        for i in range(self.pde.components):
            mesh.nodedata[f"u_{i}"] = values[:, i]
        mesh.nodedata["u_magnitude"] = bm.linalg.norm(values, axis=-1)
        if self.pde.components <= 3:
            vector = bm.zeros((len(values), 3), dtype=values.dtype)
            vector[:, : self.pde.components] = values
            mesh.nodedata["u"] = vector
        mesh.to_vtk(fname=str(filename))
```

分配节点×分量数组，用重心坐标单位矩阵采样单元顶点，再按 cell 连接表写入全局节点。浅复制网格并单独复制节点数据字典，逐分量添加字段和模长；至多三分量时补齐 VTK 三维向量。to_vtk 是 FEALPy 网格方法。
