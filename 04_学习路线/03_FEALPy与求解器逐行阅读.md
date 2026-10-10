# FEALPy 与 CDR 求解器联合逐行阅读路线

这是一份后续教学安排，不是已完成学习记录。已有学习断点仍以 [学习顺序](01_学习顺序.md) 为准。

## 每次读一条完整调用链

先说明这一小段要解决的数学问题，运行最小例子，再逐行说明变量、形状、调用对象与返回值。遇到 FEALPy 调用时，先标注“库原生”，再进入实际安装版本的对应源码；回到 cdr 文件后标注“项目实现”。不要按文件长度机械地一次读完。

每一小段回答：为什么需要它？每个变量是什么？输入/输出 shape 是什么？对应哪个数学项？执行后改变了什么？在看懂的范围做一次小改动和解释。

## 阅读顺序

| 小主题 | 项目源码 | FEALPy 原生源码/对象 | 小实验 |
|---|---|---|---|
| 导入与对象创建 | example_minimal.py、cdr/__init__.py | backend_manager | 打印对象来源，不把 import 当求解 |
| 数据与模型 | cdr/pde.py 的 CDRModel、SymbolicCDR | cartesian | 求三个坐标点上的源项 |
| 网格 | cdr/solver.py 的 box_mesh | IntervalMesh、TriangleMesh、TetrahedronMesh | 对照分段数、单元数、节点数 |
| 标量与物质轴 | CDRSolver.__init__ | LagrangeFESpace、TensorFunctionSpace | 解释 81×3=243 以及长向量排列 |
| 一个扩散矩阵 | assemble 中对角项 | ScalarDiffusionIntegrator、BilinearForm | 对照三节点手算 |
| 反应与分块 | assemble 的 i,j 循环 | ScalarMassIntegrator、BlockForm | 从 2×2 的 C 到 4×4 的反应块 |
| 载荷与边界 | at_time、_solve | LinearForm、ScalarSourceIntegrator、DirichletBC | 把零边界换成已知非零值 |
| 线性求解与函数 | _solve | spsolve、space.function | 区分系数数组与有限元函数 |
| 时间推进 | solve_time | 稀疏矩阵乘法 | 手算第一步，追踪 old/previous |
| 稳定化 | cdr/stabilization.py | 求积、basis、grad_basis、ConstIntegrator | 解释完整残差每项的来源 |
| 验证与输出 | errors、export_vtu、verify_vector.py | mesh.error、网格导出 | 比较残差、误差与图像 |

FEALPy 源码位置通过 inspect.getfile 查询实际环境，不固定为某台电脑的绝对路径。项目方法定义在 cdr/，库方法定义在安装的 fealpy/；两边使用相似风格并不意味着来源相同。

## 每次学习的留痕

记录已读函数与版本、学习者原话、解释、运行结果、小练习与尚未理解的问题；以实际交互确认理解，不预先标记完成。一次只沿当前主题展开，遇到问题停在对应语句，不用整段抽象公式代替解释。

新生先完成两份主讲义的入门例子；已有基础的学习者从真实断点继续。历史一维版本适合比较设计演变，但后续维护以当前 cdr 包为主。
