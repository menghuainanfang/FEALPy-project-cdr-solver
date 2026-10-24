# 本地 Codex 执行提示词：收敛审计、迭代统计纠正、专用网页与物理时间三重循环

> 导师专项任务：2026-10-23 至 2026-10-30；P1 收敛诊断、P2 三重循环演示。
> 这是执行任务说明而不是实验完成报告。第一轮必须先诊断旧程序，不要先重写前端。

你正在 Windows 本地科研工作区工作。主仓库为 research-workspace/移动网格，GitHub： https://github.com/menghuainanfang/research-workspace/tree/main/移动网格 。已有本地路径可能是 D:\桌面\movingmesh\移动网格，FEALPy 可能位于 C:\Users\10251\fealpy。绝不能假设这些路径、Python 版本、依赖或源码提交仍然正确，首先实际探测并记录。

【固定前提与授权范围】
毕业研究的长期核心是导师团队 Wang–Huang–Wei 2026 迹—对数泛函的三维四面体实现；但这周仅改良二维演示并探索物理时间耦合。保留已有数值实现的科学口径和实验记录。先阅读 README、docs/environment_audit.md、docs/variational_mesh_visual.md、docs/experiment_protocol.md、docs/fealpy_code_map.md、tests 和 session_handoff。检查工作树，保护未提交更改；使用功能分支/单独提交，不覆盖用户正在改的代码、不强推、不删除历史结果、不自动合并到主分支。如环境缺依赖，只在隔离 venv 安装，记录安装命令。不要修改 FEALPy 第三方源码作为默认方案。

【已通过 GitHub 预读确认的事实，仍须本地复核】
- experiments/e_variational_mesh_visual.py：默认 --duration 0.5、--frames 61；以 np.linspace(0,duration,frames) 构造展示采样；这不是61个BDF接受步。
- src/fealpy_mesh_evolution.py：solve_ivp(..., method="BDF", t_eval=output_times, jac_sparsity=..., rtol=1e-6, atol=1e-8, max_step=duration/10)。当它成功走到伪时间终点，直接回填 stop_reason="pseudo_time_end"。没有网格收敛停止条件。
- history.iteration 是输出帧索引，history.step 是两帧采样时刻之差，gradient_norm 当前实际上是 FEALPy 生成的 velocity 的欧氏范数，而不是独立求出的真实 ∇I 范数；backtracks=0 是占位字段。页面中的 BDF 右端调用次数来自 nfev，未报告真实接受步数。
- 导师观察“61次和57次均对应0.5秒伪时间”。57若来自 nfev，则是 RHS 调用次数而非接受步数，先从输出和日志证明究竟是什么；不得事先宣称已核验本地57这个数字。0.5 是无量纲/规定单位的网格伪时间参数，不能写为真实运行耗时或物理时间。
- experiments/mesh_viewer.html 是科学网格演示，物理网格/计算网格并排是必要对照；移动网格/index.html 目前混合了“研究与学习”和“交互实验”两栏。请用实际页面结构确认导师要分开的正是“研究任务主档与演示”，不要误删计算/物理双网格对比。
- 上游 FEALPy develop 源码线索（与本地不一定一致）：
  https://github.com/weihuayi/fealpy/blob/develop/example/meshopt/mmesh/scalar_burgers_mm_solver.py
  https://github.com/weihuayi/fealpy/blob/develop/example/meshopt/mmesh/rayleightaylor_mm_solver.py
  https://github.com/weihuayi/fealpy/blob/develop/fealpy/cfd/chan_hilliard_navier_stokes_fem_model.py
  https://github.com/weihuayi/fealpy/blob/develop/example/fem/allencahn_lfem_example.py

【任务 A：原版程序的证据链与数值停止条件——必须优先做】
A1. 不改变代码先在当前环境跑旧默认的 smoothing / adaptive_fixed / adaptive 三组场景，保存 git SHA、FEALPy SHA/dirtiness、解释器、依赖、命令和原始 JSON。统计表必须独立列：输出帧数/t_eval、终点伪时间、真实被接受的BDF步数、接受步长序列、nfev/njev/nlu、耗时、速度范数和能量变化、最小有向面积、求解器失败/终止理由；给出61、57的明确来源及实际运行输出，不能推测当实测。
A2. 若 solve_ivp 的 t_eval 不包含内部步长历史，可写独立的诊断驱动，使用 scipy.integrate.BDF 类 .step()、.t_old、.step_size 和计数属性采集接受步，必要时用 dense_output 在统一伪时间抽样；保留现有 solve_ivp 基线及原结果以作对照。避免把 RHS 函数调用、Jacobian 评估和线性系统LU次数混成“迭代”。
A3. 加入可配置收敛逻辑。优先使用与 FEALPy 数值速度一致的自由度投影速度/残量范数，另测网格位移相对局部网格尺度、固定能量的相对变化；阈值应量纲/尺度一致，包含绝对与相对判据、最小伪时间或连续成功次数等防止瞬时误判机制。检测 min signed area 和固定/切向边界约束。停止原因至少明确：converged、pseudo_time_cap、solver_failure、invalid_mesh、max_internal_steps（若适用）。指定容差不可通过仅提高阈值伪造收敛；“速度≈0”必须考虑时间尺度 τ、平衡量、自由度和实际残量定义。
A4. 输出 accepted BDF step traces、网格收敛曲线、冻结度量能量曲线与关键指标 JSON/CSV。界面将 iteration 改称 output_frame_index；step 改称 sample_dt；gradient_norm 若仍表示速度，必须改称 velocity_norm 并注明定义；backtracks=0 无真实回溯就删除或注释为不适用。绝不能写“已能量稳定”而仅靠存帧能量下降。
A5. 做至少以下对照（可缩小网格但保留可重现命令）：网格 n=8/16；原版固定 t_cap=0.5 与新收敛停止；初始轻微扰动和强扰动；两档求解器精度。三场景优先全部对照。记录是否达到阈值、错误步数、CPU墙钟、最小面积、插值误差；给导师一段可复述的“为什么61和57不同、0.5是什么、何时真正停止”。

【任务 B：让专门网页只做专门事情】
B1. 先检查 index.html、experiments/mesh_viewer.html、results/*.html 及研究主档引用。区分：研究导航/主档网页、单独实验导航页、E网格演示、F Poisson对照、G经典复现、后续物理时间演示。研究主档保留在文档页和独立导航入口；面向导师的演示页不要在同一屏强行塞进项目任务看板。
B2. 尽量用小改动新增独立 demos/index.html 或等价专门实验入口，维持旧URL可用（或清晰重定向/兼容链接），修正相对路径，给新旧页面加清楚的返回/导航入口。保留演示页“计算/物理网格”并排，因为这是数学比较，不是混合“任务/演示”。无需新框架或大规模视觉重构。检查静态离线可打开、交互播放可用、路径无404，并提供截图/页面检查说明。

【任务 C：真实物理时间驱动的三重循环演示】
C0. 先定位本机 FEALPy 移动 Burgers 或 CH–NS/Rayleigh–Taylor 示例的接口、更新物理解、网格更新、重映射、非线性求解和时间推进策略。绘制依赖图，标记上游已有能力、兼容性和可直接复用的代码，不得把查到源码当成跑通。
C1. 正确建模三个层次：
  外层 n：真实物理时间 t_n→t_{n+1}，固定物理时间步 Δt_phys（或自适应），推进 u 或 φ；
  中层 k：在该物理时刻/时间步内，交替更新物理解、度量 M(u)、重建物理网格，检测外循环一致性，只有需要时迭代；
  内层 j：在一次冻结旧物理网格和 M 的前提下，推进计算坐标的网格伪时间 s，按任务A的真实收敛判据或上限停止。
  可按具体算法调整实际操作顺序，但必须明确依赖与冻结量；绝不能把帧序号当三重循环次数。
C2. 首选构建小型可运行的二维 Burgers 随物理时间演化实验，复用本机求解框架；可先降低网格/时间档次，使用可靠初边值、PDE残量或制造解误差验证。每次物理时步更新场 u，产生新 metric，并从上一时刻网格继续优化而非从初始均匀网格重来。根据离散方式处理 mesh velocity、场的转移（ALE或保守/非守恒重映射），说明质量守恒或几何守恒条件；不要直接把新网格上的旧节点值当成精确PDE更新。
C3. 若本周真实移动网格PDE耦合暂未跑通，可先用明确标为“解析随时间目标函数驱动的网格更新演示”的替代版本展示三重循环架构；必须同时保留可运行的固定网格非线性 PDE 基线及真实耦合尚未完成的声明。禁止将时间变化的解析函数动画伪装成数值求得的物理解。
C4. 网页新增专门场景，分别显示 physical_time、outer_iteration、mesh_inner_accepted_steps、pseudo_time、PDE求解状态、当步 M 版本、是否达到网格和外循环收敛、最小有向面积、误差或能量。至少能按物理时间逐帧播放，且点开某物理时刻能观察其网格内迭代变化。两个时间轴必须标出不同含义；旧网格从 t_n 延续到 t_{n+1}。
C5. 对 2–4 个物理时间步保存原始 JSON 和图片，说明网格、metric、PDE 解确实改变；对静态时间输入退化为单次网格优化、对强扰动不会静默翻转。不存在的 PDE 误差不能编造。
C6. 对比导师论文方程/现有 FEALPy 算法的正确边界，实测仅称“本例验证”，不得扩大为全离散稳定性或三维结论。

【完成标准与工作输出】
提交 docs/mentor_2026-10-23/diagnosis.md（61/57/0.5的事实与数值解释）、源码变更、测试（区分原测试和新增测试）、三重循环流程图、面向导师的独立演示网页、最小重跑命令、关键性能表、docs/session_handoff_mentor_2026-10-30.md。完成至少测试：采样帧数不影响内部步统计；收敛前停止或明确到上限；不正定向必失败；二维边界约束保持；物理时间不与伪时间混淆；每个真实时步更新正确的场/metric；静态路径不退化。对需要较长运算的任务先做最小网格，记录成本和阻塞，避免无限跑。
每项证据附代码文件和函数位置；所有结论标记“本机已验证/源码可见/待验证”。对于无法在一周内完成的部分，提交最小可运行中间结果和精确阻塞原因，而非虚报完成。

【开始的第一个动作】
不要直接重写代码。先读仓库工作树与 src/fealpy_mesh_evolution.py、experiments/e_variational_mesh_visual.py、experiments/mesh_viewer.html、index.html、旧运行输出；运行一个最小默认实验并出一张“伪时间输出帧、BDF接受步、nfev及停止条件对照表”。之后按 A→B→C 实施，并在每个阶段保留可回滚提交。
