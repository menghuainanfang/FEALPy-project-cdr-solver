# 移动网格研究与数值演示

[科研工作库](../README.md) · [演示导航](index.html) · [文档分类](docs/README.md)

当前交付为二维切向移动网格的完整优化流程，以及独立的长物理时间三重循环演示。边上节点沿边移动，四角固定。研究方向与导师任务保存在[研究方向与阅读计划](docs/research_overview.md)，本页按实际完成的程序和证据组织。

## 当前演示

|演示|可以查看什么|数据与说明|
|---|---|---|
|[P1 · 网格优化全流程](results/mentor_2026-10-23/optimization_flow/index.html)|平滑、陡变带、局部尖峰；每轮计算网格全部内迭代、物理映射与外层收敛；能量、对齐、等分布和阈值|[方法与收敛结果](docs/mentor_2026-10-23/optimization_flow.md)|
|[P2 · 物理时间三重循环](results/mentor_2026-10-23/p2_v2/index.html)|静态优化、100步解析移动目标、60步数值Burgers；物理时间n → 耦合k → 内层j|[实验与复现](docs/mentor_2026-10-23/three_loop_time_demo_v2.md) · [三分钟讲稿](docs/mentor_2026-10-23/mentor_walkthrough_v2.md)|

[P1在线演示](https://menghuainanfang.github.io/research-workspace/移动网格/results/mentor_2026-10-23/optimization_flow/index.html) · [P2在线演示](https://menghuainanfang.github.io/research-workspace/移动网格/results/mentor_2026-10-23/p2_v2/index.html)

P1三例实际外迭代分别为3、19、17轮，均达到连续两轮外层判据。P2解析目标实际推进100步至t=2.0，明确未求解PDE；数值Burgers实际推进60步至t=1.2，并独立记录固定网格误差、转移质量差和PDE残量。关键物理步保存全部内层坐标，其他步保存完整诊断与初末坐标。

## 分类入口

|类别|位置|用途|
|---|---|---|
|算法实现|[src](src)|几何、变分能量、BDF监测、网格优化、Burgers求解和时间驱动|
|实验入口|[experiments](experiments)|生成数值记录和演示网页|
|验证|[tests](tests) · [验证结果](results/mentor_2026-10-23/p2_v2/validation.json)|数学性质、收敛状态、时间加密和数据一致性|
|结果与数据|[结果分类](results/README.md)|P1/P2页面、JSON、CSV、逐步压缩记录与历史对照|
|文档|[文档分类](docs/README.md)|方法、字段、验收、研究规划及历史说明|
|文献|[论文索引](论文精读/资料索引.md) · [参考实现](references/README.md)|阅读记录与来源；论文原件留在本地|

## 本地预览与复现

在本项目目录运行：

```powershell
.\.venv\Scripts\python.exe -m tools.preview_demos
```

打开 `http://localhost:8765/demos/index.html`，由P1进入P2。网页回放已保存的真实计算结果；播放速度不改变求解。

```powershell
.\.venv\Scripts\python.exe -m experiments.mesh_optimization_flow
.\.venv\Scripts\python.exe -m experiments.three_loop_time_demo_v2 --scenario analytic --physical-steps 100 --nx 8
.\.venv\Scripts\python.exe -m experiments.three_loop_time_demo_v2 --scenario burgers --physical-steps 60 --nx 8
.\.venv\Scripts\python.exe -m experiments.three_loop_validation_v2
```

Python环境与外部FEALPy源码要求见[环境说明](docs/environment_audit.md)。本轮数值使用FEALPy提交`69faa5f41e3c12633ed714520c3ec512eb1a784b`；FEALPy作为外部依赖保留。

## 验证与限制

本轮49项Python测试、7项JavaScript测试和13个页面检查通过；长时间实验的网格继承差为0，记录网格保持正定向。详细误差和性能见[P2报告](docs/mentor_2026-10-23/three_loop_time_demo_v2.md)。

当前Burgers移动网格误差高于同初始固定网格基线；采用非守恒P1重分区，质量差单独记录。网格收敛、积分器成功和PDE残量达标分别判断。当前结果不构成三维算法或全离散能量稳定性证明。

## 历史实验

[E/F及旧诊断归档](results/README.md#历史对照)保留原始数值证据。旧E、G网址引导到当前P1；Poisson对照作为历史基线保留。[导师任务原文](docs/prompts/README.md)与[阅读计划](docs/research_overview.md)分开保存，完成状态由当前报告给出。
