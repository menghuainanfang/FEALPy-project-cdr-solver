# 科研工作库

以 **高维鲁棒的能量稳定变分移动网格算法研究及其应用** 为主线，集中管理论文阅读、学习交流、数值实验和研究产出。其他科研任务与练习放在 `projects/`，每个项目独立成目录，并配套 README。

[科研导航网页](https://menghuainanfang.github.io/FEALPy-project-cdr-solver/) · [主线项目说明](docs/moving_mesh_project.md) · [论文与学习资料](论文精读/README.md) · [独立项目目录](projects/README.md) · [CDR 在线汇报](https://menghuainanfang.github.io/FEALPy-project-cdr-solver/projects/cdr-solver/05_汇报/index.html)

## 主线：变分移动网格

| 入口 | 内容 |
|---|---|
| [研究主档](docs/高维鲁棒能量稳定变分移动网格_研究主档.md) | 研究问题、目标与本机基线 |
| [研究计划](docs/research_plan.md) | 阶段目标、时间安排和待导师确认事项 |
| [论文资料索引](论文精读/资料索引.md) | 文献、阅读路线与当前学习断点 |
| [首篇论文精读总览](论文精读/briefs/几何离散与变分网格生成_本科生精读总览.md) | Huang–Kamenski 几何离散方法的阅读入口 |
| [学习交流记录](<论文精读/notes/first paper/首篇论文_学习留痕合并.md>) | 已有推导、问答、回顾与后续续学记录 |
| [项目交接](docs/session_handoff.md) | 实验进展、证据与下一步 |
| [源码地图](docs/fealpy_code_map.md) | FEALPy 源码与自有程序的对应关系 |

## 独立项目与练习

| 项目 | 类型 | 内容与入口 |
|---|---|---|
| [FEALPy CDR 求解器](projects/cdr-solver/README.md) | 独立科研任务 | 1—3 维对流—扩散—反应有限元求解；程序、数学讲义、FEALPy 导读、学习路线、汇报和开发阶段档案 |

新增项目或练习时，在 `projects/` 下新建一个目录，在 README 中说明目标、材料、运行或学习方法、产出与当前进展，并在本表和 [项目目录](projects/README.md) 中增加入口。

## 文件布局

```text
moving-mesh-thesis/
├── README.md                 # 科研工作库总入口
├── index.html                # 科研导航网页
├── docs/                     # 主线研究规划、实验说明、环境与交接
├── 论文精读/                 # 主线论文、精读笔记、学习交流与本地译稿
│   ├── README.md
│   ├── 资料索引.md
│   ├── papers/               # 本地论文原件
│   ├── briefs/               # 论文总览
│   ├── notes/                # 持续学习与交流记录
│   └── exports/              # 本地阅读交付物
├── references/               # 文献来源与公开代码索引
├── experiments/              # 主线 A—F 数值实验
├── src/                      # 主线可复用算法
├── tests/                    # 主线数学性质与回归验证
├── results/                  # 主线实验结果与网页演示
├── tools/                    # 主线演示打包工具
└── projects/                 # 独立科研任务与练习
    ├── README.md
    └── cdr-solver/            # CDR 程序、讲义、汇报与开发阶段
```

## 实验与展示

- [变分网格映射动画](results/e_variational_mesh.html)：比较固定边界与切向边界，观察能量、单元有向面积、对齐与等分布指标；[实验说明](docs/variational_mesh_visual.md)。
- [Poisson 网格精度对比](results/f_poisson_comparison.html)：三个制造解、三档规模、五种物理网格的误差与成本；[实验说明](docs/poisson_mesh_comparison.md)。
- [CDR 求解器汇报](projects/cdr-solver/05_汇报/index.html)：多分量含时求解、稳定化和收敛数据。

在本机根目录运行主线实验：

```powershell
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual
.\.venv\Scripts\python.exe -m experiments.f_poisson_mesh_comparison
```

其他电脑先按 [环境说明](docs/environment_audit.md) 准备 FEALPy 与 Python 依赖。CDR 使用自己的 [安装与运行说明](projects/cdr-solver/01_程序/README.md)。两个项目分别管理依赖。

## 保存与公开范围

本地工作目录为 `D:/桌面/moving-mesh-thesis`。公开科研工作库沿用 [FEALPy-project-cdr-solver](https://github.com/menghuainanfang/FEALPy-project-cdr-solver)，移动网格原仓库的 Git 历史与远端配置保留。

公开内容包括研究代码、研究规划、学习笔记、文献来源、独立项目和选定展示结果。本地论文原件、含原文图件的翻译工程、虚拟环境、临时文件与分发压缩包保存在本机；文献检索入口见 [参考资料](references/README.md)。资料索引中的本地论文与译稿链接供本机阅读使用。

## 当前研究边界

“二维做深、三维验证”和“固定度量下网格能量稳定优先”仍为待导师确认的工作假设。当前主线已有几何、等分布、Poisson 与二维移动网格基线；尚无三维、Burgers 或全离散能量稳定结论。详细实验条件与已记录结果见 [主线项目说明](docs/moving_mesh_project.md) 和各实验报告。
