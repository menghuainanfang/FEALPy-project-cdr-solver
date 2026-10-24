# 科研工作库

变分移动网格研究与有限元学习实践。

**[资料总索引：论文、学习记录、实验与交付归档](资料总索引.md)**

| 方向 | 内容 | 入口 |
|---|---|---|
| **移动网格** | 毕业论文主线：论文精读、学习交流、算法代码、数值实验与研究产出 | [研究说明](移动网格/README.md) · [研究网页](移动网格/index.html) |
| **学习练手项目** | 独立学习任务、编程练习与科研小项目，目前归档已完成交付的 FEALPy CDR 求解器 | [项目目录](学习练手项目/README.md) · [项目网页](学习练手项目/index.html) |

[在线首页](https://menghuainanfang.github.io/research-workspace/)

```text
movingmesh/
├── README.md
├── index.html
├── 移动网格/
│   ├── README.md
│   ├── index.html
│   ├── docs/                  # 研究主档、计划、实验说明与交接
│   ├── 论文精读/              # 文献索引、精读笔记与学习交流
│   ├── references/            # 参考文献与源码来源
│   ├── experiments/           # 可运行实验
│   ├── src/                   # 算法实现
│   ├── tests/                 # 数学性质与回归验证
│   ├── results/               # 实验数据与交互展示
│   └── tools/                 # 演示打包工具
└── 学习练手项目/
    ├── README.md
    ├── index.html
    └── cdr-solver/            # CDR 已交付归档及配套资料
```

移动网格的资料、代码与结果统一保存在 `移动网格/`。其他任务在 `学习练手项目/` 下各建一个目录，配套项目 README，并加入该目录的项目表。

## 当前移动网格成果

|分类|入口|
|---|---|
|P1：收敛网格优化全流程|[在线演示](https://menghuainanfang.github.io/research-workspace/移动网格/results/mentor_2026-10-23/optimization_flow/index.html) · [实验说明](移动网格/docs/mentor_2026-10-23/optimization_flow.md)|
|P2：长物理时间三重循环|[在线演示](https://menghuainanfang.github.io/research-workspace/移动网格/results/mentor_2026-10-23/p2_v2/index.html) · [实验说明](移动网格/docs/mentor_2026-10-23/three_loop_time_demo_v2.md)|
|代码与数据|[算法](移动网格/src) · [实验](移动网格/experiments) · [结果分类](移动网格/results/README.md)|
|验证与研究资料|[文档分类](移动网格/docs/README.md) · [研究方向与阅读计划](移动网格/docs/research_overview.md)|

## 本地使用

本地工作库：`D:/桌面/movingmesh`。双击 `打开科研工作库.cmd` 浏览首页；运行程序时进入对应项目，按项目 README 配置环境。

论文原件与含原文图件的翻译工程保存在本地的 `移动网格/论文精读/` 中；公开仓库提供来源索引、学习笔记和研究代码。

CDR 已完成交付，按归档项目保留；后续工作集中在移动网格主线和新增学习练手项目。
