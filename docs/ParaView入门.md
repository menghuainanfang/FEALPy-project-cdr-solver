# 用 ParaView 查看 CDR 结果

ParaView 是免费的开源科学数据可视化软件。FEALPy 负责计算，VTU 保存网格和数值，ParaView 把这些数据画成云图、切片、曲线或动画。PVD 将多个 VTU 按实际时间组织起来。

## 安装

从 [ParaView 官网](https://www.paraview.org/download/) 下载 Windows 版本。已有的 Python vtk 包用于读写 VTK 数据，不包含 ParaView 桌面界面。ParaView 无需安装进 CDR 的 Python 虚拟环境。

## 打开本项目

1. 启动 ParaView，选择 File → Open。
2. 打开 `results_extension/vtu/solution.pvd`。保持该文件与 21 个 VTU 文件的相对位置不变。
3. 点击左侧 Properties 的 Apply。
4. 颜色字段选择 u，显示方式选择 Surface；若想看到网格，选择 Surface With Edges。
5. 使用 Reset Camera 让网格进入视野。二维平面可沿 Z 轴查看。
6. 使用时间播放控件查看 0 到 1 的变化。每一帧间隔为 0.05。

建议把颜色范围固定为 [0,1]，便于观察峰值随时间降低。每帧自动缩放颜色会让不同振幅看起来相似。初始峰值约为 1，t=1 的峰值约为 0.368。

## 保存和进一步分析

- File → Save Screenshot 保存当前图片。
- File → Save Animation 保存时间动画或图片序列。
- File → Save State 保存当前可视化设置，下次通过 Load State 恢复。
- 二维结果可以用 Warp By Scalar 将 u 值显示成高度；三维结果可用 Slice 查看内部截面。

官方参考：[加载数据](https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html)、[颜色映射](https://docs.paraview.org/en/latest/ReferenceManual/colorMapping.html)、[保存结果](https://docs.paraview.org/en/latest/UsersGuide/savingResults.html)。

## 当前交付

项目已经生成 VTU/PVD，且用 VTK 读取器验证了顶点数、单元数和 u 数据。文件可供 ParaView 打开。当前未完成 ParaView 桌面界面的实际打开验证。

后续可以配置 ParaView 安装、整理云图与切片、固定配色、制作汇报动画，并通过 ParaView 自带 pvpython 编写批量出图脚本。
