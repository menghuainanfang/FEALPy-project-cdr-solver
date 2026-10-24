# 实验结果与数据分类

[项目说明](../README.md) · [文档分类](../docs/README.md)

## 当前演示与数值记录

|类别|页面|数据|说明|
|---|---|---|---|
|P1内外迭代|[优化全流程](mentor_2026-10-23/optimization_flow/index.html)|[三种情况及逐接受状态](mentor_2026-10-23/optimization_flow/)|[收敛判据与实测结果](../docs/mentor_2026-10-23/optimization_flow.md)|
|P2物理时间|[三重循环](mentor_2026-10-23/p2_v2/index.html)|[JSON、CSV及逐步gzip记录](mentor_2026-10-23/p2_v2/)|[方法、性能和误差](../docs/mentor_2026-10-23/three_loop_time_demo_v2.md)|
|验证对照|[科研诊断图](mentor_2026-10-23/p2_v2/long_time_diagnostics.png)|[验证摘要](mentor_2026-10-23/p2_v2/validation.json) · [性能](mentor_2026-10-23/p2_v2/performance.csv) · [加密对照](mentor_2026-10-23/p2_v2/refinement.csv)|含固定网格基线与非守恒转移质量差|

P2目录按`情况_网格规模_物理步数_时间步长`命名。每组包含`dataset.json`、`performance.json`、`states.csv`和`steps/`。`burgers_n4_steps10_dt0.02_outer1`为受控预算退出验证；不作为成功演示。`audit/`保存扩展前的三步来源。

## 历史对照

|记录|入口|定位|
|---|---|---|
|E变分网格|[历史页面](archive/e_variational_mesh_2026-10-17.html) · [原始JSON](e_variational_mesh.json)|旧版采样与边界对照；当前演示使用P1|
|G论文教学页|[历史页面](archive/g_huang2015_gallery_2026-10-23.html)|教学资料；当前数值演示使用P1/P2|
|F Poisson|[误差与成本对照](f_poisson_comparison.html) · [方法](../docs/poisson_mesh_comparison.md)|保留已发布的稳态基线|
|阶段诊断|[导师任务结果目录](mentor_2026-10-23/)|P1 v2、计数审计、短时间Burgers与旧案例来源|

旧E/G地址转到当前P1。历史数据保留原有含义，当前新增网格优化仅采用切向边界移动。
