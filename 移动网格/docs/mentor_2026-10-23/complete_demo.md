# 冻结问题诊断记录

2026-10-24；本页保存八个冻结问题的诊断结果。[当前主演示](../../results/mentor_2026-10-23/optimization_flow/index.html)展示完整内外迭代，见[全流程说明](optimization_flow.md)。冻结问题原始JSON保留。

八个切向边界移动案例均为实际BDF计算，网格规模16×16；长时间上限20，短时上限0.5。长短时配对初始网格、冻结度量、边界约束和精度相同，最大允许BDF步长仍沿用cap/10，因此不将跨上限轨迹前缀视为逐位相同。

|案例|最终伪时间s|接受步|完整判据|停止原因|
|---|---:|---:|---|---|
|陡变带，切向边界，长时间|12.276361|120|满足|converged|
|陡变带，切向边界，短时|0.500000|60|未满足|pseudo_time_cap|
|扰动平滑，长时间|4.185149|107|满足|converged|
|扰动平滑，短时|0.500000|86|未满足|pseudo_time_cap|
|初始平衡|0.422222|13|满足|converged|
|较强扰动平滑|4.336230|145|满足|converged|
|平滑达标后继续|8.000000|113|满足，首次s=4.185149|pseudo_time_cap|
|接受步预算12|0.007713|12|未满足|max_accepted_steps|

全部检查的物理单元有向面积为正，最小值约0.000200875。初始平衡案例仍等待固定窗口与连续检查，不在初始残量接近零时立即宣告收敛。达标后继续案例明确关闭早停，展示“网格已满足判据”和“运行最终到达上限”能够同时成立。

案例总览展示最终结果；主时间轴、上一步/下一步、最终结果按钮及五阶段讲解按同一个实际接受状态回放。每个案例的播放帧数为接受步数加1，初态j=0不计为接受步。点击诊断曲线或阶段按钮时先暂停播放，再定位实际网格。短时案例的阶段按钮明确标明“未达到”，窗口不足时显示“等待固定窗口”。

完整库每个接受状态记录残量、阈值、固定窗口、几何、连续计数和积分器统计；窗口就绪前的窗口值为空，表示尚不可计算。SciPy未公开的内部拒绝尝试数仍为null；它不参与收敛判定，不放在主统计表中。网格坐标每个状态保存一次，离线HTML约17.3 MiB。

重建入口：

```powershell
cd D:\桌面\movingmesh\移动网格
.\.venv\Scripts\python.exe -m experiments.p1_complete_demo
.\.venv\Scripts\python.exe -m tools.render_mesh_viewers
node --test tests/mesh_diagnostics_view.test.cjs
.\.venv\Scripts\python.exe tools/check_mentor_pages.py
```

数据与逐案例摘要位于演示同目录的`e_variational_mesh.json`、`summary.json`。原版E/H数值JSON不覆盖；历史HTML仅更新说明和入口。H默认入口采用已有v2完整监控轨迹的两步PDE演示。新实验与当前页面仅保留切向边界移动，边上节点法向受约束、四角固定。

## 验收与限制

五项Node测试通过，其中完整案例测试逐状态核对网格显示时刻、接受步编号、能量和诊断时刻一致，检查必要数值字段、几何合法性，以及各案例的实际停止状态；页面数据测试确认E/F只展示切向移动案例。十一页链接、内嵌数据、脚本语法和元素ID静态检查通过。

本地HTTP预览入口可用，导航页响应为HTTP 200；浏览器实际播放、点击和桌面/窄屏渲染尚未验收。数值、数据一致性和静态检查不替代浏览器验收。上述数值平衡判据不证明能量极小或全离散稳定性。

## 本地HTTP预览

在项目目录运行`.\.venv\Scripts\python.exe -m tools.preview_demos`，保持该进程运行，打开[实验导航](http://localhost:8765/demos/index.html)或[完整收敛案例](http://localhost:8765/results/mentor_2026-10-23/p1_complete/e_variational_mesh.html)。Ctrl+C停止本次预览；`--port`可指定端口。

服务仅绑定本机127.0.0.1。[官方浏览器说明](https://learn.chatgpt.com/docs/browser)支持启动本地开发服务后预览localhost地址。若浏览器显示权限拦截，在Settings → Browser → Agent permissions中为该站点配置浏览权限；组织管理的限制需由管理员处理。此入口使用HTTP，不要求浏览器直接读取file协议。
