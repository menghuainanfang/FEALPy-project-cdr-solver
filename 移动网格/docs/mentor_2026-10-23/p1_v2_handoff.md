# P1 v2：结果、重跑与交接

2026-10-24；分支 `mentor/p1-convergence-v2-2026-10-24`。

## 结果入口（本机数值已验证）

- [E共用判据面板](../../results/mentor_2026-10-23/p1_v2/e/e_variational_mesh.html)：提前达标、恰在上限达标、到上限未达标三种真实场景。
- [H共用面板验证](../../results/mentor_2026-10-23/p1_v2/h_preview/h_burgers.html)：现有Burgers求解器的两个物理步，物理时间/耦合轮次/网格接受步分别显示。
- [静态输入](../../results/mentor_2026-10-23/p1_v2/static/h_burgers.html)：物理步n=0，仅优化网格；PDE误差为null。
- [计数审计](convergence_counts_audit_v2.md)、[字段词典](mesh_field_dictionary_v2.md)、[P2可复用的数据契约](convergence_trace_contract_v2.md)。
- [实际接受步科研图](../../results/mentor_2026-10-23/p1_v2/e/convergence_diagnostics.png)：残量阈值、首次完整达标、分epoch能量与真实步长。

## E：为何在这里停止

|场景|保存帧|BDF接受步|nfev|首次单次全条件s|首次完整连续判据s|最终s|criterion_met|cap_reached|停止原因|
|---|---:|---:|---:|---:|---:|---:|---|---|---|
|smoothing|52|107|263|3.234146|4.185149|4.185149|True|False|converged|
|adaptive_fixed|61|103|233|4.533893|5.000000|5.000000|True|True|converged|
|adaptive|61|104|247|未达到|未达到|5.000000|False|True|pseudo_time_cap|

切向自适应末态仍有三个真实未达条件：r_RMS=0.000513114 > 5.56072e−5；窗口相对位移0.00115829 > 1e−5；窗口能量变化率5.36672e−7 > 1e−7，连续计数0/3。
平滑场景的残量先在约0.808达标，但辅助条件和连续计数完成得更晚；“残量已低于线”与“完整停止条件已成立”分开显示。
固定边界场景在s=5同时达标并到上限，记录两个true，不能写为提前收敛。
旧H t=0.02/k=1 的111步在s=3才完成旧连续判据；v2固定窗口在同一冻结问题下于s=2.886990达标。辅助判据版本改变，不能将该差异表述为BDF或泛函算法加速。

## 66组有界对照

36组核心：n=8/16 × 三场景 × 上限0.5/3/5 × 两档BDF精度；12组轻/强扰动；12组额外21/121帧采样对照；6组达标后继续到5。
核心物理初值统一为seed=42、内部扰动±0.2/n；额外扰动为0.05/n、0.3/n。自适应度量在实际初始物理单元中心冻结。
此矩阵的自适应场景也有初始扰动；E页面自适应场景保持原版均匀初始网格，因此两个表不混作同一初值。
核心36组中12组converged、24组到上限；额外扰动12组中5组converged、7组到上限。
驱动合计墙钟约91.08秒；不含JSON/CSV写入与最终科研绘图。所有接受检查的最小有向面积为正（全组最小约0.000112529）。

较严精度rtol=1e−6、atol=1e−8的核心结果：

|n|场景|上限|首次完整达标s|最终s|停止|到上限|保存帧|接受步|nfev|墙钟s|
|---:|---|---:|---:|---:|---|---|---:|---:|---:|---:|
|8|smoothing|0.5|未达到|0.500000|pseudo_time_cap|True|61|86|225|0.363|
|8|smoothing|3.0|1.808398|1.808398|converged|False|38|99|251|0.326|
|8|smoothing|5.0|1.808398|1.808398|converged|False|23|99|251|0.461|
|8|adaptive_fixed|0.5|未达到|0.500000|pseudo_time_cap|True|61|94|251|1.129|
|8|adaptive_fixed|3.0|2.113647|2.113647|converged|False|44|118|299|1.149|
|8|adaptive_fixed|5.0|2.113647|2.113647|converged|False|27|118|299|0.818|
|8|adaptive|0.5|未达到|0.500000|pseudo_time_cap|True|61|94|255|0.872|
|8|adaptive|3.0|未达到|3.000000|pseudo_time_cap|True|61|128|323|1.123|
|8|adaptive|5.0|3.692513|3.692513|converged|False|46|131|329|1.268|
|16|smoothing|0.5|未达到|0.500000|pseudo_time_cap|True|61|86|221|2.649|
|16|smoothing|3.0|未达到|3.000000|pseudo_time_cap|True|61|105|259|3.235|
|16|smoothing|5.0|4.185149|4.185149|converged|False|52|107|263|2.434|
|16|adaptive_fixed|0.5|未达到|0.500000|pseudo_time_cap|True|61|98|286|2.778|
|16|adaptive_fixed|3.0|未达到|3.000000|pseudo_time_cap|True|61|130|354|3.759|
|16|adaptive_fixed|5.0|未达到|5.000000|pseudo_time_cap|True|61|140|374|1.315|
|16|adaptive|0.5|未达到|0.500000|pseudo_time_cap|True|61|97|289|0.989|
|16|adaptive|3.0|未达到|3.000000|pseudo_time_cap|True|61|131|366|1.219|
|16|adaptive|5.0|未达到|5.000000|pseudo_time_cap|True|61|141|386|1.333|

完整[JSON](../../results/mentor_2026-10-23/p1_v2/matrix/summary.json)/[CSV](../../results/mentor_2026-10-23/p1_v2/matrix/summary.csv)同时包含阈值、初末残量/能量、njev/nlu、最小面积、步长极值/中位值及明确的mesh_status/solver_status。
每次运行的固定输入、metric、终点、输出采样与逐接受诊断另存JSON/CSV。

## 采样、步长和继续积分

21/61/121帧：12个与61帧的配对全部具有相同接受时刻序列、终点（最大差0）、停止原因与首次达标时刻。该结论是本次实际对照结果。
6组同上限、同max_step、同初值且关闭早停的旧/新驱动对照，终点最大差均为0，nfev与接受步完全一致。
不同上限仍沿用原版max_step=cap/10，所以跨上限比较也改变了最大允许内部步长；不把跨上限的前缀差异当作方程变化。

以cap=3、61个请求采样、较严精度为例：

|n|场景|接受Δs最小|中位值|最大|名义输出采样间隔|
|---:|---|---:|---:|---:|---:|
|8|smoothing|1.4817e-06|0.00183989|0.162179|0.05|
|8|adaptive_fixed|1.68706e-06|0.00376989|0.147288|0.05|
|8|adaptive|1.68084e-06|0.00462711|0.135837|0.05|
|16|smoothing|1.53564e-06|0.00253116|0.274765|0.05|
|16|adaptive_fixed|1.18853e-06|0.00468603|0.150962|0.05|
|16|adaptive|1.18272e-06|0.00465455|0.138745|0.05|

提前停止时截断展示采样并补实际终点，不强造请求帧数；BDF步长来自误差控制。

|n|继续积分场景|首次完整达标s|首次后继续到5的能量减少|后续失去完整判据次数|终点网格判据|
|---:|---|---:|---:|---:|---|
|8|smoothing|1.808398|6.23945e-13|0|True|
|8|adaptive_fixed|2.113647|3.33955e-13|0|True|
|8|adaptive|3.692513|2.14406e-12|0|True|
|16|smoothing|4.185149|1.10956e-12|0|True|
|16|adaptive_fixed|未达到|不适用|0|False|
|16|adaptive|未达到|不适用|0|False|

继续积分6组均按显式诊断策略运行到cap，4组曾完整达标且到终点未失效，能量减少约1e−12；另外两组始终未完整达标。只称这些冻结问题的观察。

## 工程验收与重跑

Python：39项通过（4.21s），包括原有30项与本次新增9项。新增覆盖固定窗口、小步长卡顿、初始零运动、cap与criterion同时发生、几何/积分失败、边界自由度、能量缩放及epoch重置。
Node：3项通过，核对旧E不将velocity伪装成残量、旧H的111接受步与8采样帧分开、未收敛案例不生成虚假讲解事件。
10个页面静态检查通过：本地链接、内嵌JSON、JavaScript语法、面板与被引用元素ID、模板占位符。
科研PNG已目视检查。包内依赖清单已加入监控器、序列化、词典与共用面板脚本。

```powershell
cd D:\桌面\movingmesh\移动网格
.\.venv\Scripts\python.exe -m experiments.p1_counts_audit_v2
.\.venv\Scripts\python.exe -m experiments.p1_convergence_matrix_v2
.\.venv\Scripts\python.exe -m experiments.e_variational_mesh_visual --duration 5 --output-dir results/mentor_2026-10-23/p1_v2/e
.\.venv\Scripts\python.exe -m experiments.h_burgers_physical_time --steps 2 --output-dir results/mentor_2026-10-23/p1_v2/h_preview
.\.venv\Scripts\python.exe -m experiments.h_burgers_physical_time --static --output-dir results/mentor_2026-10-23/p1_v2/static
.\.venv\Scripts\python.exe -m tools.export_mesh_field_dictionary
.\.venv\Scripts\python.exe -m tools.render_mesh_viewers
.\.venv\Scripts\python.exe -m experiments.p1_convergence_figures_v2
.\.venv\Scripts\python.exe -m pytest -q
node --test tests/mesh_diagnostics_view.test.cjs
.\.venv\Scripts\python.exe tools/check_mentor_pages.py
```

Node未在PATH时使用本机已核验的 `C:/Users/10251/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node.exe`；静态检查可通过`--node`指定。
阶段提交：cb46086为旧输出与计数审计，a52371b为独立监控器、66组对照与契约；共用面板和完整交接在随后提交。
原版E/H数值JSON保留；HTML采用可读取历史JSON的共用解释面板。P2可直接复用MeshConvergenceMonitor与convergence_trace，PDE状态与外层一致性仍分别管理。

## 限制

- 早期“上限0.5对应57”的日志尚未定位；当前H的57来自不同物理时刻、不同上限，不能据此认定历史截图来源。
- SciPy公开API不暴露准确内部拒绝尝试数，记null；不从RHS计数推断。
- 达到本例数值平衡判据不证明能量局部/全局极小、全离散稳定性或任意网格非退化；不同metric epoch的能量不连成单调性结论。
- 本机浏览器安全策略禁止file协议；页面交互播放、鼠标点击、窄屏实际渲染和页面截图尚待浏览器验收。静态检查、纯函数契约测试与科研PNG检查不替代该验收。
