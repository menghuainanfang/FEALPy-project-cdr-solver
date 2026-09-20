# 从微分方程到计算结果：向量含时 CDR 数学讲义

这份讲义面向学过微积分、线性代数，尚未系统学习有限元的数学本科生。我们的目标是把一个偏微分方程变成计算机能求解的线性方程组，并且知道结果为什么可信。建议顺着第 1—8 节完成一次纸笔推导，再结合[代码讲义](CDR代码与FEALPy源码学习讲义.md)运行程序。阅读完材料不等于掌握；文末的小题可用于自检。

## 1. 我们究竟在求什么

### 1.1 坐标、时间与未知量

空间维数记为 $D$，在程序中叫 `dim`；未知量的分量数记为 $m$，叫 `components`。两者是不同的概念。

- 一根细杆可以用一个坐标 $x$ 表示，$D=1$。如果同时研究两种物质的浓度，未知量有两个分量，$m=2$。
- 一个立方体需要 $(x,y,z)$ 三个坐标，$D=3$。若只研究温度，仍然可以只有一个未知量，$m=1$。

向量函数写成 $\boldsymbol u(\boldsymbol x,t)=(u_1,\ldots,u_m)^T$。这里的向量可以是一组物理量，并不一定是空间中的速度箭头。我们的方程是

$$
d\frac{\partial u_i}{\partial t}
-\nabla\cdot(a\nabla u_i)
+\boldsymbol b\cdot\nabla u_i
+\sum_{j=1}^m C_{ij}u_j=f_i,\qquad i=1,\ldots,m.
\tag{1}
$$

每个系数都可以依赖 $\boldsymbol x,t$。$a$ 和 $d$ 为各分量共享的正标量；$\boldsymbol b$ 有 $D$ 个空间分量；$C$ 是 $m\times m$ 的矩阵。$a$ 描述扩散，$\boldsymbol b$ 描述输运，$C$ 描述局部反应及分量交换，$f_i$ 是外部源项。

本项目采用全边界 Dirichlet 条件和初值：

$$u_i|_{\partial\Omega}=g_i(\boldsymbol x,t),\qquad u_i(\boldsymbol x,0)=u_{0,i}(\boldsymbol x).$$

Dirichlet 的含义只是“边界上直接给出未知量的值”。它不同于给出法向导数或流量。稳态问题删去 $d\partial_tu_i$，其余结构保留。

### 1.2 梯度和散度并不神秘

二维中

$$\nabla u=(u_x,u_y)^T,\qquad
\nabla\cdot(a\nabla u)=\partial_x(au_x)+\partial_y(au_y).$$

按乘法求导法则展开：

$$-\nabla\cdot(a\nabla u)=-a\Delta u-\nabla a\cdot\nabla u,\qquad
\Delta u=\sum_{k=1}^D u_{x_kx_k}.
\tag{2}$$

所以变系数扩散不能只写成 $-a\Delta u$。例如一维 $a=1+x^2,u=x^2$，有

$$-(au')'=-(2x+2x^3)'=-2-6x^2,$$

而 $-au''=-2-2x^2$，少了 $-4x^2$。程序的 `diffusion_gradient` 正是提供式 (2) 中的 $\nabla a$。

对流项是 $\boldsymbol b\cdot\nabla u$，不是 $\nabla\cdot(\boldsymbol b u)$。二者相差 $(\nabla\cdot\boldsymbol b)u$。本项目直接按前者实现。

### 1.3 时间前面的 d 是系数

式 (1) 中 $d$ 是容量系数；偏导符号是 $\partial$。即使 $d=d(t)$，我们仍采用 $d(t)u_t$。乘积求导则为

$$\partial_t(du)=d u_t+d_tu.$$

例如 $d=1+t,u=e^t$，则 $du_t=(1+t)e^t$，而 $\partial_t(du)=(2+t)e^t$。它们是两个不同模型。若物理题目指定后一种，必须保留多出的 $d_tu$，不能只改函数名。

### 1.4 一个真正耦合的例子

取

$$C=\begin{pmatrix}2&1\\1&5\end{pmatrix},\qquad
\boldsymbol u=\begin{pmatrix}1+x^2\\2+2x^2\end{pmatrix}.$$

第一行反应为 $2u_1+u_2=4+4x^2$，第二行为 $u_1+5u_2=11+11x^2$。若 $a=1,b=1$，稳态源项是

$$f_1=-2+2x+4+4x^2=2+2x+4x^2,$$
$$f_2=-4+4x+11+11x^2=7+4x+11x^2.$$

一个分量变化会影响另一个方程。如果把 $C$ 的非对角元去掉，同时保持这组 $f$ 和边界不变，解应当改变。这是验证代码确实实现“耦合”的好办法。

## 2. 为什么需要弱形式

### 2.1 从乘一个函数开始

先看一个分量。取一个在边界为零的光滑函数 $v$，方程两边乘 $v$ 后在区域上积分：

$$\int_\Omega d u_t v-\int_\Omega\nabla\cdot(a\nabla u)v
+\int_\Omega(\boldsymbol b\cdot\nabla u)v
+\int_\Omega\sum_j C_{ij}u_jv=\int_\Omega f_i v.$$

对扩散项分部积分：

$$-\int_\Omega\nabla\cdot(a\nabla u)v
=\int_\Omega a\nabla u\cdot\nabla v
-\int_{\partial\Omega}a\nabla u\cdot\boldsymbol n\,v.$$

边界上 $v=0$，最后一项消失。于是只要求 $u$ 和 $v$ 有一阶导数，而不再要求 $u$ 在每处都具有普通意义的二阶导数。这让连续、但单元之间导数可以跳变的分片多项式也能参与计算。

### 2.2 写成向量形式

$H^1(\Omega)$ 可以粗略理解为“函数和一阶弱导数的平方都能积分”的函数集合；$H_0^1$ 再要求边界值为零。令

$$V=[H_0^1(\Omega)]^m.$$

非零边界时，解在具有边界值 $\boldsymbol g$ 的仿射空间中，测试函数仍在 $V$ 中。对任意 $\boldsymbol v\in V$，要求

$$M_t(\partial_t\boldsymbol u,\boldsymbol v)+A_t(\boldsymbol u,\boldsymbol v)=F_t(\boldsymbol v),\tag{3}$$

其中

$$M_t(\boldsymbol w,\boldsymbol v)=\sum_i\int_\Omega d w_i v_i,$$
$$A_t(\boldsymbol w,\boldsymbol v)=\sum_i\left[\int_\Omega a\nabla w_i\cdot\nabla v_i
+\int_\Omega(\boldsymbol b\cdot\nabla w_i)v_i\right]
+\sum_{i,j}\int_\Omega C_{ij}w_jv_i,$$
$$F_t(\boldsymbol v)=\sum_i\int_\Omega f_i v_i.$$

下标 $t$ 提醒我们：系数变了，形式也随时间变化。测试分量 $i$ 决定矩阵行，未知分量 $j$ 决定矩阵列。

### 2.3 为什么不默认使用共轭梯度法

扩散项对 $u,v$ 对称，但对流项一般不对称，反应矩阵也可能不对称。因此最后的矩阵通常不是对称正定矩阵，不能仅因它来自有限元就使用需要对称正定条件的 CG。

本项目采用一般稀疏直接法；以后大规模问题可研究 GMRES 等非对称迭代方法及预条件，但那需要新的性能和稳定性验证。

## 3. 网格与有限元：用有限个数表示函数

### 3.1 一维 P1 的帽子函数

把 $[0,1]$ 划成小区间。每个节点对应一个函数 $\phi_j$：它在本节点等于 1，其他节点等于 0，在相邻单元内是直线。近似解写成

$$u_h(x,t)=\sum_{j=0}^{N}U_j(t)\phi_j(x).$$

因此 $U_j$ 正是节点值。二维三角形 P1 在每个三角形上是一次多项式；三维四面体类似。P2 则使用二次多项式，除顶点外还有边中点自由度。

| 单元 | P1 局部自由度数 | P2 局部自由度数 |
| --- | --- | --- |
| 线段 | 2 | 3 |
| 三角形 | 3 | 6 |
| 四面体 | 4 | 10 |

局部自由度属于一个单元；全局自由度把相邻单元共享的节点合并。连续有限元就是通过共享自由度实现跨单元连续。

### 3.2 单元上的三个小矩阵

在 $[x_L,x_R]$ 上令 $h=x_R-x_L$，两个 P1 基函数是

$$\phi_1=(x_R-x)/h,\quad\phi_2=(x-x_L)/h,\quad
\phi_1'=-1/h,\quad\phi_2'=1/h.$$

常系数时，按“行是测试函数 $i$，列是试探函数 $j$”定义矩阵：

$$K_{ij}=\int a\phi'_j\phi'_i
\quad\Longrightarrow\quad K=\frac a h\begin{pmatrix}1&-1\\-1&1\end{pmatrix},$$
$$B_{ij}=\int b\phi'_j\phi_i
\quad\Longrightarrow\quad B=\frac b2\begin{pmatrix}-1&1\\-1&1\end{pmatrix},$$
$$M_{ij}=\int\phi_j\phi_i
\quad\Longrightarrow\quad M=\frac h6\begin{pmatrix}2&1\\1&2\end{pmatrix}.$$

注意 $B$ 的两行相同，不是对称矩阵。反应 $c$ 为常数时贡献 $cM$，容量 $d$ 贡献 $dM$。源项 $f=1$ 时局部右端为 $(h/2,h/2)^T$。

### 3.3 两个单元的纸笔算例

求 $-u''=2$，$u(0)=u(1)=0$，网格节点为 $0,1/2,1$。每个单元 $h=1/2$，局部扩散矩阵为

$$\begin{pmatrix}2&-2\\-2&2\end{pmatrix}.$$

把共享中点上的贡献相加：

$$A=\begin{pmatrix}2&-2&0\\-2&4&-2\\0&-2&2\end{pmatrix},\qquad
F=\begin{pmatrix}1/2\\1\\1/2\end{pmatrix}.$$

边界 $U_0=U_2=0$，中间一行变为 $4U_1=1$，故 $U_1=1/4$。精确解 $u=x(1-x)$ 在中点也等于 $1/4$，但有限元解是两段直线，并不是整条抛物线。这个例子已经说明“节点正确”与“整个函数精确”之间的区别。

### 3.4 向量块矩阵

若标量空间有 $N_s$ 个自由度，两个分量按

$$U=(U_{1,1},\ldots,U_{1,N_s},U_{2,1},\ldots,U_{2,N_s})^T$$

排列。空间矩阵是

$$A=\begin{pmatrix}K+B+R_{11}&R_{12}\\R_{21}&K+B+R_{22}\end{pmatrix},\qquad
(R_{ij})_{kl}=\int C_{ij}\phi_l\phi_k.$$

反应系数即使随空间变化，也不改变块结构。容量共享且不直接混合分量，所以普通质量矩阵为 $\operatorname{diag}(M_d,M_d)$。

## 4. 数值积分和边界条件

### 4.1 为什么计算机还要做求积

一般系数的积分没有方便的原函数。对一个单元 $K$，使用求积点 $\boldsymbol x_q$ 和归一化权重 $w_q$：

$$\int_K g(\boldsymbol x)\,d\boldsymbol x\approx |K|\sum_qw_qg(\boldsymbol x_q).$$

在本项目使用的 FEALPy 网格规则中，权重与单元度量分开使用。系数、基函数、梯度都在这些点求值，再进行乘加。`q` 是积分规则的阶数参数，不应把它简单看成每个单元恰好有 q 个点。

默认组装使用 `q=p+4`，误差积分再提高两级。平滑算例中这足够支持已报告的收敛观察；厚度远小于网格的指数层可能仍需更细网格或更高阶积分。

### 4.2 非零边界必须影响右端

把自由度分成内部 I 和边界 B：

$$\begin{pmatrix}A_{II}&A_{IB}\\A_{BI}&A_{BB}\end{pmatrix}
\begin{pmatrix}U_I\\g_B\end{pmatrix}=\begin{pmatrix}F_I\\F_B\end{pmatrix}.$$

内部方程是

$$A_{II}U_I=F_I-A_{IB}g_B.\tag{4}$$

所以不能只把边界行替换为单位行而随意丢弃边界列的影响。`DirichletBC` 先扣除 $Ag$，再处理边界行列与右端。本项目每个新时间层都执行这个过程，使用该时刻的 $g$。

时间迭代历史向量包含过去时刻的边界值。这是正常的：历史项和当前边界值在差分式中承担不同角色。

## 5. 从常微分方程组到时间格式

空间离散以后，式 (3) 变为

$$M(t)\dot U(t)+A(t)U(t)=F(t).\tag{5}$$

这里是 $M\dot U$，不是 $\frac d{dt}(MU)$。后者还会出现 $\dot M U$。

### 5.1 后向欧拉

取均匀时间点 $t_n=n\Delta t$。以

$$\dot U(t_{n+1})\approx\frac{U^{n+1}-U^n}{\Delta t}$$

替换导数，所有系数在 $t_{n+1}$ 计算，得

$$\left(M^{n+1}+\Delta t A^{n+1}\right)U^{n+1}
=M^{n+1}U^n+\Delta t F^{n+1}.\tag{6}$$

泰勒展开给出导数近似误差为 $O(\Delta t)$，在足够正则和稳定的条件下全局时间误差通常也是一阶。

### 5.2 BDF2 的系数从哪里来

在 $t_{n+1}$ 展开：

$$U^n=U^{n+1}-\Delta t\dot U+\tfrac12\Delta t^2\ddot U-\tfrac16\Delta t^3U^{(3)}+\cdots,$$
$$U^{n-1}=U^{n+1}-2\Delta t\dot U+2\Delta t^2\ddot U-\tfrac43\Delta t^3U^{(3)}+\cdots.$$

取 $3U^{n+1}-4U^n+U^{n-1}$，常数项和二阶导数项恰好消掉，得到

$$\frac{3U^{n+1}-4U^n+U^{n-1}}{2\Delta t}
=\dot U(t_{n+1})-\frac13\Delta t^2U^{(3)}(t_{n+1})+O(\Delta t^3).$$

代入式 (5)：

$$\left(\tfrac32M^{n+1}+\Delta t A^{n+1}\right)U^{n+1}
=M^{n+1}\left(2U^n-\tfrac12U^{n-1}\right)+\Delta t F^{n+1}.\tag{7}$$

质量矩阵始终取新时刻。把右侧改为 $2M^nU^n-\tfrac12M^{n-1}U^{n-1}$ 会改变对时间项的离散，不再是式 (1) 的直接 BDF2。

### 5.3 第一步怎么办

最开始只有 $U^0$，没有 $U^{-1}$。本项目用后向欧拉计算 $U^1$，之后进入 BDF2。在光滑问题上，一次 BE 起步得到的第一步解误差为 $O(\Delta t^2)$，满足二阶多步方法的起步精度需要。不能把“全程使用 BE 的一阶误差”混同于“只用一次 BE 起步”。

### 5.4 稳定不是所有值都不增长

对测试方程 $y'=-\lambda y$，$\lambda>0$，BE 的放大因子为 $1/(1+\lambda\Delta t)$，幅度不超过 1。BDF2 也是 A 稳定方法；该结论针对线性测试方程的稳定域，不能直接推出任意耦合离散系统都保持正性。

如果反应系数使物理模型本来就增长，正确解也应该增长。例如 $u_t-u=0$ 的解是 $e^t$。数值算法不应把真实增长“稳定”成常数。BDF2 的方法背景可参阅 [SUNDIALS 数学说明](https://sundials.readthedocs.io/en/v6.0.0/cvode/Mathematics_link.html)。

## 6. 强对流为什么会让普通方法振荡

### 6.1 一个可以看见的边界层

考虑

$$-\varepsilon u''+b u'=0,\qquad u(0)=0,u(1)=1,\quad b=1.$$

令 $w=u'$，有 $w'=w/\varepsilon$，于是 $w=K e^{x/\varepsilon}$。再积分并使用边界条件：

$$u(x)=\frac{e^{(x-1)/\varepsilon}-e^{-1/\varepsilon}}{1-e^{-1/\varepsilon}}.\tag{8}$$

当 $\varepsilon$ 很小时，几乎全区间都接近 0，只有出口 $x=1$ 附近快速上升。变化宽度约为 $\varepsilon$。计算式 (8) 时使用指数不大于零的形式，可以避免直接计算 $e^{1/\varepsilon}$ 溢出。

### 6.2 单元 Péclet 数

定义

$$Pe=\frac{|b|h}{2\varepsilon}.$$

它比较“一个单元尺度上的对流作用”与“扩散作用”。一维均匀 P1 Galerkin 的内部行中，左右邻点系数分别是

$$-\varepsilon/h-b/2,\qquad-\varepsilon/h+b/2.$$

当 $Pe>1$，后一个系数变为正数，产生非物理振荡的风险。$\varepsilon=0.001,h=1/32$ 时 $Pe=15.625$，已经很大。整个问题有解并不代表粗网格离散会给出合理的曲线。

![强对流对比及出口放大](../results_final/layer.png)

## 7. SUPG：沿流线加入完整方程残差

### 7.1 加什么，为什么这样加

单元 $K$ 上，第 i 个方程的强残差定义为

$$R_i(\boldsymbol u_h)=d\partial_tu_{h,i}-a\Delta u_{h,i}
-\nabla a\cdot\nabla u_{h,i}+\boldsymbol b\cdot\nabla u_{h,i}
+\sum_j C_{ij}u_{h,j}-f_i.$$

在普通弱形式上增加

$$\sum_K\sum_i\int_K\tau_K(\boldsymbol b\cdot\nabla v_{h,i})R_i(\boldsymbol u_h)=0.\tag{9}$$

测试函数被沿流向调整，所以称为流线上风 Petrov–Galerkin。若代入光滑精确解，每个 $R_i=0$，新增项就消失，这叫一致性。若只添加对流部分、却漏掉非零源项或反应项，精确解一般不再满足修改后的方程。

P1 在单元内部 $\Delta u_h=0$，但变系数时 $\nabla a\cdot\nabla u_h$ 仍然存在。P2 的 Hessian 一般不为零，所以必须计算其迹来得到 Laplacian。这里用的是每个单元内部的强残差，不是把跨单元梯度跳跃当成普通二阶函数值。

### 7.2 稳定参数的两个极限

本项目采用流线参数

$$\tau=\frac{h_b}{2|\boldsymbol b|}\left(\coth Pe-\frac1{Pe}\right),\qquad
Pe=\frac{|\boldsymbol b|h_b}{2a}.$$

仿射单纯形上，令 $\lambda_i$ 为重心坐标，取

$$h_b=\frac{2|\boldsymbol b|}{p\sum_i|\boldsymbol b\cdot\nabla\lambda_i|}.$$

对一维 P1，$h_b$ 就是单元长度 h。P2 使用除以 p 的尺度选择；多维应用是流线尺度的推广，不具有一维节点拟合公式的普遍精确性。变系数时参数在求积点计算。

当 $Pe$ 很小时，直接相减会损失有效数字，使用

$$\coth Pe-1/Pe=Pe/3-Pe^3/45+O(Pe^5),$$

所以 $\tau\approx h_b^2/(12a)$。当 $Pe$ 很大时，$\tau\approx h_b/(2|\boldsymbol b|)$。速度为零时稳定化测试方向也为零，代码直接取 $\tau=0$。

### 7.3 一维为什么能拟合节点

一维常系数、无源无反应、稳态 P1 情形，新增项为

$$\int_K\tau b^2u_h'v_h',$$

相当于把扩散从 $\varepsilon$ 改成

$$\varepsilon_{\rm eff}=\varepsilon+\tau b^2
=\frac{|b|h}{2}\coth Pe.$$

离散齐次方程的相邻差值之比为

$$\frac{\varepsilon_{\rm eff}/h+b/2}{\varepsilon_{\rm eff}/h-b/2}
=e^{bh/\varepsilon}\quad(b>0),$$

与精确指数解在相邻节点间的比例一致。因此该特殊问题中节点可以精确拟合。但基函数仍是分段直线，单元内形状并没有变成指数函数。

### 7.4 SUPG 如何进入代码中的 M、A、F

定义单元测试权 $T_i=\tau\boldsymbol b\cdot\nabla\phi_i$。对标量空间的行 i、列 j，质量修正为

$$S^M_{ij}=\int_K T_i d\phi_j.$$

向量分量 r、s 对应的空间修正为

$$S^{A,rs}_{ij}=\int_K T_i\left[C_{rs}\phi_j+
\delta_{rs}(\boldsymbol b\cdot\nabla\phi_j-a\Delta\phi_j-\nabla a\cdot\nabla\phi_j)\right].$$

右端修正为

$$S^{F,r}_i=\int_K T_i f_r.$$

组装后得到 $\widehat M=M+S^M,\widehat A=A+S^A,\widehat F=F+S^F$，再代入 BE 或 BDF2。因此含时 SUPG 并不是只给空间矩阵加一个黏性项。

方法的原始来源是 [Brooks 与 Hughes 的 SUPG 论文](https://www.sciencedirect.com/science/article/pii/0045782582900718)。

## 8. 怎样证明计算结果值得相信

### 8.1 代数残差只回答一个问题

求得 U 后计算

$$r=\frac{\|AU-F\|_2}{\max(1,\|F\|_2)}.$$

很小的 r 说明线性方程组求解准确。但如果装配时扩散符号写反，错误的方程组也能有很小的残差。因此还需要边界检查、独立源项推导和精确解比较。

### 8.2 制造解：先选答案，再造题目

真实问题通常不知道解。验证代码时却可以先选择一个光滑 $\boldsymbol u$，代入强形式计算 $\boldsymbol f$，同时取相应边界和初值。此时精确答案已知。

本项目令 $w=1+\sum x_k+\prod\sin(\pi x_k)$，取 $\boldsymbol u=(w,2w+1)^T$ 验证空间阶。时间验证换成 $w=1+\sum x_k^2$ 和 $\boldsymbol u=e^t(w,2w+1)^T$。

后一个 w 能被 P2 精确表示，因此缩小时间步时不会很快被空间逼近误差掩盖。如果同时任意改变网格与时间步，就很难判断误差下降到底来自哪一个方向。

另取 $\boldsymbol u=(1+t)(w,2w+1)^T$，空间二次、时间一次，作为接近舍入误差水平的再现测试。手工源项再与独立符号微分交叉核对，避免只依靠同一段实现自我证明。

### 8.3 两种误差范数

向量 $L^2$ 误差为

$$E_0=\left[\int_\Omega\sum_i(u_i-u_{h,i})^2\right]^{1/2}.$$

梯度误差半范数为

$$E_1=\left[\int_\Omega\sum_i\sum_k(\partial_{x_k}u_i-\partial_{x_k}u_{h,i})^2\right]^{1/2}.$$

程序结果中的 `H1` 指 $E_1$。完整 $H^1$ 误差还应合并 $E_0$，成为 $\sqrt{E_0^2+E_1^2}$。

一个数是不是“小”，要结合解的幅度、单位、关注的物理量和网格成本。$10^{-3}$ 对幅度为 1 的浓度与对幅度为 $10^{-6}$ 的浓度有完全不同的意义。图上两条线重合也不等于误差为零。

### 8.4 收敛阶怎么计算

假设 $E(h)\approx Ch^r$，把 h 减半：

$$\frac{E(h)}{E(h/2)}\approx2^r,\qquad
r\approx\log_2\frac{E(h)}{E(h/2)}.$$

例如误差从 0.008 降到 0.002，则阶数约为 2；从 0.008 降到 0.001，阶数约为 3。时间阶同理，用 $\Delta t$ 替代 h。

在解足够光滑、网格规则、问题满足适当稳定性及对偶正则性的条件下，Pp 的常见预期是 $E_1=O(h^p)$、$E_0=O(h^{p+1})$。理论有前提；尖角、薄层、欠积分和粗网格都可能让观测阶数暂时偏离。

![空间误差和阶数](../results_final/spatial.png)

![时间误差和阶数](../results_final/temporal.png)

所有实测数值可在[任务报告](向量含时CDR_任务报告.md)及 CSV 中逐项查到。三维 P1 的最细层仍只有每方向 8 段，所以报告保留其接近但未完全达到理论值的实际阶数。

## 9. 看图与 ParaView 时怎样理解数据

VTU 存放一个时刻的非结构网格和节点场；PVD 把多个 VTU 及其真实时间列成序列。画面中的颜色通常是某个分量或模长，不是误差。只有把数值解与精确解作差，才得到误差图。

本项目导出顶点采样。P2 的边中点自由度参与求解，但普通顶点输出不完整保留二次形函数。因此 ParaView 显示适合观察演化和场分布，精确积分误差仍以求解器的有限元函数计算为准。

![双分量初始、终态与顶点误差](../results_final/field.png)

## 10. 适用范围与后续问题

本讲义对应线性、共享正标量 a 和 d、共享速度 b、反应矩阵 C 耦合、1—3 维盒形区域、全 Dirichlet 边界、P1/P2 和均匀时间步的实现。它不直接覆盖各向异性张量扩散、交叉扩散、非线性反应或其他边界类型。

SUPG 沿流向控制振荡，不提供任意网格上的最大值原理，也不能保证任意反应矩阵的解非负。未解析的横风内层还可能需要自适应网格、横风稳定或通量限制方法。代码不把解截断到某个范围，因此误差与节点上下界可以真实反映算法行为。

## 11. 四个纸笔自检题

1. 一维 $a=1+x,u=x^2,b=2,c=3$，计算稳态 f。答案：$f=-2+3x^2$；中间步骤为 $-(au')'=-2-4x$，再加 $4x+3x^2$。
2. 若两个网格误差为 $1.6\times10^{-3}$ 和 $2.0\times10^{-4}$，h 减半，阶数是多少？答案：3。
3. 二分量、标量空间 5 个自由度，第一分量第 3 个节点在向量中的零起始索引为 2，第二分量同一节点为 7。请写出所有 10 个位置的排列。
4. $a$ 变而 P1 的单元内二阶导数为零，SUPG 的扩散残差是否全为零？答案：不一定，还剩 $-\nabla a\cdot\nabla u_h$。

能独立完成这四题后，可在代码中找到对应位置：`source`、收敛阶计算、`TensorFunctionSpace` 和 `supg_terms`。下一步最有价值的练习是自己选一个二次向量制造解，手算 f，再运行 P2 再现验证。
