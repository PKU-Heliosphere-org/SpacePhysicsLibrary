# 本章: 特征线方法与RH关系模拟一维MHD扰动传播与激波形成演化

## 1.1 核心思想

本章梳理如何把**特征线方法** (method of characteristics, MOC) 与 **Rankine–Hugoniot 跳跃关系** (RH relation) 结合起来，用于模拟一维磁流体力学扰动的传播、非线性陡化、激波形成以及激波后续演化。

核心思想是：

$$
\boxed{\text{光滑区用特征线方法推进；间断面用 RH 条件连接。}} \tag{1.1}
$$

也就是说，MOC 负责描述扰动在光滑流场中的传播；当特征线发生汇聚并相交时，光滑解失效，必须用激波代替多值解。激波两侧状态和激波传播速度由 RH 条件决定。

## 1.2 一维 MHD 控制系统

### 1.2.1 控制方程

考虑一维空间坐标 \(x\)，但保留速度和磁场的三个分量：

\[
\mathbf{v}=(v_x,v_y,v_z),\qquad \mathbf{B}=(B_x,B_y,B_z).
\]

这通常称为 **1D3V MHD** 或一维理想 MHD。理想 MHD 守恒方程可以写成

\[
\frac{\partial U}{\partial t} + \frac{\partial F(U)}{\partial x} = 0. \tag{1.2}
\]

守恒变量为

\[
U = \begin{pmatrix}
\rho \\ \rho v_x \\ \rho v_y \\ \rho v_z \\ B_y \\ B_z \\ E
\end{pmatrix}.
\]

其中 \(B_x\) 由 \(\nabla\cdot\mathbf{B}=0\) 给出：

\[
\frac{\partial B_x}{\partial x}=0,
\]

所以 \(B_x\) 是常数，不作为演化变量。总能量密度为

\[
E = \frac{p}{\gamma-1} + \frac{1}{2}\rho v^2 + \frac{B^2}{2\mu_0}, \tag{1.3}
\]

其中

\[
v^2=v_x^2+v_y^2+v_z^2,\qquad B^2=B_x^2+B_y^2+B_z^2.
\]

为了简洁，也常采用单位制 \(\mu_0=1\)。下面很多公式会默认 \(\mu_0=1\)；若保留 \(\mu_0\)，只需把磁压项写为 \(B^2/(2\mu_0)\)，磁张力项写为 \(B_iB_j/\mu_0\)。

### 1.2.2 波族结构

一维理想 MHD 是双曲系统，具有七个特征波族，其特征速度和波型如表 1.1 所示。

**表 1.1 一维理想 MHD 特征波族**
| 特征速度 | 波型 |
|----------|------|
| \(v_x\pm c_f\) | 快磁声波 |
| \(v_x\pm c_A\) | Alfvén 波 |
| \(v_x\pm c_s\) | 慢磁声波 |
| \(v_x\) | 熵波/接触间断 |

法向 Alfvén 速度为

\[
c_A = \frac{|B_x|}{\sqrt{\mu_0\rho}}. \tag{1.4}
\]

总 Alfvén 速度为

\[
v_A^2 = \frac{B^2}{\mu_0\rho}. \tag{1.5}
\]

声速为

\[
a^2 = \frac{\gamma p}{\rho}. \tag{1.6}
\]

快、慢磁声速度满足

\[
c_{f,s}^2 = \frac{1}{2}\Bigl[ a^2+v_A^2 \pm \sqrt{(a^2+v_A^2)^2 - 4a^2 c_A^2} \Bigr], \tag{1.7}
\]

其中 \(+\) 号对应快磁声速度 \(c_f\)，\(-\) 号对应慢磁声速度 \(c_s\)。

## 1.3 特征线方法基础

### 1.3.1 基本思想

对于一般一维双曲系统

\[
\frac{\partial U}{\partial t} + A(U)\frac{\partial U}{\partial x} = 0,
\qquad A(U)=\frac{\partial F}{\partial U},
\]

设 \(A\) 有特征值和左右特征矢量：

\[
A r_k = \lambda_k r_k,\qquad l_k A = \lambda_k l_k.
\]

左乘 \(l_k\) 得

\[
l_k\cdot\frac{\partial U}{\partial t} + \lambda_k\, l_k\cdot\frac{\partial U}{\partial x} = 0.
\]

沿第 \(k\) 族特征线

\[
\frac{dx}{dt} = \lambda_k(U) \tag{1.8}
\]

有

\[
l_k\cdot \frac{dU}{dt} = 0. \tag{1.9}
\]

这就是第 \(k\) 族特征线上的**兼容关系**。因此，MOC 的基本结构是

\[
\boxed{\frac{dx}{dt}=\lambda_k(U),\qquad l_k\cdot dU=0.} \tag{1.10}
\]

如果 \(l_k\cdot dU\) 可以积分为某个标量函数的全微分 \(dR_k = l_k\cdot dU\)，则得到黎曼不变量：

\[
R_k = \text{constant along } \lambda_k \text{ characteristic}.
\]

### 1.3.2 黎曼不变量与左特征矢量的关系

需要强调一点：**\(l_k A U\) 通常不是黎曼不变量**。正确的关系是

\[
l_k A(U)\frac{\partial U}{\partial x}
= \lambda_k l_k\frac{\partial U}{\partial x}.
\]

与时间项合并后得到

\[
l_k\cdot\left( \frac{\partial U}{\partial t} + \lambda_k\frac{\partial U}{\partial x}\right)=0.
\]

沿特征线 \(\frac{dx}{dt}=\lambda_k\) 有 \(l_k\cdot dU=0\)。若存在积分因子 \(\mu(U)\) 使得

\[
dR_k = \mu(U) l_k\cdot dU,
\]

则 \(R_k\) 为黎曼不变量。因此：

\[
\boxed{\text{左特征矢量给出特征兼容微分关系；黎曼不变量是该微分关系可积时得到的标量。}} \tag{1.11}
\]

对于一维等熵气体动力学，可以得到简单黎曼不变量 \(J_\pm = u\pm \frac{2c}{\gamma-1}\)。但对于完整一维 MHD，通常不存在所有波族都简单可写的全局黎曼不变量，因此实际 MHD 特征线方法中更多使用微分形式 \(l_k\cdot dU=0\)。

### 1.3.3 光滑区的 MOC 推进

在没有激波、间断、强耗散的光滑区域，一维 MHD 可以用特征分解推进。

设当前时间层为 \(t^n\)，空间网格为 \(x_i\)，状态为 \(U_i^n\)。对每个点计算

\[
\lambda_{k,i}^n,\quad l_{k,i}^n,\quad r_{k,i}^n,\qquad k=1,\dots,7.
\]

第 \(k\) 族特征线从 \(t^n\) 到 \(t^{n+1}\) 的近似足点为

\[
x_{k,i}^{\text{foot}} = x_i - \lambda_{k,i}^n\Delta t. \tag{1.12}
\]

在旧时间层从同一光滑区域插值得到 \(U^{\text{foot}}_{k,i}\)。沿该特征线的兼容关系近似为

\[
l_{k,i}\cdot\left( U_i^{n+1} - U^{\text{foot}}_{k,i} \right) = 0. \tag{1.13}
\]

对七个波族写出七个关系：

\[
\begin{aligned}
l_{1,i}\cdot\left( U_i^{n+1} - U^{\text{foot}}_{1,i} \right) &= 0,\\
l_{2,i}\cdot\left( U_i^{n+1} - U^{\text{foot}}_{2,i} \right) &= 0,\\
&\cdots\\
l_{7,i}\cdot\left( U_i^{n+1} - U^{\text{foot}}_{7,i} \right) &= 0.
\end{aligned}
\]

这相当于建立线性系统

\[
L_i\, U_i^{n+1} = b_i, \tag{1.14}
\]

其中 \(L_i\) 的每一行是一个左特征矢量 \(l_{k,i}\)，右端项为 \(b_{k,i}=l_{k,i}\cdot U^{\text{foot}}_{k,i}\)。求解即可得到 \(U_i^{n+1}\)。这就是 MHD 光滑区 MOC 更新的一个基本形式。

### 1.3.4 MOC 的局限

MOC 依赖于一个关键假设：\(U(x,t)\) 是光滑函数。

一旦特征线相交，MOC 将产生多值解。例如某一波族的特征速度为 \(\lambda_k(U)\)，如果后方特征线比前方传播得更快，则特征线会汇聚：

\[
\frac{\partial \lambda_k}{\partial x} < 0.
\]

随着时间演化，特征线可能相交：\(x_{k,i}^{n+1} > x_{k,i+1}^{n+1}\)。特征线相交意味着 \(\partial U/\partial x \to \infty\)，这是波形陡化、激波形成的数学信号。因此：

\[
\boxed{\text{MOC 不能跨越激波继续使用；激波必须作为移动间断面单独处理。}} \tag{1.15}
\]

## 1.4 激波形成检测

### 1.4.1 特征线交叉判据

对第 \(k\) 族特征线，从 \(t^n\) 推进到 \(t^{n+1}\)：

\[
x_{k,i}^{n+1} = x_i + \lambda_{k,i}^n\Delta t. \tag{1.16}
\]

若相邻特征线发生反序 \(x_{k,i}^{n+1} > x_{k,i+1}^{n+1}\)，则说明第 \(k\) 族特征线相交，通常表示第 \(k\) 族压缩波陡化为激波。

### 1.4.2 特征速度梯度判据

也可以直接检查

\[
\frac{\partial \lambda_k}{\partial x} < 0.
\]

离散形式为

\[
D_{k,i} = \frac{\lambda_{k,i+1} - \lambda_{k,i}}{x_{i+1}-x_i}. \tag{1.17}
\]

若 \(D_{k,i} < -\epsilon_\lambda\)（\(\epsilon_\lambda\) 为小阈值），则该区域为特征汇聚区。

### 1.4.3 压缩判据

MHD 激波通常伴随法向压缩。一维中法向速度为 \(v_x\)，压缩区满足

\[
\frac{\partial v_x}{\partial x} < 0.
\]

对于快磁声压缩，还常伴随 \(\partial\rho/\partial x\) 和 \(\partial p/\partial x\) 的陡化。但仅有压缩并不一定意味着激波，必须结合特征线汇聚和熵条件判断。

### 1.4.4 熵增判据

物理激波必须满足熵增条件。对于理想气体

\[
s \propto \ln\left(\frac{p}{\rho^\gamma}\right).
\]

跨激波要求

\[
\Delta s = \ln\left[ \frac{p_2/\rho_2^\gamma}{p_1/\rho_1^\gamma} \right] > 0. \tag{1.18}
\]

### 1.4.5 激波形成位置的确定

若第 \(k\) 族特征线 \(i\) 与 \(i+1\) 在时间步内相交，其轨迹近似为

\[
x_i(t) = x_i^n + \lambda_{k,i}^n (t-t^n),\qquad
x_{i+1}(t) = x_{i+1}^n + \lambda_{k,i+1}^n (t-t^n).
\]

令二者相等，解得相交时间

\[
\tau = \frac{x_{i+1}^n - x_i^n}{\lambda_{k,i}^n - \lambda_{k,i+1}^n}. \tag{1.19}
\]

如果 \(0<\tau<\Delta t\)，则相交发生在当前时间步内。激波初始形成位置可取为

\[
x_s^{\text{form}} = x_i^n + \lambda_{k,i}^n \tau. \tag{1.20}
\]

## 1.5 Rankine–Hugoniot 跳跃关系

### 1.5.1 RH 条件的一般形式

对于守恒律

\[
\frac{\partial U}{\partial t} + \frac{\partial F(U)}{\partial x} = 0,
\]

若存在一个以速度 \(S\) 运动的间断面，则跨间断的 RH 条件为

\[
S[U] = [F], \tag{1.21}
\]

其中 \([U]=U_2-U_1\)，\([F]=F(U_2)-F(U_1)\)。这里 \(U_1\)、\(U_2\) 分别为上游和下游状态，\(S\) 为激波速度。对于一维 MHD，RH 条件仍然可写为此形式，但展开后具有明确的物理含义。

### 1.5.2 一维 MHD 的 RH 条件

在激波参考系中定义法向相对速度

\[
w = v_x - S. \tag{1.22}
\]

一维 MHD 的 RH 条件包括以下几类。

#### 质量通量守恒

\[
[\rho w] = 0 \quad\Longrightarrow\quad
\rho_1(v_{x1}-S) = \rho_2(v_{x2}-S). \tag{1.23}
\]

定义质量通量 \(m=\rho w\)，则 \(m_1=m_2=m\)。

#### 法向磁场连续

由 \(\nabla\cdot\mathbf{B}=0\) 得

\[
[B_x] = 0 \quad\Longrightarrow\quad B_{x1}=B_{x2}=B_x. \tag{1.24}
\]

#### 法向动量守恒

\[
\left[ \rho w v_x + p + \frac{B_t^2}{2\mu_0} \right] = 0,
\tag{1.25}
\]

或用相对速度写为

\[
\left[ \rho w^2 + p + \frac{B_t^2}{2\mu_0} \right] = 0.
\]

其中 \(B_t^2 = B_y^2+B_z^2\)。

#### 切向动量守恒

\[
\left[ \rho w \mathbf{v}_t - \frac{B_x\mathbf{B}_t}{\mu_0} \right] = 0, \tag{1.26}
\]

其中 \(\mathbf{v}_t=(v_y,v_z)\)，\(\mathbf{B}_t=(B_y,B_z)\)。展开为

\[
\begin{aligned}
\left[ \rho w v_y - \frac{B_x B_y}{\mu_0} \right] &= 0,\\
\left[ \rho w v_z - \frac{B_x B_z}{\mu_0} \right] &= 0.
\end{aligned}
\]

#### 切向电场连续

理想 MHD 中 \(\mathbf{E}=-\mathbf{v}\times\mathbf{B}\)，跨激波要求切向电场连续。对一维激波有

\[
[ w\mathbf{B}_t - B_x\mathbf{v}_t ] = 0, \tag{1.27}
\]

即

\[
\begin{aligned}
[ w B_y - B_x v_y ] &= 0,\\
[ w B_z - B_x v_z ] &= 0.
\end{aligned}
\]

这两个条件也可看作切向磁场感应方程的 RH 条件。

#### 能量通量守恒

总能量密度 \(E\) 由 (1.3) 给出。总能量通量在 \(x\) 方向为

\[
F_E = \left( E + p + \frac{B^2}{2\mu_0} \right)v_x - \frac{(\mathbf{v}\cdot\mathbf{B})B_x}{\mu_0}. \tag{1.28}
\]

运动激波条件下

\[
S[E] = [F_E], \tag{1.29}
\]

或者在激波参考系中用 \(w\) 写成相应的能量通量连续形式。

### 1.5.3 MHD 激波分类

RH 条件可以给出不同类型的 MHD 间断，主要包括：

1. 快激波；
2. 慢激波；
3. 中间激波；
4. 旋转间断；
5. 接触间断；
6. 切向间断。

在模拟扰动陡化时，最常关注快磁声激波和慢磁声激波。

#### 快激波

快激波通常满足

\[
|w_1| > c_{f1},\qquad |w_2| < c_{f2}.
\]

即上游相对于激波为超快磁声，下游为亚快磁声。典型特征：

\[
\rho_2>\rho_1,\quad p_2>p_1,\quad |\mathbf{B}_{t2}| > |\mathbf{B}_{t1}|.
\]

#### 慢激波

慢激波跨越慢磁声特征速度，典型特征是切向磁场减弱：

\[
|\mathbf{B}_{t2}| < |\mathbf{B}_{t1}|.
\]

在磁重联外流区中，慢激波常作为 Petschek 型重联模型中的重要结构。

#### 旋转间断

旋转间断通常满足 \(\rho_2\approx\rho_1\)，\(p_2\approx p_1\)，\(|\mathbf{B}_2|\approx|\mathbf{B}_1|\)，但切向磁场方向发生旋转。它主要与 Alfvén 波族有关，不是压缩激波。

#### 接触间断

接触间断满足

\[
v_{x1}=v_{x2}=S,\qquad p_{\text{tot},1}=p_{\text{tot},2},
\]

但密度可以跳跃 \(\rho_1\neq\rho_2\)。其中总压为 \(p_{\text{tot}} = p + \frac{B^2}{2\mu_0}\)。

## 1.6 MOC 与 RH 耦合求解

### 1.6.1 耦合思想

MOC 与 RH 结合时，计算域被划分为：

```
光滑区  |  激波  |  光滑区
MOC     |  RH    |  MOC
```

在光滑区 \(U_t + A(U)U_x = 0\) 用特征线方法推进；在激波处 \(S[U]=[F]\) 用 RH 条件连接两侧状态。因此激波相当于一个**移动内边界**，其位置为 \(x_s(t)\)，速度为 \(S=dx_s/dt\)。激波两侧状态分别为 \(U_1=U(x_s^+,t)\)（上游）、\(U_2=U(x_s^-,t)\)（下游）。例如右行激波：

```
下游区                 激波                 上游区
U2                     xs                   U1
```

### 1.6.2 下游 MOC 兼容关系的必要性

以一维气体动力学为例，未知量可以是 \(\rho_2, u_2, p_2, S\)。RH 条件给出三个方程（质量、动量、能量） ，还缺一个方程。这个方程来自下游光滑区中进入激波的特征线。对于右行声学激波，进入激波的下游特征通常是 \(C_-\)，其兼容关系为

\[
u_2 - \frac{2c_2}{\gamma-1} = J_-^d,
\]

其中 \(J_-^d\) 是从下游光滑区沿 \(C_-\) 追踪到激波处得到的值。因此，激波处下游状态不仅要满足 RH，还要与下游光滑流场相容。

对于 MHD，类似关系写为

\[
l_k\cdot(U_2 - U^{\text{foot}}) = 0,
\]

这里 \(l_k\) 是从下游进入激波的特征族对应的左特征矢量，\(U^{\text{foot}}\) 是该特征线在旧时间层的足点状态。因此：

\[
\boxed{\text{RH 负责跨激波守恒；MOC 兼容关系负责把激波后的状态接入下游光滑区。}} \tag{1.30}
\]

### 1.6.3 激波拟合的未知量与方程组

对于一维 MHD，激波处的未知量包括下游状态

\[
U_2 = (\rho_2, v_{x2}, v_{y2}, v_{z2}, B_{y2}, B_{z2}, p_2),
\]

以及激波速度 \(S\)，共 8 个未知量。若上游状态 \(U_1\) 已知，则 RH 条件提供 7 个守恒关系，需要从下游进入激波的特征关系来保证解的唯一性和相容性，从而闭合方程组。

实际 shock-fitting 中，方程组可组织为

\[
\mathcal{R}(U_2, S) = 0,
\]

其中残差包括：RH 质量守恒、法向动量守恒、两个切向动量守恒、两个切向电场连续条件、RH 能量守恒，以及下游入射特征兼容关系，同时波型约束与熵条件作为筛选条件。

### 1.6.4 上游状态的确定

以右行激波为例，下游在左、上游在右。从激波右侧光滑区插值得到 \(U_1 = U(x_s^+,t)\)。**注意：不能跨过激波插值**，求 \(U_1\) 时只能使用激波右侧的网格点。同理，左行激波的上游在左侧，只能使用激波左侧的点插值。

### 1.6.5 下游 MOC 兼容关系的确定

以右行快激波为例，下游在左侧。从激波处反向沿下游进入激波的特征线追踪到旧时间层，得到足点

\[
x_k^{\text{foot}} \approx x_s^n - \lambda_k(U_L^n)\Delta t,
\]

从下游光滑区插值得到 \(U_k^{\text{foot}}\)。然后写出兼容关系

\[
l_k(U_L^n)\cdot\left( U_2 - U_k^{\text{foot}} \right) = 0. \tag{1.31}
\]

对于不同激波类型，进入激波的特征族不同，须根据激波参考系中各特征速度与 \(S\) 的相对方向来判断。

### 1.6.6 联立 RH+MOC 求解

构造未知量

\[
X = (\rho_2, v_{x2}, v_{y2}, v_{z2}, B_{y2}, B_{z2}, p_2, S)^T,
\]

已知上游状态 \(U_1\) 和下游 MOC 兼容关系 (1.31)。构造残差 \(R(X)=0\)，包含以下部分。

#### 质量 RH 残差
\[
R_\rho = \rho_1(v_{x1}-S) - \rho_2(v_{x2}-S). \tag{1.32}
\]

#### 法向动量 RH 残差
\[
R_{mx} = \Bigl[ \rho w^2 + p + \frac{B_t^2}{2\mu_0} \Bigr]_1
      - \Bigl[ \rho w^2 + p + \frac{B_t^2}{2\mu_0} \Bigr]_2. \tag{1.33}
\]

#### 切向动量 RH 残差
\[
\begin{aligned}
R_{my} &= \Bigl[ \rho w v_y - \frac{B_x B_y}{\mu_0} \Bigr]_1
       - \Bigl[ \rho w v_y - \frac{B_x B_y}{\mu_0} \Bigr]_2,\\
R_{mz} &= \Bigl[ \rho w v_z - \frac{B_x B_z}{\mu_0} \Bigr]_1
       - \Bigl[ \rho w v_z - \frac{B_x B_z}{\mu_0} \Bigr]_2. \tag{1.34}
\end{aligned}
\]

#### 切向电场连续残差
\[
\begin{aligned}
R_{By} &= (w_1 B_{y1} - B_x v_{y1}) - (w_2 B_{y2} - B_x v_{y2}),\\
R_{Bz} &= (w_1 B_{z1} - B_x v_{z1}) - (w_2 B_{z2} - B_x v_{z2}). \tag{1.35}
\end{aligned}
\]

#### 能量 RH 残差
\[
R_E = (F_{E1} - S E_1) - (F_{E2} - S E_2), \tag{1.36}
\]

其中 \(F_E\) 由 (1.28) 给出。

#### 下游 MOC 兼容残差
\[
R_{\text{MOC}} = l_k\cdot( U_2 - U_k^{\text{foot}} ). \tag{1.37}
\]

联立求解 \(R(X)=0\)，可用 Newton 法。

### 1.6.7 Newton 求解

Newton 迭代形式为

\[
J(X^m)\,\Delta X^m = -R(X^m),\qquad X^{m+1}=X^m+\Delta X^m,
\]

其中 \(J_{ij}=\partial R_i/\partial X_j\) 可解析计算或有限差分近似。初值可取 \(U_2^{(0)}=U_L^n\)（下游最近网格点状态），\(S^{(0)} = \lambda_k(U_1)\) 或上一时间步激波速度 \(S^n\)。收敛标准可设为 \(\|\Delta X\|/\|X\|<10^{-8}\)。

### 1.6.8 激波物理性检查

求得 \(U_2\) 和 \(S\) 后，必须进行物理检验。

#### 正定性
\[
\rho_2 > 0,\qquad p_2 > 0. \tag{1.38}
\]

#### 熵增
\[
\Delta s = \ln\left[ \frac{p_2/\rho_2^\gamma}{p_1/\rho_1^\gamma} \right] > 0. \tag{1.39}
\]

#### 快激波条件
上游超快、下游亚快：\(|w_1|>c_{f1},\;|w_2|<c_{f2}\)，且 \(\rho_2>\rho_1,\;p_2>p_1,\;|\mathbf{B}_{t2}|>|\mathbf{B}_{t1}|\)。

#### 慢激波条件
跨越慢磁声速度，且通常 \(|\mathbf{B}_{t2}|<|\mathbf{B}_{t1}|\)。

#### Lax 条件
物理激波应满足信息从两侧进入激波。对于右行第 \(k\) 族激波，特征速度满足

\[
\lambda_k(U_1) > S > \lambda_k(U_2);
\]

左行激波则相反。其物理意义是：

\[
\boxed{\text{激波是信息汇聚面，而不是信息发散面。}} \tag{1.40}
\]

### 1.6.9 激波位置更新

若激波速度为 \(S\)，可采用预测‑校正法提高时间精度：

\[
\begin{aligned}
x_s^* &= x_s^n + S^n\Delta t,\\
x_s^{n+1} &= x_s^n + \frac{1}{2}(S^n + S^*)\Delta t. \tag{1.41}
\end{aligned}
\]

### 1.6.10 激波作为移动内边界

激波更新后，计算域需重新划分。若 \(x_j < x_s^{n+1} < x_{j+1}\)，则 \(i\le j\) 为下游区，\(i\ge j+1\) 为上游区。将 \(U(x_s^-)=U_2\)、\(U(x_s^+)=U_1\) 作为两侧光滑区的边界条件，后续 MOC 推进分别在下游区和上游区内进行，不得跨越激波追踪特征线。

## 1.7 完整算法与程序结构

### 1.7.1 完整算法流程

```
初始化：
    给定 U_i^0, B_x
    shock_list = []

for n = 0,1,2,...:
    1. 对每个光滑区计算特征速度 λ_k,i 及左特征矢量 l_k,i
    2. 用 MOC 预测光滑区解 U_i^{n+1,*}
    3. 检查特征线交叉，生成激波候选
    4. 若有新激波，创建 shock object（位置、类型、初速）
    5. 对每个激波：
        a. 确定上、下游侧
        b. U1 ← 上游侧插值
        c. 追踪下游进入激波的特征线，获得 U_foot
        d. 构造 RH 残差与 MOC 兼容残差，Newton 求解 U2,S
        e. 检查正定、熵增、波型与 Lax 条件
        f. 更新激波位置 x_s
        g. 施加边界状态 U2 (下游), U1 (上游)
    6. 根据新激波位置重新划分光滑区
    7. 接受 MOC 更新值，并与激波边界一致
```

### 1.7.2 推荐程序模块

若要编写教学或研究原型代码，建议按以下模块组织：

```
module mhd_state             # 保存 rho, vx, vy, vz, By, Bz, p, Bx
module mhd_characteristics  # 计算快/慢磁声速度、特征值、左特征矢量
module moc_solver           # 追踪足点、插值、求解特征兼容系统
module shock_detector       # 检测特征线交叉、估计激波形成位置
module rh_solver            # 计算 RH 残差、Newton 求解、熵检验、分类
module shock_manager        # 创建激波、移动激波、重建区域、施加边界
main_loop                   # MOC 光滑更新 → 检测激波 → 求解 RH → 移动激波 → 重建光滑区
```

### 1.7.3 最小可行模型建议

为避免一开始就陷入完整 MHD 特征矢量的复杂性，建议分层实现。

#### 第一阶段：一维等熵声学 MOC + RH
变量 \(\rho, u, p\)，使用 \(J_\pm = u \pm \frac{2c}{\gamma-1}\)。验证特征传播、交叉、RH 激波拟合及熵增。

#### 第二阶段：一维 Euler MOC + RH
处理声波、接触间断、激波与稀疏波，允许熵变化。

#### 第三阶段：简化一维 MHD
取 \(v_y=v_z=0,\; B_z=0\)，只保留 \(\rho, v_x, p, B_y\)，研究快/慢磁声扰动的陡化。

#### 第四阶段：完整 1D3V MHD
实现全变量 \(\rho, v_x, v_y, v_z, B_y, B_z, p\)，处理快磁声波、Alfvén 波、慢磁声波、接触间断、快/慢激波及旋转间断。

## 1.8 方法讨论与总结

### 1.8.1 与有限体积 Riemann 方法的关系

MOC–RH 方法属于 **shock fitting**，显式追踪零厚度激波面。现代 MHD 数值模拟中更常用有限体积 Godunov 型方法进行 **shock capturing**。表 1.2 比较了两类方法。

**表 1.2 MOC–RH 方法与 Godunov 型方法对比**
| 方法 | 激波处理方式 | 优点 | 缺点 |
|------|--------------|------|------|
| MOC–RH shock fitting | 显式追踪激波面 | 激波位置精确、物理结构清楚 | 多激波、多维、拓扑变化困难 |
| Godunov shock capturing | 用数值通量捕捉激波 | 鲁棒、适合复杂 MHD | 激波有数值厚度 |
| 纯 MOC | 只适合光滑区 | 物理直观、低耗散 | 无法处理激波 |

因此，对于教学、理论分析、单激波演化，MOC–RH 很有价值；对于复杂太阳风、日冕和湍流 MHD 模拟，通常使用有限体积 shock-capturing 方法更稳健。

### 1.8.2 本章总结

特征线方法和 RH 条件的关系可以用一句话概括：

\[
\boxed{\text{MOC 描述光滑区内信息沿特征线传播；RH 描述非光滑间断面两侧守恒跳跃。}} \tag{1.42}
\]

在一维 MHD 扰动传播问题中，基本逻辑为：

1. 把 MHD 方程写成准线性双曲系统；
2. 求特征速度和左特征矢量；
3. 在光滑区沿特征线推进；
4. 检查特征线是否汇聚和交叉；
5. 一旦交叉，认为光滑解失效；
6. 插入激波作为移动内边界；
7. 从上游光滑区获得 \(U_1\)；
8. 从下游光滑区获得 MOC 兼容关系；
9. 联立 RH 条件与兼容关系求 \(U_2\) 和 \(S\)；
10. 检查正定性、熵增、Lax 条件和波型条件；
11. 更新激波位置；
12. 重新划分光滑区并继续推进。

最终形成的计算框架是：

\[
\boxed{\text{光滑区：MOC}\;+\;\text{激波面：RH}\;+\;\text{物理筛选：熵条件与波型条件}.} \tag{1.43}
\]

这套方法特别适合用于理解一维 MHD 扰动如何从线性波传播进入非线性陡化，并最终形成快激波、慢激波或其他 MHD 间断结构。