# 空间物理知识库

欢迎来到 `SpacePhysicsLibrary`。

这个知识库基于 MkDocs Material 搭建，侧边栏层级完全由根目录下的 `mkdocs.yml` 控制。即使移动端或窄屏下侧边栏会折叠，仍然可以从本页直接进入各个主题。

## 内容导航

### 基础理论

- [磁流体力学](theory/mhd.md)
- [线性等离子体理论](theory/linear-plasma.md)
- [准线性理论](theory/quasi-linear.md)
- [等离子体不稳定性](theory/instability.md)

### 探测技术

- [原位磁场探测](observation/mag-detect.md)
- [原位粒子探测](observation/particle-detect.md)

### 数值模拟方法

- [MHD数值求解](numeric/mhd-solver.md)
- [粒子模拟PIC](numeric/pic.md)

## 维护说明

新增知识点时，只需：

1. 在 `docs/` 对应分类目录中新建 Markdown 文件。
2. 在 `mkdocs.yml` 的 `nav` 中追加对应条目。
3. 刷新或重新启动 MkDocs 服务，侧边栏将自动更新。

## 数学公式测试

行内公式：$\\omega = k v_A$

块级公式：

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{v}) = 0
$$
