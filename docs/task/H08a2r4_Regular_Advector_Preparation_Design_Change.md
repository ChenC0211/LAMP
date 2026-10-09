# H08a2r4 / 生产几何准备的实施前设计

2026-10-05。输入为 H08a2r3 冻结的 60 个实际布局、144 组边界检查、12 个收敛布局与额外 sine 检查。本项只部署准备层；生产 B、atom、characteristic、F/R/uR 和 KEP 不切换。四项材料硬约束不变，不增加 pointwise constant divergence、commuting 或 R=0。

现有 `CchAdvectorReconstruction` 实际持有不可变的 `Detail::CchRegularAdvectorGeometry`。它持有原 momentum geometry 身份、实际 active full sides 与 AMR master DAG、P2 最小范数权重、Q4 装配关系、实际 face 均值权重和修正支持边界。没有 mesh 裸引用、速度、密度、界面或时钟状态。应用只产生待消费的代数数据（每 cell 两个 25 元 Q4 数组和每实际 face 的修正幅度），不发布第二个 velocity callback。

边界仍由现有显式 `CchAdvectorBoundary` 绑定。几何准备为每个独立 datum/component 预制无 wall 和 zero-normal 两种模式；只有可能受实际边界 face 影响的 datum 才有第二模式。触发集合是与 datum 所在 master full side 相交的实际外边界 face。求值只检查声明并选权重，禁止 SVD 或重新拟合。RZ 轴奇偶条件始终施加；重复约束复用同一行。非零 wall/axis 输入拒绝，不能修正调用方数据。

P2 使用物理平移/尺度坐标 `[1,X,Y,X²,XY,Y²]` 与真实 XY/RZ face mean。从 incident cells 的一圈面邻居开始，统一扩展至满允许列秩、四圈或分区耗尽。SVD 阈值为 `128*epsilon*max(nFaces,6)*sigmaMax`；所有行采用同一最小范数定义，保留 rank、deficit、实际支持、权重范数及缺秩 null polynomial。每个真实 AMR 子面独立校正均值，支持宽度为相邻 cell 法向宽度的一半。Q4 基函数与 H08a2r3/Core 相同，先装配共享数据再计算 face mean。不把旧三系数 atom 截断近似为新来源。

固定文件预算：生产/build 5 个：`CchSpatialAdvector.h/.cpp`、`Detail/CchRegularAdvectorGeometry.h`、`CchRegularAdvectorGeometry.cpp`、根 CMake；测试 1 个：原 momentum geometry contract。短测试覆盖身份/寿命/拒绝条件和准备应用；完整矩阵通过该测试的临时审计入口运行，代码和数据归档到本任务 checks 后移除入口。

门槛沿用 TODO：原确定性 B/atom/F/R 输出不变；同原型、Core 的 C1/face mean/Gauss/可辨识 P2 缺陷 <=1e-11；原 source 值/导数阶 >=2.5/1.8 或舍入精确；保留 sine 与全部 rank 证书。报告 C++ n8/16/32 准备/应用成本、权重规模和应用无 SVD。通过也只表示准备层通过；部署实际 B/J/atom 后仍须重跑 H10o 与 H10b-v2，再进入 H08a3/H08a，H11 未启动。
