# H10j0 / 实施前运动学设计

2026-10-04。只实施 [任务书](H10j0_Spatial_Advector_Reconstruction_TODO.md) 的 actual compatible advector。H10i2 188 payload 与 canonical 已复核；前一轮为 progress。生产固定新 header/cpp 和现 CMake 单一源登记，测试只扩展现 momentum geometry contract。

prepare 持有原 momentum geometry，以实际选中压力 cell 邻接和同分量 face 形成 closed one-ring；使用 H10i 的物理面均值、归一化三列 JacobiSVD，threshold=1e-12，最小范数切向 functional。保留 rank/deficit/weights，轴 radial 不参与除 r 的拟合。evaluate 只乘 prepared 权重，复制 actual native values；保留 caller 声明 geometry/full clock，生成私有 fresh identity。声明不冒称 H14 commit。

体积坐标 X=x 或 s=r²/2，Y=y/z。side traces 保持原 A_f*u_f，B=(u_x,u_y) 或 (r*u_r,u_z)。实际 cell side piece 决定 L/R/Bottom/Top，A=(R-L)/dx、C=(Top-Bottom)/dy。B0 为 normal transfinite，减去 H_X=integral(C-meanC)dX、H_Y=integral(A-meanA)dY。每条 integral 分实际对边 knots 积分；边界 lift 精确置零，局部 divergence=meanA+meanC，对照原 D_I/V。RZ axis 正规极限取零，不除轴半径。

prepare 按 (normalAxis,plane) 对原 primal/dual segments 与 kinematic knots overlay。一份 canonical atom 必须只有同一 primalVelocity 或 ownerPressure 来源；若矛盾拒绝。primal/dual 原顺序各持有 atom 索引；negative/positive incidence 提供相反号，无材料配分。一次 line sorting O(N log N)，存储及 evaluate 工作与真实局部数量线性。

atom normal polynomial 用距离实际 begin/end 的局部 chart 表达（XY 一次，RZ axial 二次 r）；shifted integral 直接使用 fromEnd/offset/length，不重建两全局 root。cell B/原 native 同点两分量只读查询；transverse 用实际 half-dual overlap 体积权重，normal 来自唯一 primal trace。

验收先 independent affine/physical mean、D_I 与局部 divergence、所有 atom 和真实 views 加回，后非仿射独立局部积分/系数参考。16 full +48 selected 布局保留，包括原 n4 selected2/3 信息缺口。两次 evaluate 的 reconstruction storage、SVD 与 atoms 不变，identity/系数改变；来源、clock、FP/axis/nonfinite/重网格拒绝和生命周期检查。原研究数值路径保持，新增只读 actual owner 模式，不接入 phase/material/native convection。T0 必要短契约和 architecture/layout 通过后冻结，只关闭 H10j0，下一 H10j。
