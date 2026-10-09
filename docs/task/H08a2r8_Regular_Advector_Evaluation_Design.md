# H08a2r8 / 共同解析内核实施前设计

2026-10-05。上一轮为 progress；r7 冻结 119 payload、17 当前文件及 1382 源码/输入已核验。检查目录 `runs/checks/H08a2r8/20261005_151948_867`。

范围 production 2：原 `Detail/CchRegularAdvectorGeometry.h` 与 `CchRegularAdvectorGeometry.cpp`；test 1：原 `cch_momentum_geometry_contract.cpp`。不新增源码、目标、App、benchmark 模式或生产开关。

准备对象新增 `cellJet(coefficients,pressure,point)` 和 `normalPolynomial(coefficients,pressure,normalAxis,begin,end)`。`apply` 的原代数系数附带私有原几何 shared identity，查询检查身份、库存形状、索引、范围和实际使用的有限系数。该数据仍是 Detail 暂存系数，不携带时钟或可发布速度身份；生产实际 owner 后续必须统一绑定。无网格借用，查询不动态拟合或进行全网格验证。

值/J 使用原 Q4 cardinal tensor 与 δ b(xi)g(s)n 的解析导数，b=(1−xi²)²、g=16s²(1−s)²。法向支持宽度取当前 pressure 在该 face 的真实单侧半宽；只枚举原 cell.faces。支撑外贡献零，端点值/一阶导数零；边界面中心两侧取值和一阶导数与宽度无关。轴使用物理 u，不引入除 r；Gauss 检查用 2π[r div(u)+u_r]。

cellJet 返回由原 supportKnots 定义的当前分片闭包；内部 knot 选正侧，cell 上端用最后分片。normalPolynomial 必须在同一 cell 中沿 canonical 正切向，不能跨过任何严格内部支撑 knot。它用 cardinal 幂系数的解析仿射复合提取两端全部四次项，并使用同一面修正形状；不从点样拟合多项式。独立积分继续使用 r7 的 `(offset,width)` 内核。

独立测试把 Q4 原 DOFs 转为 Bernstein controls，再用 de Casteljau/差分 controls 求值和导数；修正也用已知 Bernstein controls [1,1,2/3,0,0] 与 [0,0,8/3,0,0]，不调用生产 basis/Horner。原 48 个 XY/RZ、uniform/2:1、n4/8/16、full/half/hole/disconnected 布局及三种 wall 选择共 144 组，覆盖零、允许 P2、非多项式 native 输入。P2 只在原信息可辨识范围声明精确，缺秩报告保留。

在真实子面及共同支撑分片上检查两个 incident cell 的完整四次值/一阶导数迹：将五个固定节点上的差转换为 Bernstein 系数界，节点数来自解析次数，不把任意场的几个样点称为连续性证明。另核验真实 native 均值、atom begin/end/内部点与独立积分、piece 上 GL3×GL3 Gauss 积分、轴/指定 wall、极短 width 和高次非零贡献。精度预算固定为 r3/r4 的 1e-11，不强制旧点态常散度。记录场更新/查询成本、分片与存储库存，原接受路径回归必须保持。

永久原契约只保留 n4 的 48 结构组与 32 个身份/库存/区间/FP 拒绝；同一检查函数的 n4/8/16 完整 144 组由临时入口执行，源码和 runtime 捕获后退役。完整矩阵不是新的默认长 benchmark，临时入口不留在交付代码中。

本任务只交付共同数学内核；来源部署须同步改实际 B/J/atom、piece/overlay、完整次数 C1 运动证书和来源专属参考，再复验 native/affine/H10b-v2。不能因局部查询 PASS 宣布 H08a3 或 H08a 已完成。
