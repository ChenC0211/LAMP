# H08a2r7 / 完整表示的实施前设计

2026-10-05。上一轮为 progress，H08a2r6 冻结 141 payload、12 当前交付文件及 1381 源码/输入已核验。检查归 `runs/checks/H08a2r7/20261005_145304_426`。

本轮 production 4：`include/Modules/Hydrodynamics/CchSpatialAdvector.h`、`src/Modules/Hydrodynamics/CchSpatialAdvector.cpp`、`src/Modules/Hydrodynamics/CchPhaseGeometry.cpp`、新同模块声明 `include/Modules/Hydrodynamics/Detail/CchNormalRatePolynomial.h`。test 2：现有 phase/momentum geometry contract。无 CMake、benchmark、材料或界面代码修改。

`CchNormalPolynomial` begin/end 数组增加至五项。值用完整 Horner；积分先计算 endpoint chart 在 offset 处的移位系数 b_k=Σ_{j≥k} binom(j,k) a_j offset^(j−k)，再按独立 width 积分 Σ b_k width^(k+1)/(k+1)。RZ 额外为 2π[r(offset)·原积分 + chartDirection·Σ b_k width^(k+2)/(k+2)]。乘法/求和保留低次原顺序；高次不靠 rounded endpoint 相减，不截断小项。

原 `ratePolynomial` 提取成 `Detail::normalRatePolynomial`，仅增加同模块函数声明，以便原 phase contract 调用真正生产转换。继续使用原 Core expansion primitive、长度/测度分离项和原系数表示，不另建根或几何 owner。遍历完整数组；每个 normal 的 k 次项对应 unit density 的 k 次，径向切向再包含 k+1 次。生产 prepareSpatial 调用这一个函数。

phase GL3 参考（对四次×线性半径精确）与绝对系数尺度、线性检查覆盖全部项；momentum 的固定三项零数组/非恒定检查扩展为真实库存。原 within-piece 二次运动导数测量属于当前未切换 B，保留，不把它误作新正则 B 证书。

两份原契约追加独立检查：单项式 0–4 和混合四次，两种长度、XY/RZ 固定/变化半径含轴、双端 chart、整段/子段、正负 normal；已知四重根的极短独立 width（2^-60/2^-80）；生产转换的完整次数、逐点值、Core range/positive-range 实际消费。参考以闭式幂积分或 GL3 取得，不复制生产转换。输出另做高精度重算。精度门槛事先固定 1e-12，极短非零矩按其真实量级验收。

标准 Release 构建六份受影响 Hydro 契约与既有 benchmark；架构/目录及六契约完整运行。以原冻结输出比较确定性行（仅原计时另列）；既有三个导出模式完整复跑并按列名比较，新增 c3/c4/e3/e4 在实际旧二次场中应严格为零。实际物理列的任何舍入差必须定量解释并通过原 gate，禁止为保持输出而增加次数 fallback。H10o/H10b-v2 完整窗口在实际 B 切换后必跑，本表示准备不代替该项。
