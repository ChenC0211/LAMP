# H07b / 同源正权分相体积求积与 CCH 黏性

状态：**PASS**，2026-10-03；验收见[报告](H07b_Phase_Quadrature_Viscosity_Report.md)。前置 H01–H06、[H07a](H07a_Same_Geometry_Phase_Rate_Report.md) PASS；依据[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md)、H-R5 和用户授权的[完全非守恒材料设计](H07_Pure_Gals_Hydro_Coupling_Design.md)。H07a/H07b 均已验收，H07 关闭；下一项 H08，不代表 Gate H。

## 目标

只解决一个数学事实：由唯一 FirstJet Q4 得到位于实际液气区域内部的正权物理体积求积，使两相黏性通过 H06 的同一 E^T H_mu E 分相积分并绑定 geometry/time/native 来源。

## 允许修改

- 新增 `include/Modules/Interface/HermiteGalsPhaseQuadrature.h`、`src/Modules/Interface/HermiteGalsPhaseQuadrature.cpp`：正式只读几何 owner 提供求积。
- 新增 `include/Modules/Hydrodynamics/CchPhaseViscosity.h`、`src/Modules/Hydrodynamics/CchPhaseViscosity.cpp`：typed CCH 消费并接入 H06。
- root `CMakeLists.txt`：新 CCH source 显式 owner；Interface 沿用现有模块 source 注册。
- `tests/contracts/hydrodynamics/cch_phase_viscosity_contract.cpp`、`tests/CMakeLists.txt`。
- 本 TODO、API、报告、现行路线图/索引；证据 `runs/baselines/H07b_20261003`，临时工作 `runs/checks/H07b/<时间戳>`。

预算：5 production/build + 2 tests；实现前超过预算须拆任务。

## 禁止修改

H01–H07a 实现/测试、现有 Interface 状态/输运/测度/发布、Core/Solvers、旧 Hydro/Compatibility、电场/App、cases、冻结 PDF、历史报告。禁止新 geometry、材料库存、体积拟合、剪裁/权重缩放、中心/顶点判相、按 alpha 混合全单元点黏度、人工最小黏度、第二 D_I/C、失败 fallback。H08/H09 非守恒材料缺陷与动量更新、时间推进及完整生产边界留待对应任务。

## 数学约束

- 唯一 immutable publication，`phi<0` 液相；所有点、区间、容积、时间层来自同一实际 Q4/native mesh。无 accepted commit 或主时钟。
- 每条 x/r 外层积分线用正式 MappedSegment 完整液气区间；GL4 内层 y/z 正权积分、GL8 外层自适应积分，RZ 权重含同一 `2*pi*r`。
- 短区间用独立 offset/length 计算，不用两舍入端点相减；点须有同一 scalar 的严格相证书，未分辨点/区间拒绝。
- 完整 topology 图端点/样本投影仅用于初始积分分段，不能用于判相或替代原 scalar。显式 bounded subdivision 控制缩放至原单元坐标的各相四次多项式矩；非轴 RZ 另控制 (r_c/r) 和 (r_c/r)^2 权重矩，覆盖 H06 局部 P2 hoop 能量。
- GL8 与两个子区间 GL8 的差为误差估计，不冒称严格误差上界。相容积与 publication 的独立 Green 测量检查为 `1e-8 V_c`，正的两相必须各自实际有正权节点；不能把容积差用于修改权重/界面。默认矩估计阈值 `1e-9 V_c`，初始分段按宽度分配局部预算，并核验所有接受子段估计的累计和；有限工作量/深度超限明确拒绝。
- 各相 mu≥0，H_mu=2*mu_phase*dV≥0；K_mu=E^T H_mu E 对称半正定，零黏矩阵/载荷/耗散精确零，无黏度 floor。CCH 同一原生速度支撑/轴奇偶/H05/D_I 均保留。
- typed prepared 对象保留 geometry/layout/time/material 的完整来源；选中 CCH 子域覆盖齐全，不允许给求积换候选标签或静默重用 regrid 前数据。

## 实现要求

- [x] TODO 先于源码；复核 H07a 交付并冻结源/输入/PDF/程序身份。
- [x] Interface 正式正权分相体积求积、有限误差/工作量校验、不可变寿命与完整来源。
- [x] CCH 正式相黏性适配复用 H06 E/K；保留材料值与实际相标签，无平均判相。
- [x] XY/RZ、uniform/2:1、纯相/平面/圆球/小 cut、非多项式几何，独立容积/矩/耗散检验及拒绝事务。
- [x] 标准 build/受影响 CTest、前置/App/架构回归、源范围/数值/文档/交付哈希审计。

## 验收

**T0：** local Release 无新增 warning/error；新 cpp 分别且仅归正式 Interface/CCH，真实依赖不进入旧 Hydro/Compatibility/App；架构、自检、目录、help 通过。公开头自足，无 IO/字典/时钟 kernel。

**T1：** 所有点权重>0且实际相正确；解析平面/圆球体积与二阶矩（原单元缩放坐标）相对 V_c 误差≤2e-8；小液/气相不能消失，解析平面短相体积相对自身误差≤1e-5。P2 native 速度的两相积分耗散与独立分相解析参考相对误差≤5e-8；u^TKu 与逐点积分误差≤1e-10、对称相对缺陷≤1e-12、最小特征值≥-1e-10*max(1,谱半径)。mu_l/mu_g=0、1、1e-6、1e-12 的合法组合连续缩放，无黏精确零；实际相 K 进入 H04 同源鞍点参考求解。错 geometry/time/layout/mesh/measure、负/非有限材料、非法积分设置/工作量、极端溢出/轴值明确拒绝，既有对象保持。

**T2：** 两相 mu_l=1、mu_g=.01，XY u=sin(pi*x)cos(pi*y), v=cos(pi*x)sin(pi*y)，RZ u_r=r*sin(pi*z), u_z=cos(pi*r)sin(pi*z)。n=4/8/16 uniform/2:1 的实际分相物理耗散与独立解析场/高阶体积分参考误差递减且末档阶≥1.8；参考不进入生产接受条件。只验证算子，不声称两相 NS/低黏动态稳定性或时间阶。

## 失败规则

相节点/正权/切分精度/同源/PSD/无黏/收敛不通过即 FAIL，在允许范围修复，不放宽 gate。需要改 Interface 主状态/核心几何或超预算时先拆缺失 contract，H07 保持未关闭；外部依赖真正无法解决时才 BLOCKED。不能用 alpha 混合点或原单相矩阵恒等式替代真实分相积分。

## 交付

正式同源分相求积与 typed phase viscosity、短契约、API/TODO/报告、原始数值/回归及哈希证据；列明 H08/H09、CCH App、完整 PDE/时间推进和低黏 EHD 未完成。总目标保持 active。
