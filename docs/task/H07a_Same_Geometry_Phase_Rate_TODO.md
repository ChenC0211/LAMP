# H07a / 同源 native 相测度与瞬时体积流率

状态：PASS，2026-10-03，见[报告](H07a_Same_Geometry_Phase_Rate_Report.md)。前置 H01–H06 PASS；H07 拆为本项和 H07b（真实分相体积求积与 H06 黏性接入）。本项 PASS 不关闭 H07。依据[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md)、冻结 CCH §4.1–4.2 和用户本轮授权的[完全非守恒材料耦合设计](H07_Pure_Gals_Hydro_Coupling_Design.md)。无 Interface/Hydro 守恒材料库存。

## 目标

只解决一个数学事实：CCH 原生单元/面的液气测度及瞬时相体积流率，来自唯一 FirstJet publication 的同一 geometry identity 和时间层，不再按单元中心判相。

## 允许修改

- 新增 `include/Modules/Hydrodynamics/CchPhaseGeometry.h`、`src/Modules/Hydrodynamics/CchPhaseGeometry.cpp`；只读消费正式 Interface publication。
- root `CMakeLists.txt`：显式 source 与正式 `quantum_interface` 依赖。
- `tests/contracts/hydrodynamics/cch_phase_geometry_contract.cpp`、`tests/CMakeLists.txt`。
- 本 TODO、耦合设计、API、报告、现行路线图/规范/索引；证据 `runs/baselines/H07a_20261003`，临时工作 `runs/checks/H07a/<时间戳>`。

预算：3 production/build + 2 tests。

## 禁止修改

H01–H06 实现/测试、Core、Interface、Solvers、旧 Hydro/Compatibility、Electricity、App、cases、冻结 PDF、历史报告。禁止新增 Interface M、体积校正、平滑 Heaviside、中心判相、顶点平均、第二 geometry 或 D_I、epsilon/clipping、失败 fallback。H07b 分相体积求积、H08 时间积分/质量更新、H09 dual transfer、H14 时间推进不在本项中冒称完成。

## 数学约束

- 当前程序液相为 `phi<0`；只读消费 `V_l(phi)`，`alpha=V_l/V_c`，保留同源真 XY/RZ 容积。几何质量不是守恒库存。
- 每个真实 native 子面仅消费 publication 保留的两条 vertex–midpoint trace。分别积分液/气正长度区间；RZ 使用同一 `2*pi*r`。短相区间用独立 offset/length 积分，不能相减两个舍入端点。
- `q_l=A_l*u_f`、`q_g=A_g*u_f`、`q=A_f*u_f` 是瞬时体积流率；正号沿 native 正坐标法向。液气分别积分，不能用 `A_f-A_l` 抹掉极小气相。
- 共享面只算一次；逐单元外流由保留的同一 D_I 的符号生成，内部抵消和真实 2:1 子面求和成立。零面积轴面流率严格零，轴上径向速度必须严格零。
- geometry/layout/actual mesh/measure/完整时间层一致，不能仅比较 time() 或尺寸；不同候选或 regrid 必须显式重建。结果不可变并保留来源寿命。
- 有限性、合法相范围、分割闭合不满足即拒绝；不得修正几何或流率。分割闭合允许既有浮点表示误差（相对真实 face measure 不超过 `128*epsilon(double)`），不剪裁相区间。

## 实现要求

- [x] TODO/设计先于源代码；冻结前置交付与源/输入/PDF 身份。
- [x] typed immutable geometry 和 rate 对象，来自实际 publication、真测度与 D_I；无 IO、时钟 owner 或 accepted commit。
- [x] 显式来源与时间校验，复制/寿命/regrid/selected domain 处理。
- [x] 原生 XY/RZ、uniform/2:1、纯相/小切片、双根曲面解析测度与共享流率检验。
- [x] 标准 build/CTest、前置/App/架构回归、源界限/数值/交付哈希审计；闭合本项，H07b 为下一项。

## 验收

**T0：** local Release 无新增 warning，唯一 CMake owner；依赖闭包可含正式 Interface/Core/linalg，绝不含旧 Hydro/Compatibility/App；architecture/selftest/repository_layout/lamp_help 通过。

**T1：** 解析平面和圆/球的实际单元/面测度、瞬时流率归一化误差 ≤1e-10；纯相精确零/全相（浮点测度容差）；短液/气 face interval 仍为正。面液气和误差、逐单元与边界总外流闭合 ≤1e-12。正反速度、逐面符号及真实 2:1 聚合通过；无守恒 Interface M。错 layout/geometry/time/actual mesh/measure、历史 ConservativeVolume、缺几何覆盖、非法 native 值/溢出/轴值明确拒绝，失败不改变既有对象。

**T2：** 从同一非多项式解析曲线的三个 n=4/8/16 Hermite publication 测量材料分割/瞬时相率误差，XY/RZ uniform/2:1 递减且末档阶 ≥2；独立解析参考不得参与计算。此算子验证不等于两相 NS 或质量守恒验收。

## 失败规则

同源/测度/小相/闭合/收敛失败即 FAIL，在允许范围修复，不放宽门槛。若需修改 Interface owner 或超预算，先拆缺失 contract，H07 保持未关闭。实际外部依赖无法解决时才 BLOCKED；不能因为数值工作困难而阻塞。

## 交付

生产同源测度/瞬时 rate、短契约、API/TODO/报告、原始日志及哈希证据。明确本项无时间积分、质量更新、分相体积求积或 CCH App；总目标保持 active。
