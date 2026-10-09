# H10a / 共同位置中心对流与材料动能契约

状态：**PASS**，2026-10-03，见[报告](H10a_Compatible_Center_Convection_Report.md)。前置H01–H09 PASS；依据[路线图§0](../roadmap/05_LAMP_Module_Task_Roadmap.md)及实现前的[相容中心装配设计](H10_Compatible_Center_Convection_Design.md)。上一轮H09b是已核验的进展，本轮只做H10a；下一项H10b完整无黏benchmark，H10未关闭。

## 目标

只解决一个数学事实：唯一H09共享材料流量消费共同位置的单份affine速度重构，并通过负转置配对装配同时满足内部动量抵消、常速度一致性与显式非守恒材料动能预算。

## 允许修改

- 新增 `include/Modules/Hydrodynamics/CchKepConvection.h`、`src/Modules/Hydrodynamics/CchKepConvection.cpp`；root CMake唯一注册。
- `tests/contracts/hydrodynamics/cch_kep_convection_contract.cpp`、tests CMake注册。
- 本TODO、设计/API/报告、现行路线图/架构/索引/耦合设计；证据 `runs/baselines/H10a_20261003`，临时工作位于start.json登记的 `runs/checks/H10a/<timestamp>`。

预算3 production/build +2 tests。完整研究driver与benchmark输入留在H10b，不能把长研究脚本塞进短契约。

## 禁止修改

H01–H09实现/旧tests、Interface/Core/Solvers、旧Hydro/Compatibility/App/cases、冻结PDF及历史TODO/报告。无第二材料库存/流量、压力散度、全局速度未知量、隐藏ghost、人工黏度、upwind fallback或接受事务。H11耗散、H14corrector、H15一般Hodge/生产边界不在本项。

## 数学约束

- 保留同一不可变H09 owner、三层identity、完整中点clock和mesh/layout；重构为几何派生值，不重算F。
- 内部段仅一份共同位置t，affine再现；RZ实际测度中心。native采样中点与共同位置对齐的XY段保持算术中心t0；外边界半dual段依照真实几何判定，不将所有uniform段误当成对齐段。
- 载荷严格按设计的 `0.5*F*(b*(t^T*u)-t*(b^T*u)+diag(b)*u)`；公开raw中心通量与最终相容载荷的差异，不将其混为一项。
- 外部段显式同层trace，逐共享段一次计算；boundary momentum与KE功显式报告，轴径向trace0。
- `uR`消费H09原缺陷；半离散KE预算单列 `0.5*u^2*R`，无力power与边界KE功分别报告，不冒称全离散守恒。

## 实现要求

- [x] TODO与设计先于源码；核验H09b全部38payload并冻结1,037源/输入/PDF/App。
- [x] 不可变prepared、bound boundary、只读result及报告；finite/overflow/underflow和误用拒绝。
- [x] 单份真实同分量支持与affine约束；无IO/字典/主时钟，source唯一CCH owner。
- [x] 独立segment/reference载荷、常数/affine、随机功、真实材料非零缺陷与误用/regrid检查。
- [x] Release、新契约和前置回归、来源/依赖/链接/哈希审计与证据归档。

## 验收 T0/T1/T2

**T0：** local Release无新增warning/error；公开头自足，唯一CCH owner，正式Interface只读依赖，help/architecture/self-test/layout/App及前置流体契约通过；原输入和冻结证据不变。

**T1：** XY/RZ uniform/真实2:1、全域/半域/孔洞/不连通域。共同trace常数/affine归一化误差≤1e−12，随机独立装配与动量总收支≤1e−12；KE与材料预算≤1e−12（gross rates/powers标度，生产256eps）；中心项线性、复制只读、共同位置对齐的XY段与普通中心散射一致。实际coarse–fine owner段必须检出t0共同位置偏差>1e−3*h及raw独立功失配>1e−5归一化；配对后的KE及常速度同时过门槛。非有限边界/错误尺寸、foreign owner/layout/clock/identity、非法轴、overflow/underflow、坐标/拓扑改变明确拒绝。初版关于整个uniform mesh等于对齐段的几何前提已据实际半dual边界修正，详见设计；所有数值阈值保持。

**T2：** 正式FirstJet平面/圆球平移，XY/RZ uniform/2:1、密度比1/1e3/1e6；非零R和几何转移差异显式保留。实际载荷进入无外力常速度材料更新，速度恢复≤1e−10；一般非均匀中点速度的材料动能预算独立解析/代数检查≤1e−12。实际Interface refine/projection-coarsen后重建并继续clock，旧prepared拒绝；一般Hydro速度transfer不在本项。此处是短算子检查，不代替H10b完整时间窗/空间精度/无一阶系统耗散验收。

## 失败规则

同源、共同位置、动量/能量恒等式、常速度或误用拒绝不过即FAIL，在本范围修正；不得放宽门槛、改材料密度/物理时窗、加黏度或绕旧路径。需要改前置operator或超预算，先明确缺失接口并拆后续任务，H10保持未关闭；只有无可推进外部依赖才BLOCKED。

## 交付

中心重构/相容载荷、同层边界端口、材料动能预算、短契约及API/TODO/报告/原始数值与哈希证据。H10a通过后下一单任务H10b；总目标保持active。
