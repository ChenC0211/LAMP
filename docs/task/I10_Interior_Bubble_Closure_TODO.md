# I10 / interior bubble scalar closure

## 任务 ID / 标题

I10：单元内部 bubble 体积闭合。前序 I09 的归档与交付文件已核验，上一轮为有效进展。

## 目标

固定同一候选 first-jet/mixed skeleton，仅解每个 cell 的一个 bubble 系数，使同源几何体积等于独立守恒库存 M；保持单调性、求根括区与拓扑/幅值安全界。

## 允许修改

- 新 `include/Modules/Interface/HermiteGalsBubble.h`、`src/Modules/Interface/HermiteGalsBubble.cpp`。
- `include/Modules/Interface/InterfaceCellMeasure.h`、`src/Modules/Interface/InterfaceCellMeasure.cpp`：新增不依赖 legacy InterfaceState 的完整 Q4 体积及 center sensitivity 入口，复用已有 topology/integrator；暴露既有 quadrature estimates 的逐 cell 汇总以验收小 AMR 单元。
- 新 `tests/interface_hermite_gals_bubble_contract.cpp`、`tests/CMakeLists.txt`。
- TODO/API/报告、README 与 `runs/baselines/I10_20261001/`。预算 4 production + 2 tests。

## 禁止修改

Core 数学内核、mesh、I01–I09 数值公式、旧 bubble/legacy 生产路径、三份 PDF/S08 冻结专用源、Hydro/Electricity、Apps/reader/case。不得提前做 I11 accepted publication 或 I12 整体输运，不引入 clipping/KKT/旧算法 fallback。

## 数学约束

1. `phi = phi_boundary + c*q4(xi)*q4(eta)`，`q4(t)=16*t^2*(1-t)^2`。本入口固定以现有 Q4 cardinal basis 的 zero-center extension 定义 phi_boundary，故本入口 c 恰等于 center sample；这是明确的系数基准，不把任意旧中心值通称为 bubble coefficient。
2. bubble 及梯度在边界恒零；共享 master/hanging first jet、mixed derivative、普通/粗细 C1 trace、axis/contact trace 不变。c 只为派生闭合量，prediction 不新增独立 transported state。
3. `M(c)=integral H(-phi) dV` 为非增函数，`M'(c)=-integral_Gamma psi/|grad(phi)| dA<=0`。XY 使用面积，RZ 使用真 `2*pi*r` 度量。体积、导数、拓扑来自同一试算 snapshot。
4. 重用 `nativeQuarticCenterInterval` 的 exact Bernstein 零集管道偏导证书，在整个系数区间排除临界零/新增内部圈，并界定相对 cell 位移。它是充分安全区间，不冒充最大可达区间；无法证明时明确拒绝/请求细化。
5. 保持有符号体积括区。Newton 不在括区安全内部或斜率不可用时二分；迭代/存储分辨率耗尽拒绝。最终残差由独立 volume tolerance 检查，quadrature estimate 另有门槛。
6. 原始 M 永不修改；纯相的小容限内不匹配保留并报告，宏观纯相不匹配拒绝，禁止新造未标记 bubble component。超安全可达范围的 cut 库存拒绝，不扩大 correction 空间。

## 实现要求

输入 unpublished Candidate、同 identity 的 I06 mixed closure、带 revision 标识的预测系数及具名 Config。校验 mesh/identity/time/尺寸/有限性和 M 严格容量界。仅返回不可变 Closure/Report 及 revision 绑定，不发布第二份 production geometry。试算 snapshot 为私有临时量；多个独立 cell root 可批量评估，Report 不参与数值分支。用既有 topology discovery 完整发现所有分量；不允许用 caller phase mask 隐藏零集。

## 验收

- **T0：** Release 主程序、I01–I10、旧 native state、architecture 构建和回归；无新增 warning，4+2 文件预算。
- **T1 单调/导数：** XY/RZ，uniform/2:1 AMR，安全区间内系数扫描 M 非增；center sensitivity 非正，与中心差分在尺度归一后 `<=2e-5`。纯相极限与相反修正方向覆盖。
- **T1 求根：** 人工已知 admissible bubble 产生 target M，恢复后逐 cell 相对残差 `<=1e-10`；独立 adaptive volume estimate `<=1e-11`，最终 coefficient 保持证书区间，位移界 `<=0.125` cell。检查 Newton/二分 safeguards 和工作量上限拒绝。
- **T1 几何不变：** 边值/完整 first jet 不变；AMR master/hanging trace 与 RZ axis 保持；primary M/first jets、candidate identity 与未发布状态不变。
- **T1 失败隔离：** 不可达目标、纯相超界/隐含新分量、拓扑证书失败、错 revision/mesh、非有限/错尺寸 prediction/config、迭代耗尽均拒绝，源候选不变。
- **T2：** 公开 V2 API 在 XY/RZ 圆/球及 AMR 几何上完成混合正负 bubble 闭合；解析初始体积误差 `<=1e-10`（相对总容量），闭合 residual 独立检查。不宣称完整 self-consistent GALS 时间推进或真实 breakup。
- **回归：** I07/I08/I09 完整数值输出与前序一致；旧 GALS/parametric controls、69 个 S08 文件和 3 PDF 哈希保持。

## 失败规则

数值/ownership/拓扑 guard 不满足为 FAIL 并在范围内修复；缺必要数学接口、超预算或需改变冻结方法为 BLOCKED 并说明证据。不得放宽验收、隐藏异常或把本 task PASS 等同 Gate I/总目标完成。

## 交付

源码/测试、TODO/API/报告、原始指标与全部构建/失败日志、前后哈希及可验证归档。I10 验收完成后下一项才是 I11，总体 goal 保持 active。

## 执行结论

**PASS，2026-10-01。** 4 production + 2 tests；14 项 CTest、旧 GALS/parametric/bubble 控制通过。8 组 XY/RZ uniform/AMR 闭合最大相对 residual `9.6544e-11`，独立解析截面积分确认；Newton/二分/端点与 21 次拒绝均覆盖。原 M/skeleton 不变，既有数值及冻结哈希不变。

单 cell 小闭合圈的 Core discovery/guard 仍要求细化；这个拒绝限制与早期失败完整记录在 [报告](I10_Interior_Bubble_Closure_Report.md)。未发布生产 geometry、未证明真实 breakup 或低黏稳定性。下一项 I11，总体 goal 保持 active。
