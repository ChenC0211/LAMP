# I09 / 2:1 coarse–fine reflux

## 任务 ID / 标题

I09：2:1 reflux；唯一验收点 `coarse flux = sum(fine subflux)`。前序 I08 的归档及 10 项交付哈希已核验；上一轮为有效进展。

## 目标

让每个完整 leaf side 消费其实际 native 子面的唯一最终通量，粗面严格取两个细子面之和；相库存更新、粗面视图与细面视图保持同一 authority。

## 允许修改

- 新 `include/Modules/Interface/HermiteGalsReflux.h`、`src/Modules/Interface/HermiteGalsReflux.cpp`。
- `include/Modules/Interface/HermiteGalsPhaseFlux.h`、`src/Modules/Interface/HermiteGalsPhaseFlux.cpp`：仅添加 Result 的 mesh identity 校验入口并更新注释，不改变 I08 积分、修复和更新公式。
- 新 `tests/interface_hermite_gals_reflux_contract.cpp`、`tests/CMakeLists.txt`。
- 本 TODO、API、报告、README 和 `runs/baselines/I09_20261001/` 证据。预算 4 production + 2 tests。

## 禁止修改

不改 Core mesh/几何、I01–I08 数值公式、Hydro/Electricity/App/reader/case、三 PDF 与 S08 冻结实现。不做 I10 bubble、I11 发布或 I12 完整输运；不引入 subcycling/动态 regrid、不增加第二套面积分或几何。

## 数学约束

1. 冻结方法式 (34) 与 I-R3：`Q_coarse = Q_fine,0 + Q_fine,1`，同时适用于总通量、液相通量和两个方向的 gross 通量。输入只来自 I08 最终共享面条目，包含微修复；不使用 raw diagnostics。
2. native 面以 +x/+y 定向；cell side 外法向 west/south 为 -，east/north 为 +。coarse outward 与两 fine outward 成对抵消。
3. 已缓存的 `liquidVolume = dt*Q_l` 逐子面相加，避免对粗面再次乘 dt 改变舍入。粗面没有独立通量未知量。禁止对 M clipping 或再次施加 reflux correction：I08 已按实际子面更新两侧 M。
4. 每 cell 四个完整 side，一侧只能有 1 或 2 个实际面；由既有 FaceSegmentGraph 原生拓扑归组，不由浮点坐标近邻搜索建立另一套关联。
5. 只读 Result 持有原 FluxResult，faces/M/geometry/time/identity 都保留同一所有权。不同 mesh、moved-from 等明确拒绝。Report 只诊断，不参与数值分支；内核不读写文件。

## 实现要求

提供按 cell/side 排列的完整 side 列表、其子面索引和汇总；统一/细侧一个条目、粗侧两个条目。校验子面 incidence、方向、2:1 邻层及原生身份。保留同一 FluxResult 的 updated M，只暴露一个库存 owner。报告粗细侧数量及定向 reflux residual。

## 验收

- **T0：** Release 主程序、新 contract、architecture、前序 Interface contracts 编译/运行通过；无新增 warning。I08 数值函数文本保持不变。
- **T1 拓扑/恒等式：** XY/RZ、uniform/2:1 多层网格、四方向、正反流、同面 counterflow、纯相/静止/axis；每实际面恰按真实邻接入组。粗面每个数值严格等于两个最终子面之和。coarse outward + fine outward residual `<=1e-13*scale`。
- **T1 解析：** cut 粗面通量对独立解析原函数，误差 `<=3e-11`；验证 RZ 真度量，不能面积平均代替积分和。修复 fixture 必须在 coarse–fine 子面触发，粗面使用 final 而非 raw 通量；M 数组保持原 owner/值。
- **T1 ownership：** 源 handle 销毁、只读 API、陈旧 mesh、moved-from 拒绝；失败不改输入。
- **T2：** 公开 I08/I09 API 多步解析平移跨 AMR 接口，逐 cell M 误差 `<=3e-11`，全局 ledger `<=1e-12`，库存严格有界。明确 analytic midpoint 供给只是模块算子验证，不冒充 I12 自洽 jet/geometry 时间推进。
- **回归：** I01–I08 与旧 native state/GALS/parametric control；I08 完整数值输出不变，69 个冻结源文件与 3 PDF 哈希不变。日志、早期失败、源码/runtime 身份归档。

## 失败规则

恒等式/解析/ownership 失败记 FAIL 并在本范围内修复；必要接口缺失、需偏离冻结方法或超预算记 BLOCKED，说明证据。不放宽验收，不以本任务替代全项目完成。

## 交付

实现/测试、API/TODO/报告、原始指标/构建与失败日志、源码/runtime 哈希及逐项验证归档。全部验收后才记 I09 PASS；下一项 I10，总体 goal 保持 active。

## 执行结论

**PASS，2026-10-01。** 4 production + 2 tests，13 项 CTest 通过，24 组解析 flux 最大误差 `1.11e-16`，三类定向 reflux residual 为 0；4 个粗面读取已修复的细面最终通量。前序输出、69 个冻结文件及 3 PDF 不变。

多步精确传递验收使用整面穿越发生于步端的 8 步序列；早期一般 dt 穿越显示 midpoint Heaviside 时间求积不精确，该证据和几何零迹拒绝保留在 [报告](I09_Coarse_Fine_Reflux_Report.md)，未宣称任意穿越时刻稳定/精确，也未替代 I12。I10/I11/I12 与总体目标继续未完成。
