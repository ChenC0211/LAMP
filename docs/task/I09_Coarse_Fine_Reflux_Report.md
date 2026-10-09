# I09 / 2:1 coarse–fine reflux 交付报告

**状态：PASS。日期：2026-10-01。下一项：I10。** 本次只验收粗面通量等于两个最终细子面通量之和；总体 goal 与 Gate I 未完成。

## 改动与数学关系

新增只读 `HermiteGalsRefluxResult` 与 `assembleHermiteGalsReflux(mesh, flux)`。原生 `FaceSegmentGraph` 提供 cell side 归组；粗侧两子面、细侧单子面共享 I08 的同一 final face 条目。总/液相 signed rate、正反 gross rate、已缓存液相输运量与求积估计直接相加。四方向外法向对应明确，检查实际 face incidence、2:1 层级和双向配对。

I08 本来已按实际子面更新两侧库存，故 I09 保留其同一 faces/M owner，没有第二次库存 correction。原通量微修复后的值自动进入粗面，raw 通量只留在 I08 诊断。没有新增积分、第二套 scalar geometry、独立粗面通量或 subcycling 时钟。Report 在结构检查后生成，不参与数值分支。

预算内共 **4 production + 2 tests**：

- 新 `include/Modules/Interface/HermiteGalsReflux.h`、`src/Modules/Interface/HermiteGalsReflux.cpp`。
- `HermiteGalsPhaseFlux.h/.cpp` 只增加 Result 的 `validateMesh` 并更新注释。移除这一新增方法后，I08 cpp 与归档文本完全相同（统一换行后比较）。
- 新 `tests/interface_hermite_gals_reflux_contract.cpp`、注册 `tests/CMakeLists.txt`。

API 见 [接口说明](../api/I09_Hermite_Gals_Reflux_API.md)，任务边界见 [TODO](I09_Coarse_Fine_Reflux_TODO.md)。

## 验收证据

全部原始证据：`runs/baselines/I09_20261001/`。

| 项目 | 结果 |
| --- | --- |
| T0 Release | 主程序、I01–I09/旧 Interface state、architecture 全部构建；无 compiler warning |
| CTest | 13/13 PASS，18.94 s；I09 contract 5.28 s |
| 24 组独立解析完整侧通量 | XY/RZ × uniform/level-1/level-2 AMR × 两类方向场 × 正反速度，最大误差 `1.1102230246251565e-16` |
| 网格与四方向 | 64/112/364 leaves；对应 0/16/56 coarse sides；多层 AMR 同时覆盖 0→1 与 1→2 的 2:1 交界 |
| exact 子面和 | 每侧 8 个数值字段严格等于最终子面和；独立扫描 actual face owners 确认索引、方向与 incidence |
| reflux 残差 | 总速率、液相速率、液相体积的 coarse outward + fine outward 最大值均为 0 |
| counterflow/RZ | 液相正反 gross 均非零；RZ 轴面恒零；径向/轴向无散场及 `2*pi*r` 解析积分覆盖 |
| 纯相与静止 | 12 组，纯相通量恒等式成立，M 严格不变 |
| 多步 AMR 穿越 | XY/RZ × 正反方向，每组 8 步；XY 逐 cell 误差 0，RZ 最大 `2.7755575615628914e-17`；全局 closure 最大 `1.6653345369377348e-16` |
| coarse–fine 微修复 | 12 actual faces 被修复，涉及 4 coarse sides；粗面与 raw 和最大差 `3.4694469519536142e-18`，使用 final 和且 closure=0 |
| ownership / 拒绝 | faces/M/snapshot 地址与源一致；源 handle 释放后仍有效；5 项 moved-from/陈旧或异网格拒绝通过 |

解析 band 使用 Q2（XY 或轴向）/径向偶 Q4（RZ）；轴向 counterflow 速度为 `x-a` 或 `r^2-a^2`，另有 XY 横流/RZ `(u_r,u_z)=(r,-2z)` 及其反向。解析测通量 fixture 的库存为半容量，验证通量算子；不宣称其几何与 M 已经 I10/I11 闭合。真正多步测试携带上一更新结果 M，重新提供解析 endpoint/midpoint geometry，验证跨 AMR 的守恒传递。

数值数据见 `analytic.csv`、`translation.csv`、`repair.csv`；完整输出在 `reflux_final.log` 与 `ctest_details.log`。未放宽数值误差阈值。

## 回归与冻结

- I08 完整 contract 输出与 I08 最终通过日志逐字相同。
- 旧 native state、I06a skeleton、I07 volume 完整输出相同；I03/I04 共 23 行数值相同。
- 旧 GALS 与 parametric state control 退出 0，完整输出与 I08 相同。GALS 日志 SHA256 为 `C6D8D43979F6A5998E9C54B292EA7EE38EA627A6352E66F9CB0BDBA779AAA5D1`。
- 基线 1568 个源/输入文件中 1565 不变、3 个既有文件变更，另新增 3 文件；未发现越界变更。69 个 S08 专用冻结文件、3 份方法 PDF 哈希不变。
- 本次增量构建没有 warning；之前冻结 Hydro 的两处 C4244 本轮未重新编译，不能据此称其已经消除。
- 未接线 App/Hydro/字典，未改变 baseline case。S01 `hydro_moment_xy` 已知首步失败没有被重分类为成功。

## 早期失败与验证边界

所有早期失败日志保留，`iteration_history.json` 列明原因：

1. 初始 cut-coarse 覆盖检测只看 signed net，在 counterflow 净和为零时误判；改为检查 gross phase throughput。
2. 初始 level-2 refine box 缺少 level-1 缓冲，得到非 2:1 网格并被既有 Core 拒绝；测试改为嵌套两层盒，未修改网格实现。
3. 初始 dyadic band 根与某条细面重合，触发已知正长度零迹拒绝；解析通量矩阵改用非对齐 band。此 Core 限制未解决，也未被包装为支持。
4. 初始多步平移以一般 dt 穿过整张面，却错误要求 midpoint Heaviside 时间积分精确，观察到约 `3.9e-5` 的 cell 误差。最终 reflux 精确传递测试使用步端对齐的穿越，实际跨越 coarse–fine 交界及后续细单元。该结果**不证明 off-step 穿越的时间精度/稳定性**；一般穿越须在 I12 以时间收敛与严格有界性继续验收，不能把本测试替代 I12。
5. 清理测试间接 include 后缺少 Snapshot 完整类型，full build 失败；补直接 Core include 后重建、全部 CTest 通过。生产数值实现从初次构建起未因这些测试失败改变。

本任务不引入 bubble、同源几何 publisher、动态 AMR/regrid/subcycling 或跨模块耦合。半面内部精确相切/重根、正长度零迹支持仍是 I11/I12 前的待解几何问题；任意时间跨面输运、长时间稳定性、曲率、低黏/无黏与电喷雾收敛尚未证明。下一项严格按路线图为 **I10 interior bubble scalar closure**。

## 交付与归档

`verification.json` 汇总任务证据；`delivery_manifest.json` 给出交付文件与 archive SHA256；`evidence_inventory.csv` 为 payload 清单，`archive_verification.json` 校验全部 archive entry、工作区源码副本与 runtime 副本。含构建、13 项测试、回归、早期失败、源码审计和程序身份；Git commit 为 `unavailable`。
