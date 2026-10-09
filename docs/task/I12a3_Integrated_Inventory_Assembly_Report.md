# I12a3 / integrated volume 库存装配报告

**FAIL。** 已建立直接消费同源 integrated face volumes 的库存结果，旧 I08 与新区间路径共享同一个私有装配器；但真实 RZ 跨面矩阵发现原有单向舍入修复的可达性缺陷。本报告不把代数检查或旧回归通过等同于新装配通过。按用户持续推进授权，下一项先修复该缺陷，不接 reflux/driver。

**后续状态：** 上述为初次交付的历史结论，失败 archive 保持不变。独立 [I12a3r1](I12a3r1_Bidirectional_Roundoff_Repair_Report.md) 已在原预算内修复该缺陷并通过完整新矩阵；I12a3 当前验收 PASS，I12 整体仍 FAIL。

## 变更

4 个 production 文件：`HermiteGalsTemporalPhaseFlux.h/.cpp`、`HermiteGalsPhaseFlux.cpp`、新 `Detail/HermiteGalsInventoryAssembly.h`。2 个测试相关文件：新 `interface_integrated_inventory_contract.cpp`、CMake 注册。见 [TODO](I12a3_Integrated_Inventory_Assembly_TODO.md) 与 [API](../api/I12a3_Integrated_Inventory_API.md)。

新结果保留 interval、原始/最终共享面量、库存与 ledger；不伪造 midpoint，不重算几何库存。只接受旧时刻起始、完整有序且同源的面量。提取公共代码时保留原 rate 分支的 dt 乘法与加减次序，新分支直接使用体积。未修改旧修复规则、预算、CFL 或任何其他模块。

## 实际验证

- 主程序及 21 个相关测试目标 Release 构建通过，29.74 秒，无 warning。
- CTest 共 20 项，**19 PASS / 1 FAIL**，108.55 秒，真实退出 8。失败项是新增 integrated inventory contract，未设置 WILL_FAIL。
- 前序 19 项全部通过；除 architecture 外，18 项完整输出与 I12a2 相同。architecture 为 12 targets / 379 translation units / 897 files / 3078 includes，5 条 legacy edges 保持。
- 新跨越矩阵完整通过 XY uniform/AMR 正反速度的 68 组；RZ uniform 反向相位 0–11/16 通过，在 12/16 首次失败。每组消费实际 I12a2 内部面量，外边界只提供本制造解的明确解析历史。末态参考仅用于比较，不写入库存。
- 独立 `--algebra-contract` 通过 2 个数值检查与 20 个拒绝检查，覆盖 counterflow、原内部共享舍入修复、来源/面完整性/区间/mesh/配置/CFL/散度/超预算库存等；输入和旧状态不变。
- 原 `--crossing-probe` 仍退出 1，完整输出与 I12a2 相同；完整 I12 driver 未接新装配，也未重复长程运行。
- 69 冻结文件、12 集成参考、3 PDF、I12a2 archive 身份一致。首次编译发现测试误把 `phaseVolumes` 字段写成函数，修正后构建通过；日志保留。

## 失败的数值事实

4×4 RZ、u_z=-1、dt=1/64、q=12/16，原始 cell 11：

| 量 | 值 |
| --- | --- |
| M_old | 0.327504898213584417 |
| 共享面净输入 | 0.0161067982727992332 |
| capacity | 0.343611696486383622 |
| 超容 | 5.5511151231257827e-17 |
| 两个轴向面总量 | -0.0214757310303989764 |
| 穿越面液量 | -0.00536893275759974409 |

两面 total 与 `native face.area * dt` 完全相同，排除了本例度量重算差异。超容在原 64*eps*capacity 预算内。旧修复只减小入流，会使上游原本满液单元超容；缺陷继续传向满液入流边界，最终报 `no interior flux capacity for roundoff repair`。与此同时，界面侧的内部出流仍小于 total，存在不改边界、M 或总通量的局部有界调整余量。

这要求单独扩展内部共享通量的舍入调整能力，不能在公共装配提取中顺手更改数学规则，也不能放宽门槛或改初值掩盖失败。下一独立修复须保留原始量、双方反号交换、gross subset bounds、64*eps 面修正预算、严格 M 容量和完整 ledger。

## 复现与归档

```powershell
build/windows/tests/Release/quantum_interface_integrated_inventory_contract.exe
build/windows/tests/Release/quantum_interface_integrated_inventory_contract.exe --algebra-contract
```

证据在 `runs/baselines/I12a3_20261001/`。`verification.json` 明确记录任务 FAIL；source/archive audit 的 PASS 仅指来源、回归与归档完整性。`evidence.zip` 和 `archive_verification.json` 固定此失败版本，后续修复不覆盖它。
