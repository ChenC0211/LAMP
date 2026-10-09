# S08 / 冻结三条实验基线

任务依据：[05 分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md) S08；前序 S01、S07 的原始交付和当前文件已核对。执行日期：2026-10-01。

## 目标

只解决一个架构事实：`parametricGals`、大型 `CompatibleMoment/cchMomentMidpoint`、`sharpHermiteP3` 是冻结的实验 baseline，只保留 A/B、旧 case 回归与比较用途，不再承载 V2 新功能。

## 允许修改

以下五个 production 文件只增加 `S08_FROZEN_BASELINE` 注释块，不改变原有代码：

1. `include/Modules/Interface/InterfaceParametricProjection.h`
2. `include/Modules/Interface/InterfaceParametricTransport.h`
3. `include/Modules/Hydrodynamics/HydroCompatibleMomentFlow.h`
4. `include/Modules/Hydrodynamics/HydroCompatibleMomentSpline.h`
5. `src/Modules/Electricity/SharpDielectricElectrostatics.cpp`

文档：本 TODO、`S08_Experimental_Baselines_Report.md`、`S08_Frozen_Baselines.json`、`docs/README.md` 中的本任务索引。快照、构建和验证证据保存于 `runs/baselines/S08_20261001/`；允许正常构建更新 `build/` 和 `bin/`。

## 禁止修改

- 离散算子、运行时分支、接口声明/布局、reader、默认方法、字典和容限。
- 原始 `cases/` 输入、现有测试、共享 Core/Solvers 实现、其他模块及旧运行结果。
- 新建 I01 API、接入 V2 方法名、删除旧方法或加入隐式 fallback。
- 修改 S01/S07 历史交付或把其已知失败记成通过。

## 数学约束

本任务不改变任何数学事实。Interface V2 必须是 first-jet skeleton + 独立 cell phase volume；Hydro V2 必须是 native face/cell、`C=-D_I^T` 与局部 compatible Hodge；Electricity V2 必须另建 SI-HHO-Q。冻结旧路径不豁免 I/H/E gate，也不冻结将来 V2 所需的通用 Core/Solvers。

## 实现要求

- [x] 核对前序证据及三条方法的 reader、dispatch、状态 owner 和实现族。
- [x] 建立显式冻结清单：方法名、V2 目标、标记点、路径范围及 SHA-256；区分专有实现与共享调用入口。
- [x] 在五处代码 owner/实现入口增加明确的 baseline 注释，禁止 V2 扩展与隐式回退，指向冻结清单。
- [x] 证明除新增注释外原代码保持原样；其他源码、测试、输入、构建配置与前序归档不变。
- [x] 完成 Release 构建及现有 Interface/Hydro/Electricity 短验证，原样记录已有失败。
- [x] 逐条复核 S01–S08 的证据，记录 Gate S 是否满足；本轮不开始 I01。

## 验收

- **T0：** 主程序及选定旧方法 contract 目标构建成功；五个标记及清单路径实际存在；不新增 compiler warning；reader/dispatch 不变。
- **T1：** 移除新增注释块后，五个文件的原内容与修改前快照一致（允许统一换行）；其它 production/test/input/build 配置哈希不变；冻结清单覆盖实际枚举的专有文件。
- **T2：** 复用 `quantum_interface_state_contract --parametric-state-control`、`quantum_hydro_module_contract --moment-time`、`quantum_electricity_evolution_contract`；明确是旧方法回归，不能据此宣称 V2 方法已通过。

## 失败规则

- 三条旧路径的 owner 无法确定，或冻结必须改变数值行为/扩大 production 修改范围时记 `BLOCKED`，先说明缺失接口或需拆分部分。
- 新增注释之外出现代码改动，或新增编译/回归失败由本次修改引入时记 `FAIL`，只修复本任务的非数值改动。
- 已有数值失败需记录并追溯；不得加 fallback、放宽容限或扩大本任务去修数学方法。

## 交付

五处代码注释、冻结清单、TODO/报告、构建/短验证/不变性证据与 Gate S 核验。预算：5 个 production 文件、0 个测试文件；只改变基线身份，不改变数值行为。

**最终结论：PASS。** Release 构建成功，三项预定旧路径验证均通过；5 处仅新增注释，其余 1,086 个源码/输入等文件未变。Gate S 通过，下一项为 I01；本轮未开始 I01。详情见 [S08 报告](S08_Experimental_Baselines_Report.md)。
