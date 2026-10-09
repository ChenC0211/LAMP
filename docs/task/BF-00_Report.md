# BF-00 / 方法切换与冻结基线报告

**用户决定 / 当前裁决：PASS（BF-00文档/冻结任务，档案前置已明确解除）。** 用户要求“不用管数据是否丢失，尽快切入新路线的实现”，明确允许用已冻结的当前源码/输入进入新路线。旧缺失、H13 FAIL、24/25回归及两次180秒timeout全部原样保留；此裁决不表示旧候选或BF-CSF物理通过，不豁免后续gate。下文BLOCKED为该决定前的历史检查点；下一唯一任务BF-01。

状态：**BLOCKED（原历史payload丢失）**。本轮仅BF-00，冻结与允许的检查已执行；0 production、0注册测试修改。当前25项基线回归24项通过，完整momentum-geometry原180秒门槛两次TIMEOUT，仍为未通过项；原H13保持FAIL，BF-CSF尚未实现/验收，H14/H15不启动。

用户已明确更新目标，以 [BF任务书](../roadmap/LAMP_BF_CSF_Codex_Tasks.md) 和 [BF-CSF V1](../roadmap/LAMP_BF_CSF_Method_V1.md) 为当前入口。按任务书先做BF-00，不继续旧H13切向模态方法开发；现有未采用代码及证据独立冻结，不能作为BF-CSF实现。

## 原数据缺失及状态依据

以下历史引用在当前workspace不存在；现存四份H13 evidence.zip的成员也未包含原20261006 long-wave payload。用户已明确答复“丢失”。

- `runs/benchmarks__cch_capillary/20261006_062652_194_H13wave_9b9a692a/run.json` 及原运行目录。
- `runs/checks/H13/20261006_061039/wave_audit.json` 及原审计目录。
- H12原报告引用的 `runs/baselines/H12_20261006/README.md` 原冻结目录也缺失。

任务书BF-00规定：“任何引用文件缺失/版本不明记BLOCKED，不得以模型推测代替源码”。因此不能把历史报告中的数字当成重新复算的原payload，也不能依靠当前n8的不同候选结果改判原45.2792倍能量缺陷。历史H12 PASS与H13 FAIL均保持；本轮结果不表示总体goal被标blocked。

## 当前可审查的冻结成果

输出位于 `runs/checks/BF-00/20261009_135037`。当前checkout没有Git元数据，版本采用实际文件SHA256树身份，不虚构commit。`source_inputs.zip`封存当前src/include/tests/cases/tools/CMake与方法/报告文档；`numerical_before.json`和后续 `numerical_diff.json`核验本轮零数值变化。

`historical_archive_audit.json`记录四份现存H13 evidence.zip的hash/成员目录；003903、035647、045122三份receipt及1,307个payload member验证通过。233555没有同名receipt，另以其原 `evidence_manifest.json` 验证123个成员，结果见 `oldest_archive_member_audit.json`；不虚称其存在receipt。当前可复核的20261009静滴72几何/1728解与n8完整256步原始输出仍在原归档；它们属于后续H13候选，必须与丢失的20261006原long-wave分开。新的BF-CSF与冻结native H12/H13路径做A/B；不把独立q/Hodge开发候选混入CSF color/material/力的职责。

当前1,397份数值源码/输入树SHA256为 `4b43dfb3f47ce53121bfdd64b18cec4b0bf640de37751365131822461e657483`；归档另含556份方法/报告文档。首次构建因原default H12契约exe被占用而LNK1168失败；改为本任务D:隔离testbin后Release全部15个相关可执行构建成功（262.60秒），未改源码、未终止其他进程、未降低编译/求解门槛。工具链为VS2022/MSBuild17.11.9、MSVC19.41.34120.0、CMake3.29.5-msvc4，x64 Release。

25项原CTest注册按原命令/环境/timeout映射到隔离可执行；只有临时工作目录和exe输出路径改变，原参数保留。首轮23/25通过（766.02秒）：架构检查的target graph输入被临时脚本误改为不存在的输出路径，修正本任务临时脚本后，同一实际架构检查通过（1.31秒）；该入口错误及首次输出保留。另一项完整momentum-geometry在180.03秒TIMEOUT，之后相同OMP2/180秒门槛、单独调度再次180.05秒TIMEOUT，不能判PASS或通过延长timeout消除。最终覆盖为24/25通过、1项TIMEOUT；包括所有7个原时空材料/XY-RZ AMR crossing检查通过。首轮/修正/单独复测分别归 `baseline_ctests.log/xml`、`architecture_corrected.log/xml`、`momentum_isolated_repeat.log`/`momentum_repeat.xml`；注册与程序hash见 `selected_ctests.json`，原错误入口版本保留为 `selected_ctests_initial.json`。

momentum超时前日志已完成原H09a/native质量/source检查和全部四邻域rank、8个modal source、4个含原始物理对流的modal evolution、2个真实移动metric断言；随后只完成部分末尾Hodge组合。它属于冻结的旧H13未采用开发候选，不是BF-CSF实现产生的新失败；仍完整保留该CTest性能未通过和未完成覆盖，不称全绿，也不修改数值文件/注册或豁免此失败。额外只读 `--kinetic-hodge-only` 诊断八组合全部通过（22.64秒），日志见 `kinetic_hodge_diagnostic.log`；它不能替代完整注册契约PASS。

## 符号与冻结门槛

实际 `CchPressureLoad::build` 构造 `-D_I.transpose()`；原H12工厂返回 `F_balance=C(J0*alpha)`，现有研究driver和契约以 `rhs=-F_balance` 进入 `A*u+C*p=rhs`。FirstJet凸液滴 `kappa>0`、`J0=-sigma*kappa`，故物理RHS为 `+sigma*kappa*C(alpha)`。BF-CSF将alpha替换为仅用于力的派生color时，必须直接返回 `+sigma*kappa*C(c_eps)`，不再额外取负。压力为 `p=sigma*kappa*c_eps+p0`，液内压力高；fixed rows由原wall reaction承担，不能要求外边界 `C1=0`。

原门槛从H12/H13 TODO与BF-V0–V12逐项冻结，完整表见 [gate_matrix.md](../../runs/checks/BF-00/20261009_135037/gate_matrix.md)。现存72个几何、864次actual analytic saddle的532,224行pressure和86,520个自由面独立复算：全部液内压力更高，压力/kappa误差最大1.91e-13，跳跃误差2.19e-14，自由面Cp-S尺度差最大1.09e-16，均保持原门槛。详见 [符号核验](../../runs/checks/BF-00/20261009_135037/rhs_sign_audit.json)；该数据检查不是BF新color/force实现。解析恒曲率机器平衡不证明variable-curvature能量/谱通过；BF-05全部PASS前不进入BF-06，BF-06没有全范围证据不允许H14。

## 供用户审查的恢复方案（尚未采用）

鉴于原payload已丢失，建议仅修订BF-00的归档前置：明确保留原H12/H13报告、原FAIL数字与缺失标记；采用本轮SHA256冻结的实际native源码/输入作为**新A/B基线**，后续用完全相同输入/周期/参数重新生成H12/H13与BF-CSF对照。新结果不能冒称恢复原run，不改变全部原频率>=1.8、能量/全谱、72/1728、密度/黏性、epsilon/时间、AMR/RZ门槛。旧切向模态扩展保留为未采用候选，BF-01–BF-05沿冻结geometry-only native source路径；不删代码、不改旧结果、不增加q依赖。

此方案只解除无法恢复旧payload的档案前置，不授权降低物理验收。需用户明确决定该BF-00方法/证据规则修订后，才能将其作为当前任务的通过依据并进入后续BF-ID。

本轮交付范围为五个现行方法/状态索引注记与BF-00 TODO/Report；源码、注册测试、输入及H12/H13既有报告逐字保持。`numerical_diff.json`确认added/removed/changed均为空；source_inputs.zip中1,397个数值成员已独立逐一hash复验。许可边界和完整源码版本见 `delivery_scope.json`、`source_receipt.json`。所有文件位于D:。
