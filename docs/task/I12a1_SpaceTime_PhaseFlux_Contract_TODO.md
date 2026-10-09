# I12a1 / Space–Time Phase Flux 区间契约与来源身份

## 任务 ID / 标题

I12a1：时间区间积分相体积的契约与来源身份。直接任务依据为用户编写并明确要求执行的 [任务书](I12a1_SpaceTime_PhaseFlux_Contract_TASK.md)。用户本次授权取代此前对 I12a1 方向的待确认状态；不要求再次确认，不推定获准实施 I12a2 或库存推进。

## 目标

只建立 integrated face volume 的 authority、真实时间区间、可追溯来源、派生平均 rate 和受检区间组合。先审查任务书要求的现有稳定身份能否复用；若触发任务书 §14，停止实现并给出指定格式的 BLOCKED 报告。

**执行中的用户补充（优先于任务书细节）：** “主要按照其中心思想，细节方面和现有程序有冲突的话以现有程序路线为主”。因此旧 velocity/tracer 未暴露数值 ID 不再机械触发 BLOCKED。采用程序已有 `shared_ptr<const Data/Revision>` 的真实不可变 owner 身份；复制共享来源，独立构造的不同来源不相等，不使用计数器、地址哈希或几何 ID 冒充速度身份。保留实际 geometry identity 和 `NativeMeshIdentity` 校验。

## 允许修改

- 按现有 V2/legacy 分离路线，新增 `include/Modules/Interface/HermiteGalsTemporalPhaseFlux.h`、`src/Modules/Interface/HermiteGalsTemporalPhaseFlux.cpp`，最多 2 个 production 文件。任务书建议的 `InterfaceTemporalPhaseFlux.*` 已承载旧数值积分，不在其中扩展 V2。
- 测试最多 2 个相关文件（含 CMake 注册）；优先新小型 interval contract，避免改变旧 time-flux 测试职责。
- 本 TODO、报告、文档入口和 `runs/baselines/I12a1_20261001/` 的只读审查证据。
- 使用持有实际 endpoint velocity views 的只读来源对象，以及持有实际 Domain/mesh 身份和固定现行 RK2 方法语义的 characteristic 来源对象；不可为空或仅为任意标签。区间只保留旧 Hermite handle，不创建第二份 snapshot。

## 禁止修改

不实现事件检测、时间求积、空间时间分割、crossing 求根；不修改库存、bubble、AMR/reflux、Core characteristic 数值公式、first-jet、Hydro、Electricity、字典/输出、容差/CFL。不得删除或隐式回退旧 midpoint operator。不得用计数器、地址哈希或几何 ID 冒充 velocity/tracer source identity。

## 数学约束

四个非负 integrated total/liquid positive/negative volumes 为基础量；signed volume 作差，average rate 仅由 volume/dt 派生。时间有限且 end>begin。同一旧 accepted authority、实际 mesh revision、velocity/characteristic source 必须保持；组合只在面方向/来源一致且区间首尾相接时成立。half/full 保持不同真实区间；失败不改变 accepted state。

## 实现要求

复用现有 HermiteGalsIdentity、NativeMeshIdentity 及模块 owner；先完整查明 velocity 与 tracer 是否已有可区分来源的稳定身份。若满足前提，按任务书 T1.1–T1.7 建立最小契约；若不满足，则记录具体声明/实现、可区分性反例、最小所需接口修改及超范围原因，不扩大 2 文件预算。

## 验收

- T0：身份前提审查；源码/配置/冻结文件身份核验。实现时要求目标编译和 architecture；未实现不得声称新目标构建成功。
- T1：按任务书合法/非法 interval、volume/rate、非负 gross、区间组合/错误组合、half/full 全部验收；缺少稳定来源身份时这些项目保持未实现，不把只读审查等同于契约通过。
- 回归：在未触碰源码的 BLOCKED 分支，核验原 executable 与源码身份，并运行既有 Hermite trace contract 验证现状；不声称该测试覆盖新的时间区间契约。
- 不修改 I12 既有失败期望，不要求 crossing 或 Gate I 本任务通过。

## 失败规则

严格执行用户任务书 §14：缺少可复用 geometry/mesh/tracer/velocity identity、必须更改 Core 算法、旧 phase-flux/reflux/库存 consumer 或必须作新方法决定时停止并报告 BLOCKED。报告逐项给出 `BLOCKED AT / WHY / MINIMAL REQUIRED INTERFACE CHANGE / FILES THAT WOULD NEED TO CHANGE / WHY THIS EXCEEDS I12a1`。

以上身份缺口按用户补充允许适配既有共享只读 owner 模式；若这种适配仍不能保持中心语义才 BLOCKED。其它不改数值公式、不接库存/事件积分、每次只做本任务的边界保持。

## 交付

本 TODO、仅针对 I12a1 的报告、实际检查命令/原始证据和变更范围。无论 PASS/FAIL/BLOCKED，任务结束停止，不启动或提出 I12a2。

## 执行结果

**PASS。** 按用户补充适配现有只读 owner 身份，2 个新增 production 文件 + 2 个测试相关文件完成契约。20 个相关可执行目标构建通过，18 项 CTest 通过；15 项原始输出完全一致，architecture 新文件计数单独核验。原 crossing 仍两次拒绝。没有实现事件/积分/库存接线。详见 [报告](I12a1_SpaceTime_PhaseFlux_Contract_Report.md)，本任务结束，不启动 I12a2。
