# H08a2 / 实际分片速度的特征线来源与跨界导数

状态：**FAIL，2026-10-05**；见[报告](H08a2_Piecewise_Characteristic_Source_Report.md)。真实跳跃/离散差分通过，但24族时间精度有7族失败；候选已退役，生产恢复H08a1。下一[H08a2r1](H08a2r1_Partition_Consistent_Characteristic_TODO.md)，父H08a尚未通过，H11未启动。下面及[实施前设计](H08a2_Piecewise_Characteristic_Source_Design_Change.md)的原范围、门槛保持，未因失败改判。

目标：在同一 actual `CchSpatialAdvector` 分支上定义可追溯的步内运动与其Jacobian，显式处理cell和AMR内部knot的切向跳跃，使后续FirstJet步内几何与材料积分有可消费的同一运动来源。

允许修改：先形成实施前设计，固定唯一时空运动owner、文件清单及≤5 production/2 test预算。归属Hydro运动学，复用H08a1局部jet；若涉及Core通用事件或Interface的来源入口，须按真实职责和原API契约明确接线，不能在回调中藏第二个trace算法。短契约在tests，完整时间研究在cases/benchmarks；临时定位在本任务时间戳runs/checks。

禁止修改：现有B重构/共享atom的数值定义、F/R、dual/KEP、压力/黏性/电场、FirstJet primary/几何拟合。无库存/bubble、独立平滑速度、任意callback伪造来源、历史守恒driver或失败fallback。不提前宣布H08a或H11完成。

数学约束：仅唯一F、内部反号/AMR子面和、来源一致、KEP/skew四项材料硬约束。局部Hdiv场有真实切向跳跃；不得强加全局切向连续作为新材料约束。事件处理必须给出同一**离散**轨迹映射的导数，说明时间定位、单侧选择与误差预算；仅积分局部J不能自动覆盖跨跳跃事件。连续ODE的saltation公式也不能未经验证就等同于离散RK事件导数。

实现要求：保留native mesh/selected domain、实际分支/axis/boundary及完整时钟；明确冻结速度和时间变化源的支持范围，端点速度必须来自同一prepared重构和明确来源。跨cell、内部knot、grazing、多事件同点、selected-domain缺口及regrid必须有确定契约和拒绝条件，不靠位置扰动或任意最小时间步掩盖。不得把三个独立几何快照插值当作步内几何来源。

验收：实现前固定T0/T1/T2的可执行矩阵和阈值，先通过一个真实非零切向跳跃的横穿反例及离散映射独立差分，再测试XY/RZ uniform/2:1、轴、selected/disconnected、单/多事件、时空来源拒绝和真实regrid。常速/affine与原连续Core契约的相容范围须明确；来源/数值异常不得静默切换。原H10o/H08a1和材料/KE回归保持；父任务的高密度比、moving/crossing相积分与R时空收敛不在本项改判。

失败规则：无法在实际B上给出一致运动与导数、来源混淆、跨界反例不通过即FAIL，保留证据和下一个明确数学动作；仅真实外部依赖才BLOCKED。不把“识别到跳跃”诊断本身算作本任务PASS。

交付：实施前设计、实际来源/运动接口与实现、独立导数和事件验证、API及报告。完成后再接同源步内FirstJet几何和相积分，以及primal/dual共同材料积分与R/uR/KE；不扩写单个segment trace补丁。
