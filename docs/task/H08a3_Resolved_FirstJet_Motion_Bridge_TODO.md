# H08a3 / 实际分片运动到FirstJet的消费桥接

状态：**历史 FAIL 保持，不再作为下一独立任务**，见[历史报告](H08a3_Resolved_FirstJet_Motion_Bridge_Report.md)。2026-10-05 用户指定 r11 为基础设施最后关口；通过后必要的真实 FirstJet 运动消费接入并入 [H08a 时空相通量积分](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)。下文保留旧任务范围与证据要求供实现参考，不能据此再拆 source/Jacobian/低秩先决任务；只有可复现且直接阻断 H08a 的数学反例才允许例外。H11 未启动。

目标：让同一原始Pure FirstJet在指定步内物理时间消费H08a2r1实际运动和Jacobian，建立不可混淆来源的只读/候选桥接；只解决运动到FirstJet的消费事实，不同时实现相积分与材料提交。

允许修改：启动先审计 `HermiteGalsTrace`、现行FirstJet pullback和CCH motion公开边界；在实现前设计中固定≤5 production/build、≤2 test/build的具体名单。若桥接与全DOF候选发布不能落在此预算，按实际owner拆任务，先固定契约与来源。中立契约归Core或Interface现有消费边界，Hydro适配器可依赖它；禁止Interface反向依赖Hydro。

禁止修改：H08a2r1原运动与事件门槛、实际B/atom、压力/黏性、电场、材料F/R/dual/KEP、库存/bubble、逐cell体积拟合或质量修复。现有规定速度Core RK2仍是其原入口，不能把CCH分片速度伪装为符合连续迹前提的任意callback。

数学约束：phi与grad从同一原始FirstJet和同一离散映射取得，内部pullback为J的转置作用于原grad；实际source实例、两端点advector/geometry身份与完整clock必须一起保留。几何完全独立，材料仅唯一共享F、内部反号/AMR子面和、来源一致、KEP/skew四硬约束，真实R/uR不强制为零。

实现要求：先检查真实顶点/edge DOF在piece面、AMR knot和r=0上的消费语义，不能只测离面内部点；合法单侧选择须有公开依据，含糊/同时/外部边界历史缺失则明确拒绝，不伪造J。任意步内时间应来自同一source及原FirstJet，不独立推进三个snapshot后插值；保留copy/lifetime、mesh/regrid与失败不修改原状态的语义。来源结果和消费结果都不得冒充完整耦合接受事务。

验收：代码前固定T0/T1/T2、真实DOF覆盖和可执行误差门槛。T0为Release、目录/架构与原九项回归；T1含来源/time/instance/regrid拒绝、与公开同一trace的逐DOF消费对照、真实跳跃/单侧/轴/外部边界及失败原状态不变；T2为独立解析FirstJet输运的XY/RZ、uniform/真实2:1、静态/时变、多步内时刻的值/梯度精度与时间收敛。保留原反例，不把同一实现互相比对当解析验收。完整时空相积分与R收敛仍由H08a后续任务验收。

失败规则：混用来源/clock、绕过实际运动、任意选支或伪造内部J、任一预设精度门槛失败均FAIL；失败候选归档退役，无隐藏RK2/平滑速度/守恒库存fallback。只有真正外部依赖阻断可BLOCKED。

交付：预算与实施前设计、唯一typed消费契约/实现、独立逐DOF证据、API/报告/源码与运行冻结及未支持边界。任务通过不代表H08a材料时空积分、R收敛或完整CCH App已通过。
