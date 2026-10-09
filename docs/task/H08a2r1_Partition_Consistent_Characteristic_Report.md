# H08a2r1 / 分片运动与Jacobian联合时间精度报告

2026-10-05，**PASS（本项运动来源）**。[任务书](H08a2r1_Partition_Consistent_Characteristic_TODO.md)、[实施前设计及补充](H08a2r1_Partition_Consistent_Characteristic_Design_Change.md)、[API](../api/H08a2r1_Cch_Characteristic_Source.md)。原H08a2的7个失败族及完整24族均通过；没有改变原窗口、初值、步长或门槛。H08a时空材料积分仍未完成，H11未启动。

## 最终实现与失败尝试

本项在同一实际CCH重构上公开cell/knot只读分支，新增prepared空间/邻接索引和不可变端点运动来源。位置与Jacobian使用同一Ralston三阶阶段、三次连续扩展及事件时间导数；事件切分后仍回到原固定时间节点。原Core RK2、FirstJet输运及空间材料F/R、dual/KEP不变，没有旧积分器fallback。

最初三阶候选仍在每个几何面重启，原七族中5个时变RZ族仍失败，最差J末档阶−2.58。此失败版本保存在 `cubic_restart_candidate`，不能据误差变小改判。最终方案在每个实际共享piece面、两个实际端点上核验原volume-coordinate B及其一阶导数：三点重建沿面二次Bernstein系数，导数按实际piece尺度归一，整面差以原双方系数工作量的128倍机器epsilon为舍入预算。零场预算为零，无幅值下限或可调容差。

通过舍入C1核验的连接只更换原分支查询，不重启步内曲线；真正不光滑连接仍定位并微分事件。每次穿越都检查实际连接，不以绕路连通替代面核验；所有速度/J仍查询原piece，未平均、平滑或拟合替代场。这个浮点核验是积分解析尺度，**不是任意幅值/病态多项式的有向区间精确C1证明**。

步末交点允许原64epsilon根定位舍入带：外部命中保留实际边界点且不返回内部J，避免常速恰好到边界时因一ulp误差无限缩步；透明内部连接按实际终点所在闭包取分支，精确在面上时选择后继。真实掠切、同时事件、方向/后继不唯一仍拒绝。

## 独立时间与分区验证

原24族 = XY/RZ × uniform/真实2:1 × n4/8/16 × 静态/线性时变，T=0.25、20/40/80步、每窗32原初值：**72窗口、2,304逐点结果、2,010实际几何穿越，24/24 PASS**。从保存的point/J独立Python重算RMS、单调性和阶，与C++裁决一致。所有位置/J RMS逐档下降，最低末档阶分别 **2.99968867 / 2.99883296**，最细最大scaled RMS **8.596475e-10**（原要求阶≥1.8、误差≤1e−4）。

原七个失败族先单独验证，再执行完整矩阵；七个原witness的最大径向—轴向伪导数为 **8.23e−18**。原uniform n8时变点由2.0488e−6降至3.2174e−18；原最细AMR点由2.2079e−6降至−4.9108e−18。这来自消除人工分区重启，未强行置零Jacobian。原witness完整映射差分最大误差6.04e−10。

24族共7,472条连接全部通过舍入C1核验，最大差/原工作量1.69e−16，所有仿射轨迹真实split数为0。对相同场/步预算跨六种网格分区比较全部2,304点，最大point/J差 **2.45e−16**（门槛1e−12）。并不据此宣称任意非光滑场分区无关或精确仿射。

## 跳跃、来源与拒绝

两个原真实cell/AMR内部knot切向跳跃仍为0.00933620 / 0.01484309，split事件 **1 / 2**；独立分段积分及解析映射J最大误差 **3.87e−17 / 5.33e−17**（≤1e−11）。一般native场288条轨迹、240事件，完整离散映射中心差分最大scaled误差 **4.70e−10**（≤5e−7）。真实速度跳跃没有被当作透明连接。

精确有理数脚本核验四个三阶dense条件；无事件仿射单步与独立Ralston稳定多项式一致，原Core RK2只作两种方法的解析误差对照。新增全零、1e−20和单位幅值、XY/RZ及uniform/2:1的12组连续性检查，含只有一个端点有跳跃的反例；小幅真实跳跃仍被识别。

检查覆盖：两端点/full clock、reconstruction及boundary law实例错配、原geometry身份与实际regrid、copy/lifetime、FP模式、RZ轴、XY/RZ全域/半域/孔洞/不连通边界、零时长、显式单侧/含糊起点、初始外向命中、步末和中间固定时间节点精确到面、透明及真实跳跃的nextafter极小横穿、掠切/同时事件和工作预算拒绝。外边界结果没有内部Jacobian；没有伪造边界历史导数。

## 回归、范围及交付

标准Release及**9/9短CTest通过**；最终补充极小真实跳跃断言后，仅重建/重跑受影响的characteristic契约（3.82s），仍通过。原五个数值契约 **4,727行确定性输出**与H08a1相同（排除原有两个计时字段）。1376份启动源码/输入中1372份不变；改动固定为5 production/build、2 test/build，其中新增3文件。未重跑H10o完整28窗口，因为其原算子/输入与确定性回归保持。

生产owner为 `lamp_hydrodynamics_cch`；唯一新实现 `CchCharacteristic.h/.cpp`，原advector只加分支查询。代码不读字典、不写文件、不拥有主时钟或接受状态。新测试在 `tests/contracts/hydrodynamics`，临时脚本、失败尝试与结果在 `runs/checks`，未向tools/cases散布调试入口。

**能力限制：** 当前是固定native mesh上的只读运动来源，完整clock/geometry标签仍为调用者声明；不是完整耦合提交。只有来源明确且事件可解析的轨迹返回J，selected外边界、掠切和同时事件不提供伪内部导数。尚未将此运动交给FirstJet，特别是实际顶点/edge DOF位于分支面时的消费语义仍需下一任务验证；尚无同源步内相积分、材料R收敛、完整CCH App/时间推进或EHD验证。动态准备效率和任意sharp-cut PDE不获普遍保证。

用户四硬约束保持：唯一共享材料F、内部反号/AMR子面和、来源一致、KEP/skew；真实R/uR保持。没有新增quarter、精确pressure→dual交换、R=0或界面库存/bubble要求。下一单任务[H08a3](H08a3_Resolved_FirstJet_Motion_Bridge_TODO.md)定义并验证实际运动到FirstJet的消费桥接，之后才接时空相积分。

检查：`runs/checks/H08a2r1/20261005_105119_173`；冻结：`runs/baselines/H08a2r1_20261005`。`verification.json`含独立72窗口裁决、分区比较、原确定性行对照及runtime哈希；`compile_sources.json`/`accepted_source`绑定最终7文件；`cubic_restart_candidate`保留失败证据。最终characteristic exe SHA256：`A5604C416A8FE2E0F5783216C8D8A5943DD557BFD26C34DD6FE7B0555104038A`。
