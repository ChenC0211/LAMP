# H08a2r15 / 闭域真实材料动能参考

2026-10-05，**PASS（消费者契约，production 0 / test 1）**。原四个 closed streamfunction 输入已用真实 q/R/uR 与 KEP 动能关系验收；旧“native 散度为零便要求所有 dual 零功”约束退出 gate。当前生产速度来源未切换。[任务书](H08a2r15_Closed_Material_Energy_Reference_TODO.md)、[实施前设计](H08a2r15_Closed_Material_Energy_Reference_Design.md)。下一唯一任务为 [r16 非线性运动参考](H08a2r16_Resolved_Nonlinear_Motion_Reference_TODO.md)。

## 参考与范围

只修改既有 `cch_kep_convection_contract.cpp`。新增独立参考直接读取同一 begin/end 几何相体积、密度和唯一共享 F，在真实 dual 四个侧面以几何符号累加 q，再构造 `R=Δm/Δt+q`。先差分相体积，避免相减两份大的混合质量；production 的 cells/report 仅被核对，不用作参考真值。

分别对两个速度分量验证：

```
Σ u C = 1/2 Σ q u² + P_boundary
Σ u(-C+uR) − 1/2 Σ (Δm/Δt)u² = 1/2 Σ R u² − P_boundary
```

旧 dense skew/load 参考继续运行。原四 closed 的 streamfunction、n4 XY/RZ × uniform/AMR、材料、probe、边界和 clock 均保持；三层几何质量逐项相等、材料质量差为零，原显式墙/轴的完整 normal 多项式及相/材料率严格零检查保留。旧 `CLOSED` 零功数值列保留为诊断，正式验收改为同一非零参考关系，没有清零或补偿 F/R。

本项没有引入 R=0、点态常散度、pressure→dual 交换、材料库存或 bubble；Pure FirstJet 不受反向修复。四硬约束仍为唯一共享 F、内部反号/真实 AMR 子面和、来源与完整 clock 一致、KEP/skew。

## 原失败的直接复验

先逐字复制冻结 r11 四文件候选，链接当前 r14 准备库，运行原测试，准确复现 **4 个 closed 失败**。然后仅替换测试消费者，候选全 KEP 契约通过，`MIGRATION_RESULT=0`。

| 原场景 | 实际对流功 | 真实 q 参考功 | 预算缺陷 |
| --- | ---: | ---: | ---: |
| XY uniform | −0.0003851941065879088 | −0.0003851941065879079 | 9.22e−19 |
| XY AMR | −0.0002357180460041972 | −0.0002357180460042321 | 3.49e−17 |
| RZ uniform | 0.004779169037928207 | 0.004779169037928206 | 1.74e−18 |
| RZ AMR | 0.005109499126268294 | 0.005109499126268258 | 3.56e−17 |

候选最大 dual 净出流 **0.0188984**，原真实动能身份仍成立。未把此非零项当成应被材料或几何修复的误差，也未宣称其精度收敛已经验收。

## 覆盖、负控和回归

新增独立材料参考在当前生产和隔离候选上各覆盖 **124 个原检查场景 / 6,144 个 DOF**，质量变化、q/R、uR 载荷、对流功和材料动能最大尺度化缺陷均不超过 **1.87e−16**，保持原 **1e−12** 门槛。按分量检查，不能靠两个分量的误差抵消通过。

全部原矩阵保持：16 组 algebra、12 组 bounded ports、24 组 crossing（密度 1/1000/1e6）、4 closed、原 44 个拒绝及 4 次 regrid；1,356 个 trace、1,516 个 native moment/3,504 个 pair 检查没有削减。paired-port 原独立矩阵和预算继续验收，新增物理参考不冒充 periodic 装配的新实现。

遗漏 uR、把 +uR 改成 −uR、误设 q=0 三种检查侧负控均实际调用同一误差门槛并被拒绝。两种来源上的显著缺陷分别约 **0.154/0.309/0.0895**；不修改接受结果、不添加生产来源分支、不放宽容差。

标准 Release 构建及 **9/9 CTest PASS**，197.69 s；其中 KEP 16.42 s，为本机测试耗时而非性能结论。剔除既有计时列和新增九行检查后，原 **10,773 行输出逐字保持**，包括 KEP 的原 484 行和四个旧来源零功诊断；候选新旧测试的原 SUMMARY 也逐字相同。所有生产、其他测试、driver/原输入保持，1,383 项源文件/输入仅预算内一个测试变化。

检查：`runs/checks/H08a2r15/20261005_200733_233`；冻结：`runs/baselines/H08a2r15_20261005`。`verification.json`核验来源哈希、原输入构造、完整矩阵、原行保持、负控及四闭域功；`candidate_negative/positive.log`保留前后结果与运行程序绑定。临时 CMake/source/runtime 仅归档，不是长期入口。上轮 r14 的 422 个冻结 payload 保持。

本项未部署正则来源，也未重跑 native/affine 和完整 28 窗口/84 gate，因为生产方法未改。r11 尚需 r16 给出可信跨分片 point/J 参考并裁决原精度，之后才可恢复部署与完整研究。H08a 同源时空积分、真实 FirstJet 桥接和 R 收敛未完成；H11 未启动。
