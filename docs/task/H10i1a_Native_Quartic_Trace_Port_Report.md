# H10i1a / 同源 native Q4 相迹接口验收

2026-10-04；**PASS（只读 Core 来源能力）/ engineering PASS**。[TODO](H10i1a_Native_Quartic_Trace_Port_TODO.md)、[实施前设计](H10i1a_Native_Quartic_Trace_Port_Design.md)、[正式 API](../api/H10i1a_Native_Q4_Trace.md)。4 production/1既有短契约/1研究驱动，无新工具或 CMake target。上一轮为 progress；H10i1 的108 payload开工前复核。下一唯一任务[H10i1b闭合裁剪矩](H10i1b_Closed_Clip_Moments_TODO.md)。

## 实现事实

`LevelSetMappedSegment` owner 新增显式 `nativeQuarticTrace` factory，返回独立 move-only immutable owning query。普通 value Range 与零 Interval 严格分型，零区间没有 negative 标记；每侧按实际 native incident 叶的首个非恒零 normal Taylor restriction 认证相区间。首项阶数1–4及负侧奇偶号显式保留，缺外侧 unavailable。exact coefficients 全零才是零迹，未覆盖/二维恒零/预算不足等明确拒绝。原 constructor、legacy bool相/variation 仍拒绝零迹，无隐式half政策。

内部使用同一个 QuarticLeafData 的 restrictLocal，复用原 bounded root/square-free/dyadic/refinement；fixed native chart/order layout 使每个子查询/细化重新限制同阶偏导。ordinary 区间仍按原程序产生；没有零区间时提供原 legacy view。保留 rational charts 与独立 side events，未按 rounded global根合并；refined 是失败保持原对象的复制事务。normal始终正坐标轴，traversal反向不交换两侧。active支持不冒充selected支持；本接口无速度、密度、clock、材料 F、接受状态或全局DOF。

生产修改只在 `include/Core/Geometry/LevelSetMappedSegment.h`、对应 Detail、`Detail/LevelSetMappedCellLayout.h` 和 `src/Core/Geometry/LevelSetHermiteGeometry.cpp`。测试只扩展既有 `interface_root_enclosure_contract.cpp`，原 tiny-cap 验收正文/容限完整保留；研究只增加真实来源的被动查询和输出。

## 实际验收

短契约112个合法 XY/RZ uniform/2:1布局，504个零 chart；涵盖内部/原生面、正负相、非仿射正factor、normal阶数1–4、偶阶同相、反向与真实子段加回。RZ使用合法偶径向Q4和轴向Q4，保留r=0规则；轴上 `phi=r²*positive_factor` 的首normal阶数2、外侧不可用。另有16个真实coarse/fine chart，逐native盒证明2:1相邻、两侧相与覆盖；没有用“AMR=true”替代真实hanging证据。

12个两侧normal导数事件及12次同阶细化通过；包含跨零系数、exact multiple、端点身份、2^-48小相的相对精度、near-zero但非零、普通横穿/相切/多重根和混合ordinary/zero片。旧view遇零拒绝、错误轴/非轴对齐/null/wrong source、owner寿命、moved-from、FP环境、inactive/外侧与失败事务核验。新测度门槛1e−12实际尺度，原tiny-cap相对门槛不变。

真实Pure FirstJet已有出版八组 XY/RZ uniform/2:1正负平面：**1,408段新port均成功，1,368 ordinary段、40 zero段、40零chart、80侧相区间**。研究读取真实CchMomentumGeometry段和publication snapshot；usesSource/native active/order逐项核验。独立saved-data verifier从端点、平面符号与XY/RZ物理测度复算，relative physical error及side覆盖误差均0；ordinary ranges与原查询逐字段一致。这里只是相迹/面测量，不是closed volume。

零段在原生网格线上结束时，可以与仅接触输入端点的另一native闭叶相交。保留这类原point-only value incidence；实际研究逐一要求pointOnly、normalOrder0、retained input端点身份。它不是正长度零迹的有限根序列。首次saved-data verifier错误地要求零段events数全部0，修正为显式point-only核验，未改变数值门槛/相区间或root算法；修正前观察与脚本保留。

## 回归与来源

标准Release及16/16受影响CTest通过：architecture/self-test/layout、root/shiftedGreen、publication/phaseRegion、四项FirstJet、H07 phase、material/momentum/dual/KEP。追加显式hanging检查后重新构建、运行唯一修改的root契约通过；其余测试与Core生产源码未再改变。原32行CAP、KEP SUMMARY/NATIVE_MOMENTS与旧归档一致。

H10i的64普通source布局重新运行，**22份CSV逐字节一致**（16份数值表、6份仅有原表头），包含157,810个查询行、原GL3/共享与one-ring来源结果，以及原材料affine/compatibility/bounds/source输出；本轮未运行interval模式。研究原数值路径逆向移除新observer后与H10i1全文一致。原3,760行affine CSV逐字节一致，内部最大1.5849193387e−13，dof911反例7.4162898045e−14。源码catalog1,042项中仅4 production/1 test/1 driver变化；输入346项中仅driver/README变化。无新增生产/测试/CMake文件，App和冻结CCH PDF不变。

普通source/segments在增加point-only诊断列之前完成，使用其实际source/exe和`regular_source_run.json`独立归档；最终`trace-contract`与`audit`使用最终重建exe。两个研究版本移除observer后都与前版同一数值路径一致，生产Core期间无变化。所有实际退出码0、source/exe hashes与各自输出对应，不将后改列的source假称为早期64布局实际程序。

开发中RZ夹具先违反轴规则/径向域要求，inactive夹具带入了外部midpoint；按既有数据契约修正，未放宽生产拒绝。一个near-endpoint squared解析profile经浮点节点采样后不再是exact squared Q4，产生真实near-critical retained scalar；新旧API保持相同bounded拒绝，测试同时保留此限制，exact multiple用可精确表示的Q4另验。失败构建/测试与verifier预检均保留，非成功证据。

## 未完成范围

原16场景的`trace_contract.csv`/`trace_clip.csv`与H10i1逐字节一致：八个dual位置质量仍拒绝、16个clip盒仍拒绝；八个native位置publication仍拒绝，未运行的下游不算成功。只有新Coreport通过，未修复实际闭合质量、native拓扑/FirstJet派生volume账本或材料连续性。legacy敏感度不接受重合迹，高阶normal测量也不使临界界面可发布。

没有新材料F、库存/bubble、体积拟合或质量修复；仍只有唯一共享F、内部面反号/真实AMR子面和、来源一致、KEP/skew四硬约束。真实R/uR与原PDE/收敛门槛保留；clip内部迹、无速度诊断迹与材料因果迹分别消费，不加回quarter/精确commuting/R=0。

本轮未重跑原完整H10b或audit32；原14项完整研究失败保持，H10仍未关闭。H10i1b–f继续迁移真实来源消费链，再复核H10i、实施H10j及H08/H09唯一材料F/R，最后收口最终native affine/原完整H10b。H11/H08a未启动；port/source PASS不等于sharp-interface时间阶、无黏稳定或完整NS/EHD通过。

冻结证据见[baseline](../../runs/baselines/H10i1a_20261004/README.md)、[port核验](../../runs/baselines/H10i1a_20261004/checks/trace_port_verification.json)、[scope](../../runs/baselines/H10i1a_20261004/scope_verification.json)、[平面manifest](../../runs/baselines/H10i1a_20261004/delivery_manifest.json)、[独立归档核验](../../runs/baselines/H10i1a_20261004/verification.json)。goal保持active，无外部BLOCKED。
