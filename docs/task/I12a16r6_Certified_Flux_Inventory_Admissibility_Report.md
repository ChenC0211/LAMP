# I12a16r6 / 认证通量区间内的库存可行性报告

**FAIL（完整修复未完成）；新算子契约PASS，原四个第二步half通过。I12 / Gate I仍FAIL。** 见 [TODO](I12a16r6_Certified_Flux_Inventory_Admissibility_TODO.md) 与实施前 [设计](I12_Certified_Inventory_Circulation_Design.md)。独立M、几何精度预算、原64epsilon舍入限额保持，不恢复逐cut cell机器精度体积匹配，不clip M。

## 成因证据

先仅增加benchmark诊断，原full除新诊断外与r5输出逐字相同。第二步half共有41个旧库存为0或C的单元超出舍入预算；最大越界／舍入预算约6.0992e7。所有独立单元的认证净更新区间都与[0,C]相交，故原每面独立代表值不适合直接宣称有界，但局部相交也不能证明存在全局共享选择。

独立C# BigInteger审计将所有binary64化为2^-1074整数单位，精确重算每个cell的上下界和原代表值更新。首次41行、最终82行均向外包围且名义更新逐位一致，未使用C++ exact expansion或其求和实现作为参考。见 [诊断验证](../../runs/baselines/I12a16r6_20261002/diagnostic_verification.json)。

## 实现与范围

5production+2tests。Core新增中立 [BoundedCirculation](../../include/Core/Discrete/BoundedCirculation.h)，用现有exact expansion加减/符号运行确定性最短增广路，严格验证可行circulation；不消除小残量。独立TU使用strict FP及无LTO。返回区间内binary64流量及其舍入后节点残量包围，不把非零亚正规残量登记为0。

Interface新增只读 `HermiteGalsCertifiedInventoryResult` 和批量 `assembleHermiteGalsCertifiedInventory`。只接受完整有序同源的真实BoundedFaceIntegral及同一旧库存snapshot。内部liquid±作为原认证区间内的共享变量；total±及物理边界流量固定，stock弧表达独立M更新和[0,C]。原代表值已处于旧舍入域内时沿用原装配，结果保持。

需要选择时使用Core可行流，再交给未改动的原库存装配/reflux，最终重新验证每个liquid±仍在原认证上下界，total/boundary未变。原证书保留在结果中，代表值改变仍受原四width误差界约束。没有第二个几何或M owner，actual tracer、速度evaluate、原库存装配、I10/I11及原算例参数未改。

## 算子验证与原驱动

8个小网络覆盖固定流、跨多边传递、分流、非零下界、反向残量、亚正规量、精确可行但binary64转换留下eta残量、已可行不变。独立已知流与割反例核验，包括所有节点各自区间可行但全局割不可行。12个XY/RZ、uniform/AMR、M/C=0/.25/1真实积分批次与原装配逐项相同；14项网络/FP/source/顺序/缺失库存拒绝通过。原integral全部输出保持完整前缀。

原四个smoke通过，除新增selection诊断外与r5数值输出逐字相同。原30/400步full四个第二步half均成功选择、装配、reflux并闭合：依次改变18/38/18/36个共享面，最大相对改变量4.9027062e-7，仍在原认证区间内。AMR vortex的第二步full也通过，改变37面；accepted时间推进到.002。

完整修复仍失败，原因已分开：

| 原算例 | 最后accepted时间 | 当前拒绝 |
| --- | --- | --- |
| rotation uniform | .01 | 第二步full，原舍入修复将face51的L+从认证上界2.7423358365237588e-26增至2.9784304699670792e-26 |
| rotation AMR | .01 | 第二步full，face191的L-从上界7.1956383242208442e-25增至7.1974406814916764e-25 |
| vortex uniform | .001 | 第二步full，face101的L-从上界2.6139143517257206e-26增至2.6152316729886983e-26 |
| vortex AMR | .002 | 第三步half，边界face29仍有1个未认证来源panel；41panels达到depth40，E约2.575e-19，未因E小于预算而放行 |

前三例被拒绝后的accepted误差诊断与独立smoke第一步逐项相同。没有接受越过证书的通量，也没有将任何full标为PASS。见 [驱动核验](../../runs/baselines/I12a16r6_20261002/driver_verification.json) 和 [完整日志](../../runs/baselines/I12a16r6_20261002/driver_full.log)。

## 验证与归档

首次诊断构建4.755s、诊断full15.451s；首次新算子构建14.013s。最终修复构建6.587s、直接contract .977s、full21.828s；全部39targets重建64.895s，smoke23.040s；26CTest350.98s，**26PASS**。记录本次耗时，不据此声称性能改善。23项旧测试完整输出相同，原integral完整前缀、原smoke数值及旧crossing probe真实exit1/输出保持。

架构12targets、384TU、910files、3121includes、5条exact legacy edges，增量只来自中立Core算子与两处生产依赖。81frozen、3PDF、r5archive身份保持；新增warning0。首次完整审计PASS。

保留开发失败记录：亚正规残量测试初始符号期望写反；随后将branch决策与Report分离时局部变量声明位置错误导致一次编译失败。脚本当时继续运行了旧exe，相关两次运行已标为 [排除](../../runs/baselines/I12a16r6_20261002/stale_binary_runs.json)，不计当前验证；最终修复后先检查构建成功再运行，完整回归及所有最终核验来自当前源码。

证据：[验收状态](../../runs/baselines/I12a16r6_20261002/verification.json)、[算子契约](../../runs/baselines/I12a16r6_20261002/integral_verified.log)、[回归比较](../../runs/baselines/I12a16r6_20261002/regression_comparison.json)、[严格编译](../../runs/baselines/I12a16r6_20261002/compiler_verification.json)、[源码范围](../../runs/baselines/I12a16r6_20261002/source_audit.json)。交付7源码/构建文件、TODO/设计/本文/入口、39runtime和原DLL，归档逐项核对live/payload SHA。

## 下一项

让**原同一个库存舍入修复kernel**在消费认证积分时同时遵守每面的liquid上下界，仍受原64epsilon限制；无证书的旧API运算顺序和输出保持。需要修改原私有库存装配kernel，超出本项预先限定的5production范围，因此单列下一task。随后核实第三步face29的旧库存权重与边界history差异，再完善该处来源认证；当前仅确认来源分支未被认证，不能仅凭小E推断其具体原因。原full、时空/AMR/曲率gate仍需实际通过，不进入Hydro。
