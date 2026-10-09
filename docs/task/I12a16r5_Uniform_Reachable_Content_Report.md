# I12a16r5 / 可达来源恒定含液量报告

**PASS（本项）；原四个smoke第一步全部PASS；原30/400步full及I12 / Gate I仍FAIL。** 见 [TODO](I12a16r5_Uniform_Reachable_Content_TODO.md) 与实施前的 [设计](I12_Uniform_Reachable_Content_Design.md)。独立库存守恒、几何精度预算及其分离保持，不恢复逐cut cell机器精度几何体积拟合。

## 修改与证明范围

4production+1test。Core新增原owned polynomial velocity在闭物理域、完整old→arrival时间上的有限速／浮点操作可达来源包围。复用r3速度符号及endpoint时间收紧，计入原128子步、midpoint、target、event插值的向外舍入界；保持原domain/CFL/FP约束，不接受opaque覆盖。RZ新能力限严格零径向速度。

Interface复用原native Q4闭leaf覆盖及同一库存wl/wg，纳入所有可能物理侧的同源显式history；联合含液量上下界**严格相等**时提取常因子L±=cT±，单独记录uniformReachableQueries。几何gas但库存非零时保留真实权重，包括1e-100尾量。不合并相近常数，不按小积分误差丢弃pending。

这是输运常量的代数求值分支，不声称证明任意原tracer的终止或唯一边界事件。原多壁／切向拒绝和first-jet B/J代码保持。原immediate boundary分支优先，非恒定来源继续既有完整回溯／仿射相积分。原速度evaluate、benchmark、库存装配/reflux、I10/I11及所有已有容限没有修改。

## 契约验证

24个XY/RZ轴向、uniform/AMR、gas/liquid积分场景通过，涵盖旋转、cubic速度、时间反向与角点，均与独立总通量原函数乘c比较。800条实际成功回溯（364 boundary、436 interior）均落在新来源包围中，boundary没有J，物理time在同一历史内。64个原map拒绝仍保留：22立即多壁、8event bracket多壁、34切向／非入流事件；新分支只给出唯一常量输运量，没有将其登记为成功map。

30个额外来源点覆盖零／非零壁坐标、nextafter、1e-100距离、1e-300及最小正double速度，含精确零坐标分量。13项拒绝检查覆盖非零gas库存、不同边界相、隐藏闭液泡、不完整history、foreign库存、1e-100库存尾量、opaque、RZ非零径向速度、CFL、时间范围、零duration、域外box及错误舍入模式。旧57主积分及其后r1–r4全部输出作为完整前缀保持。最终直接contract约0.898s，exit0。

保留开发失败证据：测试最初漏列已有bracket多壁拒绝；全域时间包围过宽后沿用原r3收紧（未增CFL）；测试库存容限1违反原严格小于1检查，改为.3；隐藏直线夹具违反原active-boundary声明，改为域内闭Q2液泡。没有为这些夹具放行production guard。

## 原驱动结果

smoke约11.990s，四个第一步全部通过。每个均处理全部native faces、半步和全步积分、库存装配、reflux及唯一几何发布；uniform每interval144faces，AMR248faces。8个装配/reflux报告的守恒闭合与reflux残差均为0；每个新state的M逐项等于原装配结果。原预算capacity*1e-6、4096panels、depth40不变。

full约14.769s，仍请求rotation30步、vortex400步。四个算例第一步成功，第二步的全部半步face积分成功，随后**原库存装配**拒绝“inventory violation exceeds roundoff budget”：uniform cell19、AMR cell18。rotation最后accepted时间.01，vortex.001；被拒绝步骤保留的accepted诊断与独立smoke第一步逐项相同。没有POLYNOMIAL_FACE_FAILURE，没有接受第二步或声称完整30/400步通过。

## T0与归档

首次两targets构建27.929s；全部39targets重建71.154s，均exit0，本项warning0。26CTest约170.70s，**26PASS**，原driver smoke由FAIL转PASS。其余24个旧测试完整输出逐字相同；原integral完整前缀逐字相同；旧crossing probe真实exit1及输出保持。架构12targets、383TU、908files、3118includes、5条exact legacy edges不变；81frozen、3PDF和r4 archive身份保持。

首次审计将closure总行数误写成8，遗漏4行initial closure；按原日志和驱动修正为12，完整复审PASS，保留首次失败。归档交付5源文件、TODO/设计/本文/入口、39runtime及原DLL，逐项验证archive/live SHA。

证据：[验收](../../runs/baselines/I12a16r5_20261002/verification.json)、[契约](../../runs/baselines/I12a16r5_20261002/integral_verified.log)、[驱动核验](../../runs/baselines/I12a16r5_20261002/driver_verification.json)、[smoke](../../runs/baselines/I12a16r5_20261002/driver_first.log)、[full](../../runs/baselines/I12a16r5_20261002/driver_full.log)、[回归比较](../../runs/baselines/I12a16r5_20261002/regression_comparison.json)、[范围](../../runs/baselines/I12a16r5_20261002/source_audit.json)、[审计](../../runs/baselines/I12a16r5_20261002/audit_verified_result.json)。

## 下一项

在原第二步记录旧M、各共享gross通量及认证误差、单元净更新与64epsilon预算，辨别真实离散输运／库存区域不相容与积分代表值选择造成的越界，再做与证据相符的共享通量有界性修复。不得clip M、提高舍入修复容限或恢复机器精度几何拟合。原full及时间/空间/AMR/曲率gate未通过前，I12继续FAIL，不进入Hydro。
