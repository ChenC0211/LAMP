# I12a15 / 有误差上界的面—时间积分报告

状态：**FAIL**。实现完成并通过 56 个积分场景，但预设均匀网格剪切场景在 4096 panels 内未达到误差预算。新 contract 保持 exit 1，未改成预期失败测试；I12 / Gate I 仍 **FAIL**，未接 Hydro。依据 [TODO](I12a15_Bounded_SpaceTime_Face_Integral_TODO.md)、[设计](I12_Bounded_SpaceTime_Integral_Design.md) 与 [API](../api/I12a15_Bounded_SpaceTime_Face_Integral_API.md) 归档本次实现、成功范围及剩余反例。

## 实现

对实际 native face × 时间区间，分别包围 total+/total-/liquid+/liquid- 四个 gross 体积，不把反号速度先平均。总量来自同一 owned polynomial 数据的解析 monomial moments；速度符号不确定时用点值范围与净积分关系包围 gross。液量使用同一旧 Q4、原离散 RK2 表达、完整 closed native leaf 覆盖和严格相符号，显式库存重载保留同一 r18 wl/wg owner。XY 使用 ds dt，RZ 使用真实 2*pi*r ds dt。

Core 在现有 strict/no-LTO TU 完成 outward 四项区间运算、合并/交集、误差宽度和有序代表值选择。Interface 负责来源与实际面、二分树及工作队列。返回私有只读结果和可由原装配器消费的单份 transfer；相未确定的最终叶允许计入误差界，但 branch 未认证的非零通量叶绝不允许返回成功。零通量入口仍验证域/opaque coverage 和轴正则条件。

为消除时间混合重复依赖，增加 `encloseInteriorRk2TightTime`；复用同一原 RK2 区间表达，将速度/导数范围与 endpoint hull 加明确浮点混合误差带相交。原 a13 查询保持数值不变。原速度数据、实际求值、RK2 运算顺序、CFL、轴、边界和缩步规则没有改变。typed unresolved 异常只允许细分尚未认证的有效 box；非法数据/非有限算术/浮点环境仍直接拒绝。

自适应方向使用候选区间误差与旧标量变化的工作启发；标量采样不进入积分值/上下界，不用来证明事件不存在。最终接受依据始终是四项区间宽度的 outward 和 E<=预设绝对体积预算，并且 pending branch 为零。该预算不能扩展 M 的 64eps 舍入修复额度。

修改范围为 4 production 文件（Core slab h/cpp、Interface temporal flux h/cpp）+2 test-related 文件（新 contract、tests/CMakeLists.txt）。原 M 装配、库存权重、几何闭合/发布、reflux、benchmark 参数和 gate 均保持。取消逐 cut cell 机器精度几何体积硬匹配的授权路线不变，几何与独立 M 仍分别按各自约束验收。

## T1：成功范围与剩余反例

57 个主要积分场景中 **56 成功、1 失败**。成功结果逐项核对解析值或独立参考及其明确误差带、四项子集关系、代表值与 transfer 一致性。时间半区间拼接、较紧容差下细化、重复调用的 bounds/work 一致性另作辅助检查。

| 项目 | 结果 |
| --- | --- |
| Core 面积分矩 | 8 组非单位坐标基、三次切向多项式、正/负速度、XY/RZ 横/竖面；另查反号 gross、区间合并与拒绝 |
| tight-time 原实际 B/J | 864 次正时长原 tracer 调用通过，含平移、时间反号、旋转、非线性和 subnormal 速度 |
| XY/RZ × uniform/AMR | 四组正反平移、两个 crossing 时刻、速度反号、薄层、纯相、零通量全部通过 |
| 反号速度 | 零净总量仍保留两向 gross 液量；参考每向 L 为对应 T 的 1/4 |
| 隐藏薄层 | 4 组通过；三个常用时间节点均为气相，节点之间的正液量没有漏掉 |
| 独立库存尾量 | 8 组通过；纯气相含液尾、纯液相缺额均保留；不改变 M |
| 非线性纯相 swirl | uniform/AMR 两组通过 |
| 旋转 cut circle | 57 panels，E=1.9668659661225273e-6 <= 2e-6；L 参考 1.0986319743114593e-3 |
| 非线性 cut swirl | 20 panels，E=1.5072711074137893e-7 <= 2e-7；L 参考约 1.9311904907136536e-5 |
| AMR shear | 3326 panels，E=2.4408823264994990e-6 <= 2.44140625e-6 |
| 来源/预算/边界/环境拒绝 | 24 项通过；包括 huge 容差下 pending branch 仍拒绝、opaque、foreign inventory、moved owner、stale mesh、轴及 rounding/FTZ |

统计中的成功主场景共 4508 个最终叶、17864 次分块评价；不含辅助重复/半区间/粗容差调用及失败场景。失败用例的最大 panels 与既定配置一致；本报告没有把达到预算的低阶近似解释成原 I12 时间阶或曲率收敛。

**未通过的固定反例：** uniform 4×4 网格，face x=0.5、y∈[0.25,0.5]，旧平面 x=0.484375，u_x=0.5*y、u_y=0，t∈[0.25,0.375]。原离散 RK2 在该场景有解析位移；精确 T+=0.005859375，L+=0.001953125，两个负向 gross 为 0。要求 E<=T*1e-3=5.859375e-6，深度上限 40、叶上限 4096。

最终在 4096 panels 时 E=**6.676659722917087e-6**，pending=0，程序明确抛出预算耗尽，未返回 transfer，也未改 M 或 geometry。因此该失败是现有区间积分在规定工作量下精度不足，不能宣称已通过，更不是需要用户补充输入的 BLOCKED。

旋转参考按原离散 RK2 半径因子 `(1+omega^4*d^4/4)` 计算，并以有限级数及余项误差比较。非线性 swirl 使用原离散代数的半径关系和单调根括区/时间上下 Riemann 和，另留明确浮点比较带。没有以连续 ODE 圆形不变替代离散 B_h，也没有在步内重置解析界面。

## 失败过程与 T0

首轮 production 编译完成，新测试因两个 lambda 未捕获局部 radius 而编译失败；修正显式值捕获。首次数值运行在上述 uniform shear 上耗尽预算，E=1.1030169844208343e-5。后续工作方向启发把薄层所需分块由约 3250 降至约 60，但没有独自解决剪切。tight-time 包围将剪切 E 降至最终值，仍超要求。一层候选前瞻曾得到 6.6754856789439359e-6，改善不足且增加成本，最终未保留该额外复杂度。所有初始构建与六次提前失败的原始日志均归档，未提高精度容差、panel/depth 上限或变更失败夹具。

最终 contract 改为汇总全部主要场景后返回总体失败，以完整记录已通过的旋转、涡旋、AMR 与 guard；没有将 uniform shear 重新归类成“预期成功的拒绝”。最终直接运行 3.7062877 s、exit 1；CTest 中该测试仍正常注册并报告 Failed。

37 个 runtime executable 目标完整重建成功，63.5772388 s；全部本项构建日志无 compiler warning。Core 算术与新 test 的 strict/no-LTO 选项从生成 vcxproj 核验。23 项相关 CTest 中 **22 通过、1 失败**，160.1565024 s；唯一失败为上述新增测试。21 项旧数值输出逐字相同，架构输出也相同：12 targets、383 production translation units、907 files、3114 includes、5 exact legacy edges。原 trace、a12/a13/a14、M/reflux、既有平移/跨面输出未变。

81 个冻结文件、3 份方法 PDF 和 a14 archive 身份保持；已知旧 crossing probe 仍为 exit 1、输出相同。归档包括 37 个 executable 与原运行 DLL，所有 delivery/evidence/runtime payload 逐条核验 SHA256；archive 完整性 PASS 不代表本任务数值验收 PASS。

证据：[I12a15_20261002](../../runs/baselines/I12a15_20261002/)；[最终验收](../../runs/baselines/I12a15_20261002/verification.json)、[算子与反例](../../runs/baselines/I12a15_20261002/operator_verification.json)、[旧回归比较](../../runs/baselines/I12a15_20261002/regression_comparison.json)、[源码范围](../../runs/baselines/I12a15_20261002/source_audit.json)；归档身份记录为 `archive_verification.json`。目录保留 configure/build/contract/CTest 初始失败与复跑脚本。

## 下一项

先单独写 TODO，解决同一 uniform shear 的认证误差预算反例，保持原 4096 panels 与原容差；可从区间依赖与相边界切分的数学收紧入手，不再反复增加前瞻或调整预算掩盖失败。该反例通过前不宣告本项通过。其后才推进一般边界历史/库存可行性和完整 driver 的 rotation/vortex、时间阶、动态 AMR 及原 I12 全矩阵。低黏/无黏 Hydro 尚未进入验收。