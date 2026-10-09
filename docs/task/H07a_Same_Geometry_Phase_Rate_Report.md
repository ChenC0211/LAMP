# H07a / 同源 native 相测度与瞬时率验收

**状态：PASS，2026-10-03。** [TODO](H07a_Same_Geometry_Phase_Rate_TODO.md)、[API](../api/H07a_Cch_Phase_Geometry.md)。H07a 完成，下一单任务 H07b 真实分相体积求积/H06 黏性接入；H07 和 Gate H 尚未关闭。

用户本轮选择“完全接受非守恒材料输运”，已固定[设计](H07_Pure_Gals_Hydro_Coupling_Design.md)，并更新总纲、数值/架构/诊断约束与 H08/H09 路线图。材料直接来自当前 GALS 几何，无独立守恒材料库存或反向界面修复；后续明确保留材料连续性缺陷及其动量/动能一致性项。本轮只完成同源瞬时测度/率，未实现上述时间离散项。

## 实际改动

新增 `CchPhaseGeometry.h/.cpp`，归独立 `lamp_hydrodynamics_cch`，保留正式 FirstJet publication 与 H05 同一 metrics/layout/D_I。选中单元按同源体积测量得到液气容积与 alpha；真实 native fine subface 直接使用 publication 两条半面 trace。XY/RZ 按 offset/length 分别积分液气面积，保留短相，轴面面积与合法径向流率精确零。

新增不可变 `CchPhaseRate`，核验速度来源 geometry identity、完整时间层和实际 mesh。按真实面常值 normal velocity 得到瞬时液气 rate；逐单元外流消费原 D_I 的符号，内部共享与 2:1 子面求和没有第二散度。结果保留来源寿命，不持有主时钟或 accepted commit。

root CMake 显式注册新 source 与正式 Interface 依赖，测试 CMake 注册一个短契约。预算为 3 production/build + 2 tests，符合 TODO。

## 验证证据

证据目录：[H07a_20261003](../../runs/baselines/H07a_20261003/start.json)。开始时验证 H06 的 18 个交付与 37 个证据 payload；本轮全部 1,020 个原 source/build/test 文件中 1,018 个未改变，只修改 root/tests CMake，新增三个源码/测试文件。312 个 case 资产、H01–H06 实现/测试、旧 Hydro/Compatibility/App、LAMP.exe 和冻结 PDF 未改变。H01–H05 原数值 summary、H06 的 summary/MMS/order/spectrum 逐行相同。

真实 CMake graph 中新 cpp 唯一归 CCH；依赖闭包为 CCH、compile options、Core、通用 Hermite geometry、Eigen、正式 Interface、linalg、SuperLU。无旧 Hydro/Compatibility/Solvers/App 依赖。最终 Release build/rebuild 无新增 warning/error。

新契约 **1/1 PASS，4.00 s**；相关回归 **11/11 PASS，78.50 s**，另单独执行 architecture self test **1/1 PASS，0.49 s**，合计 12 项相关回归。覆盖 help、架构/目录、FirstJet App、H01–H06 和旧 Hydro gradient load；没有把 77 个默认 CTest 全部运行。

实际新契约覆盖 40 组、2,432 次单元测度检查、5,418 次面检查、224 次真实粗细侧聚合、20 项拒绝和 12 行收敛数据。包含 XY/RZ uniform/2:1、selected subset、纯液/气、圆/球、同一 native face 的两根、复制/寿命/调用者数组隔离、Initial/restart/Midpoint/End；极小液/气区间各 5 个，宽度 `2^-40≈9.095e-13` 仍为正。

拒绝检查含不同 geometry/layout、同一物理时刻但不同 clock fields、非法时间、缺选中单元覆盖、历史守恒 publication、XY/RZ 混用、坐标变更/regrid、非法数据、非零轴径向速度和流率溢出；失败不修正既有来源。

| 检查 | 最大误差 |
| --- | ---: |
| 解析单元容积，相对 V_c | 7.772e-16 |
| 解析面测度，相对 A_f | 9.437e-16 |
| 瞬时率解析/乘积检查 | 2.220e-16 |
| 面相率与逐单元/边界闭合 | 3.109e-15 |

非多项式 `exp((x-cx)^2+(y-cy)^2-radius^2)-1`（RZ 中 cx=0）的零集圆/球提供独立解析交集参考。对 n=4/8/16，全部四组相测量误差递减；误差为单元/面绝对误差总和除以对应总容积/总面积，未用解析参考修正几何。

| 坐标/网格 | 容积分割末档阶 | 面分割末档阶 |
| --- | ---: | ---: |
| XY uniform | 5.325 | 4.724 |
| XY 2:1 | 4.816 | 5.189 |
| RZ uniform | 5.414 | 4.754 |
| RZ 2:1 | 5.414 | 4.770 |

这些是该短算子序列的观察阶，不代表一般两相流或时间输运的阶。原始日志、独立数值复算见 [verification.json](../../runs/baselines/H07a_20261003/verification.json)、[source audit](../../runs/baselines/H07a_20261003/source_audit.json)；工程审计脚本归档于该 baseline 的 `engineering_scripts`。

初次实现误将 pure cell 的 raw Green ledger 当作精确容量，RZ 纯液因舍入差拒绝；已改为消费 publication 现有纯相认证容量，cut 仍用严格有界 raw 测量，没有新增 clipping/容差。构建中的测试调用与 target 名错误、缺覆盖 fixture 的 master 支持/显式接触点/图数前提错误均保留失败日志，未修改 Interface 或放宽 gate。

## 未完成范围

H07b 的液气正权体积求积与真正分相黏性未实现。当前 rate 不含 dt，不是时空相通量或几何质量守恒更新。H08/H09 的非零材料缺陷与动量一致性、H10/H14 的动能/时间离散、生产边界与 CCH App、完整 NS、低黏高密度比及电喷雾均待后续验收。总目标保持 active。
