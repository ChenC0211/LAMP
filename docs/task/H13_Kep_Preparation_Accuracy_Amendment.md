# H13 / 原KEP准备残差精修的限定范围补充

状态：2026-10-09 / 081136，**已获联合方案授权，原两行补丁落地，标准重测待完成**。同一H13；H13整体仍FAIL，H14未启动。以下独立范围请求保留历史，现纳入[24文件联合预算](H13_Phase_Cut_Kinetic_Hodge_Method_Change.md)。

本补充请求只增加一个production文件 `src/Modules/Hydrodynamics/CchKepConvection.cpp`，将H13联合任务总预算由四份增为五份。两份已有regular/traction测试预算保持；原KEP契约已经覆盖失败场景，无需修改它。原[联合提案](H13_Native_Flux_Port_Method_Change.md)明确“不得修改characteristic积分/Jacobian公式、容差、材料F/R/uR、KEP或accepted FirstJet状态”，故本补充未获明确授权前不修改生产文件。

## 为什么需要这一处修正

原XY AMR非多项式平移圆、密度比1e6、常速度0.4的必要更新门槛1e-10失败，实际误差1.57403257095e-10。当前流量/权重/skew配对全部导出后，以80位Decimal精确复算原二进制系数：prepared系数引起的速度误差1.55043727005e-10，占实际约98.50%；apply求和仅约1.89394e-12。补偿运行时求和不足以关闭反例。

仅隔离恢复旧P2顶点消费、保持相同其余源码/依赖/原生输入/材料几何：旧语义误差6.57261e-11，通过旧门槛，但prepared残差仍贡献约6.83940e-11。该对照是明确标记的语义控制，不是旧字节档案恢复，不作为production fallback。新source的合法通量变化触发了原准备精度的薄弱处；不能靠改变新source来寻找偶然通过。

## 具体两行补丁

只在 `constrainOpaqueMoments` 内将两个内部残差精修停止系数从32改为8：

```cpp
// 原值
const double stopping=32*std::numeric_limits<double>::epsilon()*std::max(1.,maximum(target));
if(error<=32*std::numeric_limits<double>::epsilon())break;
// 候选
const double stopping=8*std::numeric_limits<double>::epsilon()*std::max(1.,maximum(target));
if(error<=8*std::numeric_limits<double>::epsilon())break;
```

8 eps来自同一文件已有bounded hard-row polish精度；不是按波频率、能量、网格、密度或case选参数。目标矩阵、rhs、局部支撑、SVD阈值、CG/最小范数方程、共享skew配对、最终256 eps验收、最大迭代数、材料F/R/uR、M/C及物理门槛全部保持。不改apply、已接受速度或FirstJet状态，不增加自由度或生产分支。

数学方程不变，仅更准确地满足已有constant/P1/P2约束；prepared系数会在舍入量级改变，不能宣称与旧系数逐位相同。求解可能增加迭代，须如实记录准备成本；不得降低精度补偿成本。

## 已完成的具体验证

证据目录 `runs/checks/H13/20261009_003903`。隔离代码 `CchKepConvection_refinement_candidate.cpp` 从当前文件复制，只作上述两处替换；production源SHA仍为 `5BEC243E003F8F060E1BC1EA8E72A94F29D0755002630D9DBACA2BF8E48CB30A`。

- 同一204行/三密度原场景全部通过；1e6峰值更新误差降为1.84324222552e-11，原1e-10门槛保持。Decimal复算prepared贡献降为1.87375e-11；相消舍入余项保留。
- **完整原KEP契约退出0**，没有修改原测试：124场景、24实际crossing、6144个DOF、1516原native矩、3504配对、57受控耗散场景、全部材料动能恒等式及三项负控；`MIGRATION_RESULT,0,PASS`。常速度最大误差1.84324e-11，原原子/矩/能量断言全部保留。
- 隔离构建无warning；原生产代码、source、载荷、永久测试及旧失败证据未改。该结果不是H13静滴/完整毛细周期/材料28窗口验收。

## 获授权后的实施与验收

仅把两行补丁落入上述一个现有文件，标准Release构建、architecture/repository_layout及全部受影响原KEP/dual/material/source契约重跑。通过后继续原H13物理验收，仍须建立新source对应的独立可信轨迹参考；本补充不授权修改第三份永久测试或任意替换参考值。

任何原门槛失败保留FAIL，不扩大文件或改变方法/容差。当前原四文件联合授权保持；此第五文件例外只在用户明确同意后实施。
