# H07b / 同源正权分相求积与黏性验收

日期：2026-10-03。**H07b PASS；H07a/H07b 共同关闭 H07，下一项 H08。** 按[任务书](H07b_Phase_Quadrature_Viscosity_TODO.md)完成静态同源相求积及 H06 黏性接入；采用用户授权的[完全非守恒材料路线](H07_Pure_Gals_Hydro_Coupling_Design.md)。总目标与 Gate H 保持未完成。

## 实现与范围

[公开 API](../api/H07b_Cch_Phase_Quadrature_Viscosity.md)新增正式 Interface 只读求积 owner 和 typed CCH 黏性消费者。唯一 FirstJet Q4 的 MappedSegment 提供液气完整区间；GL4 内层、GL8 自适应外层产生实际相内正权节点，每点有严格相证书，RZ 使用同源 `2*pi*r`。原单元 P4 矩及离轴 reciprocal-r 矩的估计、独立 Green 相容积检查和有限工作预算共同约束求积。

CCH 复用 H06 同一 native strain/E 和 `K_mu=E^T H_mu E`，每个求积点按真实相取黏度；绑定 geometry、完整时间层、layout 和实际 mesh。零黏度精确消失；无材料库存、alpha 黏度混合、几何修正或权重缩放。选中子域、轴规则和 immutable 来源寿命保留。静态求积与几何体积的一致性不恢复材料守恒。

## 独立数值验收

平面/圆/球连续相区的参考由独立解析切片和高阶 GL24×GL12 积分计算，圆端点用角度变换处理平方根。参考不查询生产 roots/求积权重，也不进入生产接受条件。

| 检查 | 实测 |
| --- | --- |
| 相区场景 / 实际点相证书 | 22 / 75,632；XY/RZ uniform/2:1、纯液/纯气/切分/小相/选中半域 |
| 容积相对单元容量最大误差 | 2.2071e-13 |
| 独立缩放二阶相矩最大误差 | 8.3266e-10 |
| 多项式 native 场两相积分耗散最大误差 | 2.6311e-9；门槛5e-8 |
| `u^T K u` 与逐点积分最大相对缺陷 | 1.0952e-13 |
| 对称相对缺陷 | 7.7360e-17 |
| 4个约束轴后的矩阵谱 | 最低值−3.4428e-15，PSD至舍入误差 |
| 66组材料/缩放检查 | 单侧0、两侧相同、0/1e-6/1e-12同比缩放；缺陷1.0300e-18 |
| 4个实际相K进入H04参考求解 | 速度/压力最大绝对误差2.6042e-14 |
| 极小液、气相 | XY各2^-40；RZ各pi*2^-40，均有正权节点，满足自身相对误差1e-5 |
| 明确选中子域 | XY/RZ各56个CCH单元，覆盖112个publication单元的半域 |
| 明确拒绝 | 30；材料、geometry/time/layout/mesh、积分设置/工作限额、非FirstJet、溢出及非零轴值 |

拒绝不改变已有 prepared 对象。副本仍持有同一 snapshot；坐标改变和实际 regrid 使旧对象失效。两相黏度全零时矩阵nonzeros=0、载荷/耗散精确零，无黏度下也拒绝非零径向轴值。

两相 `mu_l=1, mu_g=.01` 的非多项式耗散 MMS 使用实际 native 面中心值和 `exp(q)-1` 隐式圆/球；参考是独立连续圆/球和解析应变积分。n=4/8/16，预定末档阶≥1.8保持：

| 序列 | n=4相对误差 | n=8 | n=16 | 末档阶 |
| --- | ---: | ---: | ---: | ---: |
| XY uniform | .569563 | .191057 | .0513169 | 1.89650 |
| XY 2:1 | .277934 | .0563475 | .0135797 | 2.05290 |
| RZ uniform | .298399 | .0898901 | .0253179 | 1.82801 |
| RZ 2:1 | .241434 | .0689926 | .0192590 | 1.84091 |

粗网格误差仍大；通过的是本任务规定的收敛门槛。没有把算子耗散误差或离散制造鞍点 RHS 当成连续两相 NS 的解精度。

## 工程验证与证据

5 production/build +2tests：四个新header/cpp、rootCMake和新短契约/testsCMake。前置1023个文件有1021保持，仅两个CMake修改，新增5个文件；312个case资产、H01–H07a所有实现/旧测试、Core与旧生产代码及冻结PDF保持。两个新cpp分别且仅属正式Interface/CCH；CCH依赖闭包不进入旧Hydro、Compatibility或App。

- local Release定向构建通过，无新增warning/error。标准构建重新链接 LAMP；exe身份从`D67FDA82...`变为`7F0C61DD...`，完整值在source audit。App源码保持，原五文件App契约仍通过12 smokes/49拒绝。
- 新契约 **1/1 PASS，30.23秒**；受影响回归 **13/13 PASS，76.24秒**：help、architecture/self_test、layout、FirstJet App、H01–H07a及旧pressure-gradient load。H01–H07a原数值汇总/收敛行与冻结前证据相同。
- 最终文档后再次通过repository_layout、本地链接和交付哈希验证。未跑全78项CTest、Debug或完整物理benchmark。

首轮误差估计漏掉曲线转折附近的小相区，独立体积检查明确FAIL；补充同一保留图的括区转折分段后通过，未修正权重或放宽门槛。原失败日志保留。随后修正测试的自由面边界payload、零参考压力的相对均值退化及拒绝消息预期，旧H04代码/判据保持；最终鞍点使用非零零均值参考压力。误差估计不是严格界，复杂未分辨几何或工作超限仍明确拒绝，未宣称任意cut均可用或已完成规模性能验证。

证据：[verification](../../runs/baselines/H07b_20261003/verification.json)、[operator detail](../../runs/baselines/H07b_20261003/operator_detail.log)、[regression detail](../../runs/baselines/H07b_20261003/regression_detail.log)、[source audit](../../runs/baselines/H07b_20261003/source_audit.json)、[交付清单](../../runs/baselines/H07b_20261003/delivery_manifest.json)。开始时核验H07a的23项交付/40项归档payload；TODO原身份、源码/输入清单、失败定位和最终审计脚本均归档。

## 未解决范围

H08 时间积分材料流量与非零质量连续性缺陷、H09 dual质量/动量一致性、H10/H14动量与动能材料项、完整边界/Hodge、CCH App和时间推进及低黏/EHD仍待后续任务。H07设计不设独立守恒材料库存；当前实现没有材料时间更新，也没有提前证明高密度比、无黏输运或电喷雾稳定性。
