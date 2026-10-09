# H10i / 同源空间材料通量设计验收

2026-10-04。结论：**设计FAIL（真实重合相迹/裁剪source契约不完整），engineering PASS**。[任务书](H10i_Spatial_Material_Flux_Transfer_Design_TODO.md)、[未实施设计](H10i_Spatial_Material_Flux_Transfer_Design.md)。本轮0 production/0 test/1 research；没有新phase-rate、材料F、convection系数或接受状态。下一唯一任务[H10i1](H10i1_Coincident_Phase_Trace_Clip_Contract_TODO.md)先明确同源只读clip/相迹契约；[H10j](H10j_Spatial_Phase_Rate_Owner_TODO.md)等待来源关闭与设计复核。H10仍未关闭，H11/H08a等待。

## 授权与空间设计

按用户指定，材料耦合仅四硬约束：唯一共享F_rho；内部面反号共享/真实AMR子面求和；来源一致；KEP/skew及实际R/uR。quarter分配、pressure→dual精确交换、映射质量与严格R=0降为比较/诊断。固定mesh上保留`E_sp=q_d-T*q_p`、`G=Delta_m_d-T*Delta_m_p`，实际缺陷为`R_d=T*R_p+G/dt+E_sp=Delta_m_d/dt+q_d`。不向Pure FirstJet施加库存、bubble、体积拟合或质量修复。

H10h实际纯气体dual（XY 2:1、n16、rho=1e6、dof103）三层液体体积/质量变化为零，但旧q=R=12499.9875。H10i据此先设计与密度无关的native compatible advector，再按真实相区间积分唯一材料owner在primal/dual面的视图。方案保留原native总量`A_f*u_f`与D_I，使用一般one-ring physical-mean functional、体积坐标transfinite延拓及显式零normal divergence lift；不强制材料散度按pressure体积映射。公式推导详见设计，尚未数值实施。

RZ使用`s=r²/2`及`(r*u_r,u_z)`，正面积面functional按真实物理平均；axis正则限制不通过除零实现。一般RZ point/mean差异不能同时强制exact；原benchmark点参考与门槛保持。rank deficient时统一最小范数问题，显式报告未识别切向信息。canonical atom共享与prepared advector identity拟作为正式接口；新advector stamp不冒充H14 coupled commit。

设计预先拆H10j空间phase-rate、H10k迁移原phase/H08入口与来源、H10l唯一dual材料F/R/E_sp三个实施范围；但现有来源port不完整，**三项均未启动**，预算须按H10i1正式接口复核。真实质量始终来自同一publication，不由材料通量反推。

## 64个普通布局的来源核验

现有研究driver只增加被动`segments/segments32/coincident`模式。前两项消费实际移动圆形/球形Pure FirstJet，n4/8/16/32、XY/RZ、uniform/中央真实2:1共16网格组合；每组原full、半域、中心孔与不连通selected布局，共64。使用publication的owning geometry view alias调用Core正式`LevelSetMappedSegment`，不复制snapshot或重建rounded Q4。

| 来源观察 | 实测 |
| --- | ---: |
| 原native查询行 | 38,770 |
| 实际dual segment查询行 | 119,040 |
| 总查询行 / 独立物理段 | 157,810 / 152,350 |
| 相同物理段的重复查询一致 | 5,460 |
| XY/RZ正物理相测度闭合最大值 | 5.5228236329e−16 |
| shifted polynomial与独立GL3核验最大值 | 5.7858152345e−16 |
| rank deficient / 切向functional未识别 | 72 / 56 |
| full16布局最大切向缺口 | 1.1240151438e−15 |

实际range保留`fromEnd/offset/length`，核验局部多项式0–4次矩；RZ乘真实半径测度，GL3参考在每个真实相区间内计算。独立脚本从保存数据再核验全段闭式矩、共享来源/值与汇总，容差1e−12实际面尺度。此PASS只覆盖测度和几何信息，**没有验证新advecting velocity、F或材料PDE**。

56个切向未识别仅在n4的中心孔/不连通两类selected、8组布局；无法从实际单列/细连接数据要求二维affine exact。设计显式标注精度限制，保留这些行的来源/有限值/共享核验；不将其删除或伪标exact。独立脚本核验记录及库存，不冒称独立重算了SVD或证明其全局精度。H10h pure-gas witness全部三层21条真实边段液相均0，两方向物理闭合均0；这支持该box边界确为纯气体，不证明新材料率已经实现。

## 合法平面publication的重合面失败

增加被动真实能力反例：n4，XY/RZ uniform/2:1四组，平面offset分别0.375/0.125，位于pressure cell内部dual面，未落在native面上。**原Pure FirstJet publication、CchPhaseGeometry与dual拓扑均成功**。随后直接调用原`CchMomentumMass::build(mesh,geometry,phase,1000,1)`，四组均拒绝；逐段调用正式source查询也找到同因失败。

| 布局 | publication | 真实dual质量 | 查询段数 | 正长度零迹拒绝 |
| --- | --- | --- | ---: | ---: |
| XY uniform | 成功 | 拒绝 | 128 | 4 |
| XY 2:1 | 成功 | 拒绝 | 224 | 6 |
| RZ uniform | 成功 | 拒绝 | 128 | 4 |
| RZ 2:1 | 成功 | 拒绝 | 224 | 6 |

704段中的20段为正长度`phi=0` trace；拒绝原因严格核验为`positive-length coincident zero trace is unsupported`。Core当前在`include/Core/Geometry/Detail/LevelSetMappedSegment.h`的`isolated.identicallyZero`处分明拒绝；不是迭代失败或外部依赖。其余段查询未出现其他异常。observer退出0只表示反例记录完整，独立`design_verdict.json`仍为**FAIL**。

这同时说明两种契约须分开：裁剪体积积分需要区域内部的一侧相迹及重合free/fixed边界归属，防止体积虚增/双计；共享材料面通量在正测度零迹上需要全体消费者共用的规范。严格`phi<0/phi>0`在此不能覆盖全测度，两侧各自使用内部相也不能提供一份共享通量。不能直接用half trace、微扰、epsilon采样或cell材料平均填补。现有diagnostic rounded local dofs也不是另造精确几何的授权。H10i1须审计并设计公开认证port，不改变primary几何或输运路线。

因此，常规64布局通过不足以让空间source设计PASS。该结论限定于当前接口完整性；不否定compatible lift，不证明四硬约束或任意二阶方法不可能。

## 工程与历史研究

标准Release构建通过，architecture、architecture_self_test、repository_layout、原cch_kep_convection_contract **4/4通过**。生产/测试未改；CTest在最后增加被动coincident helper前完成，最后标准构建及全部4个观察模式再运行均exit0。原contract SUMMARY和NATIVE_MOMENTS逐字一致。

1,042源文件只有既有研究driver变化；346输入中仅driver/README变化。独立删除被动helper/流/调用/dispatch后与H10h前版本全文一致。原affine/compatibility/bounds/sources两档共8份材料CSV逐字节一致。原n8/16 affine audit共3,760行，内部最大误差1.5849193387e−13，XY 2:1 n16 dof911误差7.4162898045e−14。App与冻结CCH PDF哈希保持；Interface/H08/H09/KEP生产与全部测试源码保持。

未重跑full/audit32、未构造新phase/kernel/operator、未运行H08a/H11。原完整H10b 28窗口、15,116步、84gate中的10空间+4donor失败保持。rho=1、常速度时新旧材料率可能相同，材料输入改进不能代替原homogeneous最终native空间精度收口。后续须先内部affine精确/二阶，再原完整H10b；原空间阶≥1.7、最细L2≤0.05、Linf下降及KE/时间/donor门槛不改。

## 可核查证据

- [冻结目录](../../runs/baselines/H10i_20261004/README.md)：before/pretrial修订、当前delivery、最终源与exe、原工程日志及全部观察记录。
- [来源独立核验](../../runs/baselines/H10i_20261004/checks/source_verification.json)与[设计裁决](../../runs/baselines/H10i_20261004/checks/design_verdict.json)：64常规布局PASS、4个真实mass拒绝/20零段、设计FAIL。
- [scope核验](../../runs/baselines/H10i_20261004/scope_verification.json)：protected、原数值路径与8份材料输出保持。
- [平面manifest](../../runs/baselines/H10i_20261004/delivery_manifest.json)和[独立归档核验](../../runs/baselines/H10i_20261004/verification.json)：逐payload、当前canonical与H10h predecessor核验；manifest不包含自身及verification。

下一唯一任务H10i1；本轮结束但长期goal保持active。该缺口是内部可继续推进的数学/API工作，非BLOCKED。
