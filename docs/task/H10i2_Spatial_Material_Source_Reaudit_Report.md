# H10i2 / 来源与空间材料设计独立复核

2026-10-04。结论：**设计/来源 PASS，engineering PASS**。原H10i设计FAIL保留历史；本独立任务证明当前几何来源已闭合，并固化可实施接口及预算，**没有实施新advector、causal phase/material F或最终native operator**。本轮0 production/0 tests/1既有research driver，0新target/tools。[任务书](H10i2_Spatial_Material_Source_Reaudit_TODO.md)、[实施前复核方法](H10i2_Spatial_Material_Source_Reaudit_Design.md)、[现行API](../api/H10i2_Spatial_Material_Owner_Contract.md)。下一唯一任务[H10j0实际advector准备](H10j0_Spatial_Advector_Reconstruction_TODO.md)，本轮不执行它。

## 1. 实际source链已经闭合

开始前独立重哈希H10i1f的210 frozen payload/87 canonical delivery，以及e268/d189/c175前置归档。逐项检查现行源码：H07 `traceHalfFace`来自原publication，geometry/full clock/selected cell检查保持；RegionSource由原owning view alias持有同一Snapshot，FirstJet调用complete native trace topology；MomentumMass检查原layout/D_I/region scalar和实际overlap native cell/bounds。wholeleaf纯证书与side-attached graph区别明确，未把侧边界重复计入体积。

已有16个原研究场景重新运行，全部完成publication/H07/真实dual mass。既有driver只追加三份被动坐标表，不修改原数值路径、Case输入或identity次序。独立后处理**不使用expected列作参考**，从actual端点/盒界与原平面位置/相号重算XY长度、RZ面积、液气体积及rho=1000/1的质量，并核验overlap加回dual。

| 重新核验 | 实测 |
| --- | ---: |
| 实际source场景 | 16 / 16成功 |
| actual spatial face / coincident face / zero half | 864 / 40 / 80 |
| true dual rectangle / pressure overlap | 864 / 1,472 |
| whole RegionSource / clip box | 16 / 32 |
| 独立物理/体积/质量及原普通矩最大误差 | 5.5228236329e−16 |
| 原CSV完整字节保持 | 39份（16 trace +22 ordinary source +1 affine） |

新 `source_chain_faces.csv/source_chain_duals.csv/source_chain_overlaps.csv` 保存实际bounds、原Snapshot一致、完整Initial clock、selected incident cell/kind/侧可用性和无库存。正长度zero面两侧分别核验；空间对称平均只解释现有几何诊断，不能据此赋予材料率。

## 2. 普通来源和能力边界保持

64个原移动圆形/球形布局重新运行；157,810查询（38,770 primal/119,040 dual）、5,460同物理重复行按source及0–4次闭式物理矩重新核验。原GL3/closure库存保持，最大独立误差5.53e−16。

fresh研究SVD仍有72 rank deficient和56 tangent未识别，后者全部局限于n4 selected2/3。后处理独立核查库存/范围和保存结果，**没有独立重算SVD**，不据此宣称新advector affine已通过。所有行保留；后续full16/n>=8 selected的可识别再现不能套用极粗信息不足的豁免。

H10h纯气体witness（XY 2:1、n16、rho=1e6、dof103）三层21个真实边段液相仍全零，三个revision独立，当前旧q仍12499.9875；独立读取其冻结质量/缺陷行核验三层液体体积、质量变化均0及旧R=q。这证明来源与旧材料转移错误的区别，不表示新compatible F/q/R已生成。

必要Release短契约 **7/7通过**：architecture/self-test/layout、phase geometry、momentum geometry、dual material、KEP。原source/time/domain/regrid/FP/lifetime/失败事务、两normal方向/XY–RZ/三clock、ordinary→side/multiple component/exterior/selected short结果保持，未认证RZ域角点/critical/allzero拒绝范围未改变。原82行flux数值由前置冻结日志/源码身份核查，本轮没有重跑flux；没有把缺省旧算法拒绝改写成新source失败。

## 3. 数学复核与明确限制

体积坐标上取 `A(Y)=(R(Y)-L(Y))/dx`、`C(X)=(T(X)-B(X))/dy`，则transfinite场 `div B0=A(Y)+C(X)`。面总量给出 `d=mean(A)+mean(C)`；原显式lift等价于

`H_X=integral(C-mean(C))dX`，`H_Y=integral(A-mean(A))dY`。

因此两normal边界lift为零，`div(B0-H)=d`。实际子面断点在两个方向分片积分即可；不能在一个Gauss区间中忽略对边trace断点。XY affine使A/C常数，lift恒零；RZ在s=r²/2、B=(r*u_r,u_z)上使用同式，axis H_X/B_r=O(s)、u_r=O(r)。原轴正规族与even-r² axial族可再现，一般physical-mean与point参考差异保留O(h²)。这不是全部smooth/PDE二阶证明，也不保证切向跨cell C1或已接FirstJet characteristic/Jacobian。

zero causal规则由右时间characteristic而来：`Q_l=integral(a^+ chi_minus+a^- chi_plus)dA`。实际a必须按同一prepared polynomial的变号区间处理；u_f或total Q为0不能替代a恒零。正反normal交换side/速度号得到率反号，独立子区间加回；缺外部流入侧明确拒绝。短range的offset/width、root opaque chart和independent事件保持，不从全局rounded点排序/重造root。两相分别积累，避免净量消除小相。

仅公式后处理作30组XY/RZ causal profile、normal反向/子区间加回/独立右特征triangle reference、60 missing-inflow拒绝规范，90个含明确非零R的shared-F/KE identity，以及75点zero-normal lift系数代入；最大误差3.9730263109e−16。另核查width=1e−25无法表示全局内部endpoint的shifted积分仍为正，RZ point/mean差异按四档加密恰降四倍。它们是数学规范自洽证据，**不是新生产factory或native PDE测试**。

`E_sp=q_d-T*q_p`及`G=Delta_m_d-T*Delta_m_p`直接给出`R_d=T*R_p+G/dt+E_sp=Delta_m_d/dt+q_d`；无需quarter或精确commuting。`N=diag(q_d)/2+K, K^T=-K`的材料KE identity保留真实R/uR与显式边界功。共享F抵消不证明local相准确或最终affine；可变材料所有cut行exact affine与KE的不相容条件不重新加回。

## 4. 原工厂签名的复用职责已补实

原H10j草案从 `phase/momentumGeometry/actualValues` 每次工厂调用同时准备几何和速度，却要求跨更新不重做SVD/atom图；没有公开prepared重构对象来明确保证复用。完整kinematics、atom来源与实际速度point/polynomial都塞入phase两文件也会混合独立职责。此处在实施前拆出必要单一接口，未隐含修改Core或扩大H10j预算。

- **H10j0**：新的CchSpatialAdvector.h/.cpp与现CMakeLists单一owner登记，3文件/1既有momentum test。原geometry preparation和actual velocity evaluate明确分开、共享immutablestorage，freshadvector identity与declaredgeometry/full clock保留；不读phi/密度、不生成phase F。
- **H10j**：现CchPhaseGeometry.h/.cpp，2 production/1既有phase test，`prepareSpatial(mesh,phase,actualSpatialAdvector)`只消费同一owner；root/phase积分与来源认证属于此阶段，不能重新拟合V。
- 后续H10k/H10l分别迁移源stamp/原入口、唯一材料F的dual视图与E_sp/R。当前调用库存只含现phase/material/dual/KEP四类短contract及既有research，无production App调用；每项开始前仍按实际调用复核4 production/2 test预算，不隐藏多文件迁移。

canonical registry必须加入kinematic piece断点；几何规模/速度更新局部线性，一次line排序可为O(N log N)，不能称全部准备严格线性。新cpp必须仅属于现CCH源表并从legacy glob排除，不新增target。typed来源可证明immutable输入与声明一致；CchNativeValues没有物理clock，新identity不冒称H14实际coupled commit。

## 5. 范围和交付

1,042源/346输入目录成员保持：仅既有research driver改变，输入仅driver/README改变；生产、测试、Core/FirstJet、F/R/uR/KEP/D_I/C、App/PDF不变。删除三段被动observer后原driver全文保持。独立math/source核验、必要短契约、actual最终driver/exe/DLL/退出和CSV、当前接口/固定预算归 [冻结目录](../../runs/baselines/H10i2_20261004/README.md)，[manifest](../../runs/baselines/H10i2_20261004/delivery_manifest.json)与[独立归档核验](../../runs/baselines/H10i2_20261004/verification.json)可重查。

本PASS关闭来源/实施设计复核，不关闭H10。仅唯一共享材料F、内部反号/实际AMR子面和、来源一致、KEP/skew作为材料硬约束；quarter/commuting/R=0仅诊断，Pure FirstJet完全独立、无库存/bubble/质量修复。最终native operator仍须joint constant/affine空间矩及material KE，内部反例先通过，再原H10b全部28窗口/84gate；原14失败保持，未重跑完整研究，H11/H08a等待。下一唯一任务H10j0，长期goal保持active。
