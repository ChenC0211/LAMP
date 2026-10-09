# H10i1 / 同源零迹与裁剪契约验收

2026-10-04。**PASS：只关闭数学/API契约设计；engineering PASS**。[TODO](H10i1_Coincident_Phase_Trace_Clip_Contract_TODO.md)、[完整设计](H10i1_Coincident_Phase_Trace_Clip_Contract_Design.md)。本轮0 production/0 test/1既有research driver，未实现新trace port、clip、材料通量或接受状态。**实际重合面source能力仍FAIL**，H10i历史FAIL保持。下一唯一任务[H10i1a](H10i1a_Native_Quartic_Trace_Port_TODO.md)实现Core只读查询；H10j等待完整来源链，H10/H11/H08a顺序不变。

## 设计关闭了什么

同一Q4在正长度面上恒零时，必须分清两种迹。几何clip按区域内部一侧的相迹积分，free curve与fixed clip重合部分只计一次；材料面用实际normal速度的backward-characteristic流入迹：`Ql=integral(a+ chi_minus + a- chi_plus)dA`、`Qg=Qtotal-Ql`。速度变号按实际多项式分片，a=0两相率均0，normal反转交换两侧且率反号。该因果规则来自输运方程的右时间极限，设计第6节给出独立平移质量导数反例；对称half材料迹会产生有限瞬时缺口，因此不作为新材料规范。

无速度的对称面测度如有需要，只作显式标记的空间诊断；它来自对称normal mollification，不能乘u冒称材料因果率。旧evaluate遇重合kind明确拒绝，直到独立H10j新factory提供正式求值。snapshot active支持不等于Hydro selected支持，实际外部流入须有边界相数据，否则拒绝，不暗用外侧scalar或ghost。

用户材料四硬约束保持：唯一共享F_rho、内部面反号/真实AMR子面求和、来源一致、KEP/skew。任意按上述规范得到的共享F仍可定义真实q/R、保持skew材料能量身份和uR一致性；无R=0、quarter或pressure→dual精确commuting要求。Pure FirstJet primary/Q4/transport/clock保持独立，无库存/bubble/体积拟合/质量修复。

Core私有同源`QuarticLeafData::restrictLocal`和`QuarticCell::homogeneousRestriction`已支持0–4阶exact偏导，因此不需复制Q4或添加另一几何引擎。新typed NativeQuarticTrace的恒零判据为所有restriction exact coefficients为0；逐真实侧的首个非恒零normal Taylor系数，结合原有有界root算法认证相区间与奇偶符号。derivative refinement也限制同一阶偏导；小相/近零非零仍走原严格range/root路径。缺侧、二维开区域全零、未认证退化及超预算明确拒绝。原legacy bool相Range、根和敏感度不能被零迹悄悄改义。

闭合Green账本也单独明确：几何内部迹在重合面两侧可不同，`sum(fixed)=exterior_fixed+coincident_interface_moment`。新typed几何报告须显式记录重合界面矩，不能冒用原共享迹的零assemblyClosure。真实clip区域内部相迹与唯一材料面的反号抵消不冲突。原敏感度接口仍拒绝拓扑不光滑的零迹，不能伪造0 Jacobian。

## 实际source能力审计

H10i的132 payload/32 canonical在开工前独立核验。研究新增被动`trace-contract`模式，n4四种网格×两位置×正负phi号共16场景。它调用原Pure FirstJet出版、CchPhaseGeometry、CchMomentumMass、正式段查询和HermiteGalsPhaseRegionSource；闭式矩形参考只用于观察，不进入几何、phase取值、材料率或更新。

| 实际位置 | 例数 | 原publication | 原真实质量/clip |
| --- | ---: | --- | --- |
| dual内部面，uniform offset0.375、2:1 offset0.125 | 8 | 全成功，H07 phase也成功 | 8个mass均拒绝；1,408段中40个零迹拒绝；两侧16个clip盒均拒绝 |
| native面，uniform offset0.25、2:1 offset0.5 | 8 | 全在geometric inventory阶段拒绝 | 尚未进入H07/mass；没有把未运行项计为成功或重复40个零段 |

所有实际拒绝均包含`positive-length coincident zero trace is unsupported`，无隐藏其他失败。16盒独立参考capacity为XY 0.125、RZ pi/8，正负phi号和内部侧确定liquid为0或capacity；saved reference独立复算误差0，但**实际clip没有成功输出**，其异常记录不能当作新port验收。observer exit0只表示能力库存完整。

源码将缺口定位到：旧MappedSegment禁止零restriction；NativeInterfaceChain实际边库存也实例化同一严格query；MappedFacePartitions只缓存legacy半面；完整volume discovery和publication派生测量链均消费它。NativeGraphChain还要求普通graph非端点样本在cell内部。因此不能只让query跳过异常，也不能把原生面遗漏出新路线。

设计明确正规native-side graph的真实支撑/端点来源、唯一edge arc、完整zero-set normal partial认证与inside-cell phase证书；不以采样认定纯相或伪造无限零根ordinal。原非正规多重/交叉零的publication拒绝保留。附着边的新几何frame只在认证过的真实closed-cell chart上评估，不取消普通graph的interior/正规规则。

## 分阶段实施

[设计第5节](H10i1_Coincident_Phase_Trace_Clip_Contract_Design.md)给出实际文件与消费者清单。H10i1a（4 production/1现短测试）实现Core同源trace；H10i1b（5/2）实现closed clip moments与真实Region消费；H10i1c（5/2）补正规native-side component/frame；H10i1d（2/1）补只读半面trace cache；H10i1e（5/2）迁移FirstJet派生volume/publication账本；H10i1f（3/2）迁移Region/H07几何consumer并核验全部16场景。每项独立TODO、固定gate后实施，本轮只登记下一H10i1a，未启动任何实现。

随后重验H10i 64普通source布局并复核其空间设计，才实施H10j causal phase-rate及H10k/l材料迁移；最终native内部affine先通过，再原H10b 28窗口/84gate。仅改善材料输入不能代替homogeneous PDE空间失败的收口。H08a时空积分与H11仍等待最终空间通过。

## 工程核验与限制

标准Release与7/7受影响CTest通过：architecture、self-test、repository layout、原root enclosure、phase region、momentum geometry、KEP。源码1,042项中仅既有research driver变化，生产和测试均保持；输入346项中仅driver/README变化。独立去除新增被动helper/include/dispatch后与H10i driver全文一致。原3,760行affine CSV逐字节一致，内部最大1.5849193387e−13，dof911反例7.4162898045e−14；原KEP SUMMARY/NATIVE_MOMENTS、App/PDF保持。

未重跑64普通来源/48材料/full/audit32：本轮没改production或原数值路径，原数据保留且新观察输出与原affine保持。设计推导与真实能力审计不是新port/trace/flux实现测试，也不是全Q4可发布、sharp-interface时间阶、无黏稳定性或完整NS/EHD验收。原H10b的14个失败与H10i FAIL不改判。

初次运行复制误用bin根路径，未启动数值进程；其null exit记录与一次PowerShell参考脚本语法错误留在`launch_preflight_failure`，不计为运行成功。最终correct runtime/source metadata、两个真实模式exit0及独立verifier为本轮接受证据，容差与参考未修改。

证据见[冻结目录](../../runs/baselines/H10i1_20261004/README.md)、[能力库存核验](../../runs/baselines/H10i1_20261004/checks/trace_inventory_verification.json)、[scope](../../runs/baselines/H10i1_20261004/scope_verification.json)、[平面manifest](../../runs/baselines/H10i1_20261004/delivery_manifest.json)、[独立归档核验](../../runs/baselines/H10i1_20261004/verification.json)。设计PASS与当前source能力FAIL分别记录，goal保持active；无外部BLOCKED。
