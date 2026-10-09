# H10d / 最终算子与成对端口联合设计（实现前）

2026-10-04。**本设计试验失败、已退役，不是当前生产API。** 源码/可执行与实际失败记录保留在`runs/baselines/H10d_20261004/rejected_paired_candidate`及该任务checks归档，结论见[报告](H10d_Final_Operator_Consistency_Report.md)。下文保留实现前推导与可识别性修订，不能解读为已通过方案。

用户要求直接设计最终native算子，保持constant/affine矩与material-KE/skew。独立H10c逐行审计证实外部误差与AMR内部三次余项均需处理；不改Interface、H08/H09质量/F/R或冻结方法。

## 已识别的数据缺口

旧研究把两外部segment的trace取同一平均，但prepared operator不知道配对、坐标image或native重复点。normal velocity的相反边界native点在周期空间是同一物理位置，却拥有两个半dual和独立导数；分别散射端口不能保证二次矩。共同trace平均只保证配对功，不定义最终空间拓扑。H10c仅内部矩与未闭合端口共同求解，还留下加密时增长的内部三次余项。

新增输入为显式segment配对及平移image。它是几何/边界关联数据，不是问题名称、解析参考、另一材料通量或新primary。仅支持可识别的matching translated segments，所有非零外部F必须配对，RZ禁止径向平移；非法或不相容关系明确拒绝。完整wall/open/backflow及一般不匹配周期mesh仍归H15。

## 单一prepared算法

旧opaque trace接口保持原内部矩能力；新增成对端口重载提供明确的线性端口拓扑。二者使用同一原K、同一one-ring最小改动算法，不设置数值fallback或AMR分支。

成对外部段贡献完整skew项：一份outward F在i,j间给+F/2、−F/2，diagonal是全部q/2。配对前核验同分量、反向关联、真实测度/平移及唯一F一致。实际native pressure邻接在配对的normal-component faces处闭合，one-ring仍只来自该真实关联。坐标差由声明的period images计算，不从case名称猜测。

若两个native采样点在平移image中重合，以E记录相同速度约束，两个半dual质量相加；不删除原DOF/体积/诊断，不增加独立primary。最终约束对象为`E^T N E`及显式端口，所有active quotient row的常数/一次/中心二次矩同时施加。J在共享quotient pairs上反对称，J*1=0。结果系数仍引用原native representative，不能称另一组物理速度。

所有pair在闭合one-ring支持上共同最小改动，旧rate尺度/局部SVD/稀疏联合求解/256eps残差拒绝保持。只能在独立wrapped最终矩、系数余项和完整原研究均通过后称整体一致性通过；拓扑修正不是预先保证PASS。

## 真实质量、R/uR与KE

duplicate class c的速度为U，M=sum(m_mid,i)，q_i为原唯一H09外流。原全native N含同一trace、skew J及q/2。对同一class计算

`a_c = sum_i(-N_i + U*q_i)/M`。

返回每个原native行`N_final,i=U*q_i-m_mid,i*a_c`，`load_i=m_mid,i*a_c+U*delta_m_i/dt`。所以class总N/动量载荷与原式相同，`load=-N_final+uR`保持；常速度a=0时逐行load=U*delta_m_i/dt，不把真实质量变化按比例重新分配。直接形成load避免大q相消，不能以修改端点材料质量保常数。

允许状态必须满足声明的相同速度关系，apply明确校验；不同值不能静默平均。因U相同，class的u^T N、动量总和和原材料KE/skew预算保持。pair端口提供同一中心trace，配对边界功按原公式抵消。该推导不重定义材料R，也不冒称一般非守恒KE严格恒定或全离散H14已经实现。

## 检查顺序

先独立逐native矩阵/真实apply与原误差分解，再新的production成对端口API、独立quotient dense reference、所有wrapped矩与旧来源/常数/KE/非零R/uR/regrid。旧解析外边界affine反例使用原opaque端口，必须保持。原n=8/16/32四组完整28窗口/15,116步、全部84原门槛最后复验，所有原native质量/行进入误差，不删边界行。

研究driver只提供实际配对/平移并调用生产apply，不执行自己的alias载荷投影或更新修正，不注入解析状态。研究donor对照保持原算法/结果，原控制FAIL不改判。若联合矩不相容、系数余项不具备正确h尺度或完整空间失败，本任务FAIL；不提前进入H11/H08a。

## 联合可识别性检查后的修订（实现前）

实际XY 2:1 n=8的全部quotient中心P2硬约束不相容：1,440行、3,864共享pair、秩1,422，独立COD最小二乘最大残差约7.01e−4；仅常数/一次的720行秩714，最大残差为0.002441×256eps。均匀n=8/16已能闭合。失败不是迭代预算不足，不能以加大迭代或放松硬矩解决。

用户要求的constant/affine reproduction和skew继续严格共同施加。旧opaque接口仍保持H10c已通过的内部中心P2。新增成对接口改用硬C/P1与受控P2余项的单一约束最小改动：在率尺度pair未知量x上最小化`|x|²/2 + sum_i |B_i x-b_i|²/(2 lambda_i²)`，满足`A x=a`。B是同一最终中心二次矩，`lambda_i=(h_i/L)^2`，L为真实native域最大跨度（包括声明的周期长）。因惩罚松弛引入的二次残差为O(h^4)，其对导数的贡献应比O(h²)目标更小；原不相容残差的实际h尺度仍必须独立测量，不能预先声称二阶。

求解`(C C^T + D)y=b`、`x=C^T y`，D仅在二次行有lambda²；局部SVD白化时同步变换D，禁止污染常数/一次硬行。收敛按完整线性系统和原硬矩同时检查，报告原二次余项而不将其称精确再现。lambda公式在数值试验前固定，四组合/各加密统一，不按AMR或参考误差选择参数。若实际二次/三次余项尺度或原完整研究失败，保留FAIL证据并拆后续任务，不改变H09质量定义。
