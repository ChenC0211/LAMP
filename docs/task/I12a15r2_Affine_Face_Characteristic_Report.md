# I12a15r2 / 面—时间仿射 RK2 包围报告

状态：**PASS（本项）**。新增原浮点 RK2 的相关性包围；依据 [TODO](I12a15r2_Affine_Face_Characteristic_TODO.md)、[设计](I12_Affine_Face_Characteristic_Design.md) 与 [API](../api/I12a15r2_Affine_Face_Characteristic_API.md) 验收。**a15/r1 的 uniform shear 积分仍 FAIL，I12 / Gate I 仍 FAIL**；本项尚未将新包围接入相积分或完整 driver，未进入 Hydro。

## 数学与实现

原 face×time 用相同的 q_s/q_t∈[-1,1] 参数化。每个 departure 分量返回 `c+a_s*q_s+a_t*q_t+R`，R 包含严格 outward 非线性和舍入残差。一次项保留空间—时间相关性，不将所有点都替换为整个 departure box。

内部 affine arithmetic 用 interval 表示常数/两项系数/剩余项。乘法保留一次项，将变量部分乘积及常数×旧 residual 向外累计到新 residual；各原浮点操作另外加入 eps*maxabs(value range)+denorm_min 的 RN 误差带。除法、精确零/一、反号与系数代表值转换都按设计逐项处理。复用原 polynomial component/Horner、`rk2Midpoint`、`rk2Target`，没有另写回溯法。

时间混合的 `.5` 分支另包围 `(a+b)*.5` 的舍入误差；dt/2 舍入为零时，arrival-only 与 mean 两种目标运算的差异用 dt<=denorm_min 下的明确误差带覆盖。零端点的混合保留精确零，因此 RZ 轴分量的 c、两项系数及 residual 仍全部为零。入口先验证原 a15 tight-time 完整 branch witness，不改变 CFL、轴、physical boundary、opaque coverage、缩步或 wall hit 规则。

物理点查询按真实端点向外反算 q，并用区间端点直接除法避免 tiny 分母倒数溢出。已知正的 binary64 端点差至少为 denorm_min，这一事实用于分母下界；已验证原物理点属于域后，可将参数范围与[-1,1]相交。点查询最后与同一原 branch box 相交。它没有平移、裁剪或修复实际轨迹/库存。

这些包含性依据来自原算术表达的逐操作归纳及严格区间内核；下面的样本是实现验证，不是对所有输入的概率性替代。Core 返回不借用对象的值记录，不冒充具有 Interface geometry/velocity identity 的发布结果。后续 Interface 消费者仍须保持来源与实际 native face 的一致性。

修改范围为 `PolynomialVelocitySlab.h/.cpp` 两个 production 文件，新 `interface_affine_characteristic_contract.cpp` 与 `tests/CMakeLists.txt` 两个 test-related 文件。root CMake、原实际速度/RK2、a13/a15 已有查询、Q4、M、库存权重、几何发布/reflux 和 benchmark 参数均未改变。

## T1

最终新 contract：**63 个区域、9,129 次原实际正时长回溯、514 个旧时刻恒等映射检查、18 类拒绝检查通过**。每次正时长回溯都核验实际无 hit、单 substep、Jacobian 可用，B 落在相同物理点的 affine 查询范围内，包含性容差为 0。

| 组别 | 覆盖 |
| --- | --- |
| XY uniform/AMR | 两种实际面方向；zero、constant、rotation、reversing rotation、shear、cubic/reversing cubic、subnormal velocity |
| RZ uniform/AMR | 两种实际面方向；zero、axial、strain/reversing strain；另查真实 r=0 轴面 |
| 非 dyadic | 两个面方向，非单位坐标基与非整齐端点/time；对比原实际 tracer |
| 极小时间 | 到 3*denorm_min 的区间，逐个 0/1/2/3*denorm_min 查询，覆盖 half-time arrival-only 分支 |
| 独立解析参考 | 常速、shear、rotation 的原离散代数；dyadic 配置中运算精度足以表示这些参考，非连续 ODE 替代 |
| 守卫 | 非法/退化/反向 face、非法/非有限 time/point、无效 domain/default proof、opaque、inflow、axis、CFL、rounding、FTZ |

常速点查询宽度小于 1e-10，验证结果确实保留 q 的相关性。RZ 轴不仅最终 box 为零，原输出 affine form 的所有径向项也逐项为精确零。

固定 shear 的同中心 face/time 同时减半得到：

| 细化级别 | 最大 departure 残差 |
| --- | --- |
| 0 | 7.8125000000013167e-3 |
| 1 | 1.7089843750003101e-3 |
| 2 | 3.9672851562527203e-4 |
| 3 | 9.5367431640881712e-5 |
| 4 | 2.3365020752202973e-5 |

最小相邻缩减比 4.0816326530，符合局部二次剩余项的预期；还核验 residual 不小于该解析双线性映射必需的交叉项幅值。这是新包围随输入区域缩小的局部验证，**不是完整 Interface 时间收敛或 I12 原剪切积分已通过**。

## 初始失败与 T0

首轮 production 编译遇到 interval namespace 的 magnitude ADL 同名歧义，改用明确的 affineMagnitude 名称。首次数值运行发现非负舍入半径经 outward 算术后，其下端可在 underflow 处略小于零；对称带只需要已核验的非负上端，因此修正该多余断言，没有减小误差带。第二次发现零 endpoint 混合仍附加了非零的通用 half-time 舍入带；补上实际代数 `0` 混合恒为 `0` 的精确分支，恢复 raw affine 轴零。原日志和退出状态均保留。

第三次新 contract 通过（13 guards），随后增加非有限/无效 domain/default proof 等守卫到 18 项；最终直接运行 0.0758688 s、exit 0。

全部 38 个 runtime executable 目标重建成功，65.9122114 s，所有本项构建日志无 compiler warning。Core 与新测试从生成 vcxproj 核验 strict/no-LTO。24 项相关 CTest 用时 161.3568426 s：**23 PASS、1 FAIL**。唯一失败仍为旧 `interface_bounded_spacetime_flux_contract` 的同一 uniform shear；新增 affine contract 通过。

原 23 个测试的完整输出逐字相同：22 个数值输出（包含原已知失败的积分 contract）加 architecture。架构仍为 12 targets、383 production translation units、907 files、3114 includes、5 exact legacy edges。原 trace、a12/a13/a14、原 M/reflux/平移/跨面未变；uniform shear 仍是 E=6.6763159339198168e-6 > 5.859375e-6、4096 panels、pending=0。

81 个冻结文件、3 份方法 PDF 和 r1 archive 保持身份；旧 crossing probe 的 exit 1 与输出也未变。第一次审计因 PowerShell 将 negative 分词换行而未识别已保留的失败文字；仅调整该诊断匹配的空白处理，原数值输出比较仍逐字执行，初次审计结果保留。

归档包含 2 production+2 test、5 份文档，38 个 executable 与原 DLL；所有 payload 逐条核验 SHA256。证据：[I12a15r2_20261002](../../runs/baselines/I12a15r2_20261002/)、[验收](../../runs/baselines/I12a15r2_20261002/verification.json)、[算子](../../runs/baselines/I12a15r2_20261002/operator_verification.json)、[旧回归](../../runs/baselines/I12a15r2_20261002/regression_comparison.json)、[范围审计](../../runs/baselines/I12a15r2_20261002/source_audit.json)。

## 下一项

将同一 old native Q4 的权威多项式施加于相关 affine departure，保留 phi(B) 的一次参数项和严格残差；跨叶继续完整覆盖，库存仍使用原有 wl/wg。随后认证分块中明确液/气的子区域并计算其真实 XY/RZ gross 通量界，以解决原 fixed-budget shear，保持 4096 panels 和原误差预算。一般边界历史/库存可行性、完整 rotation/vortex、时间阶、动态 AMR 和低黏/无黏 Hydro 仍需后续独立验收。