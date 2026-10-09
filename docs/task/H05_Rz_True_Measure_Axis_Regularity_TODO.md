# H05 / native CCH 真测度与轴正则

**状态：PASS，2026-10-03；验收见[报告](H05_Rz_True_Measure_Axis_Regularity_Report.md)。** 前置 H01–H04 PASS；下一项 H06。

依据：[路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) H05、冻结 CCH §3.2/§6.1/§12.1 与 H-R1/H-R2。目标只解决一个数学事实：由唯一积分散度和真实体积得到物理 cell-average divergence，RZ 轴线以明确正则条件处理，无经验性1/r补项。

## 目标

实现只读 prepared CchPhysicalDivergence：保留同一 D_I/layout，消费 Core 真 V_c/A_f，返回 (D_Iu)_c/V_c。径向轴面必须 u_r=0，零面积不能隐藏非法轴速度；提供径向速度奇、轴向速度及压力偶的定型轴反射规则，使轴向 Neumann 条件可由偶延拓实施。

## 允许修改

- `include/Modules/Hydrodynamics/CchPhysicalDivergence.h`、`src/Modules/Hydrodynamics/CchPhysicalDivergence.cpp`：真实测度端口、同源体积归一化散度、轴面验证与反射规则。
- root `CMakeLists.txt`：加入唯一 CCH source list。
- `tests/contracts/hydrodynamics/cch_physical_divergence_contract.cpp`、`tests/CMakeLists.txt`：独立解析测度/通量/轴正则/MMS/拒绝契约和注册。
- 本 TODO/API/报告、现行索引/路线图；证据在 `runs/baselines/H05_20261003`，临时脚本在 `runs/checks/H05/<时间戳>`。

预算：3 production/build + 2 tests。

## 禁止修改

H01–H04 实现及原测试、Core/Solvers、旧 Hydro/Compatibility、Interface/Electricity/App、case/字典和冻结方法。禁止重装 D_I/C、人工轴面积、径向epsilon、速度clipping、经验1/r修正或fallback。H06 的 strain/hoop 黏性、H15 的 AMR pressure Hodge/完整边界、流体时间推进和耦合不在本项实现。

## 数学约束

- V_c=π(r_+²-r_-²)Δz，径向 A_f=2πr_fΔz，轴向 A_f=π(r_+²-r_-²)；XY 保留单位厚度测度。测度来自 retained Core domain，不能另建网格。
- D_I=B A_f 保持唯一；归一化调用同一 D_I 后只逐 cell 除 V_c，C=-D_I^T 保持积分含义，不改为归一化散度的负转置。
- 轴面 A_f=0 且 u_r=0 精确；即使其零通量使 D_Iu=0，也必须拒绝非零轴速度，不能自动改零。
- 轴反射为 u_r(-r,z)=-u_r(r,z)、u_z(-r,z)=u_z(r,z)、p(-r,z)=p(r,z)，对应轴上 u_r=0、∂r u_z=0。反射是派生边界值，不是新的 primary ghost unknown。没有选中轴面的域/XY 不允许调用轴反射。
- 不直接计算轴上 u_r/r；H06 再实现 strain 的正则极限。本项证明连续性及轴边界规则，不能据此声称黏性/曲率/完整 RZ PDE 已验证。

## 实现要求

- [x] TODO 先于源码；冻结 H04 交付身份、源码/测试/程序/PDF/case。
- [x] prepared 接口保留同一 D_I/layout，真实 V_c 与 A_f、轴 DOF 均只读；无额外归一化矩阵/散度循环/主时钟。
- [x] 体积归一化和轴奇偶规则具名，非法来源/非零轴速度/非有限反射或归一化明确拒绝，失败不污染来源或已返回结果。
- [x] 独立 edge endpoints/cell bounds 验证 XY、RZ 轴上/离轴、2:1 AMR、selected/不连通/带孔域测度和解析体积平均 divergence。
- [x] 验证径向线性流的轴邻域有限极限、奇偶延拓与轴向 Neumann、正则无散流；缩小物理半径不引入额外1/r噪声。
- [x] 均匀 RZ 三档 n=4/8/16 的非多项式 divergence MMS，物理体积加权相对 L2 递减、末档阶≥1.8。
- [x] 标准构建、前置/架构/目录/App 回归、来源与交付审计；更新 API/报告和下一项 H06。

## 验收

**T0：** local Release 无新增compiler warning；公开头自足，来源唯一CCH owner、Core/Eigen依赖不扩大，architecture/self_test/repository_layout/lamp_help通过。

**T1：** 独立端点测度相对缺陷≤1e-13；解析 polynomial 与正则无散场逐 cell 归一化误差≤1e-13（参考尺度 max(1,|准确 divergence|)），全局真实体积加权积分与独立边界通量闭合。轴奇偶参考在实际正负半径比较，轴向中心导数精确零；缩半径径向线性场 divergence=2a 有限。错layout/拉伸/顶点/AMR来源、非零轴面（含极小非零）、无轴反射/非法分量/非有限数据和归一化溢出拒绝；复制、独立等价mesh及析构后寿命正确。

**T2：** u_r=r sin(πz)，u_z=cos(πr)sin(πz)，准确 divergence=2sin(πz)+πcos(πr)cos(πz)。production 用 native face-center 值，独立参考按真实旋转体积分得到 cell averages，不读取 D_I/C 构造参考；n=4/8/16 物理体积加权误差递减、末档阶≥1.8。另以准确 face-average polynomial 通量检查测度的舍入级再现，不能把再现当作一般 MMS 零误差。H01–H04、FirstJet App、旧 gradient load回归保持。

## 失败规则

测度/轴条件/归一化/来源/解析或MMS阶失败即FAIL，范围内修复，不调容限或样本。需改Core/前置算子、完整Hodge/strain或超出预算时拆缺失contract后记BLOCKED，不以纯代数恒等式替代上述解析/PDE测度验证。

## 交付

prepared物理散度/轴规则、独立短契约、API/TODO/报告、标准构建/定向回归与原始误差/阶/哈希证据。完整时间推进、RZ pressure/velocity PDE、AMR Hodge、黏性和低黏/EHD仍待后续gate。
