# H03 / 从唯一散度生成负转置压力载荷

**状态：PASS，2026-10-03。** 结果见[验收报告](H03_Negative_Transpose_Pressure_Load_Report.md)，下一单任务 H04。

依据：[现行路线图](../roadmap/05_LAMP_Module_Task_Roadmap.md) H03、H-R2 和冻结 CCH 方法 §6.2 式(6.4)–(6.6)。前置 H01/H02a/H02b 已 PASS。本项只解决一个数学事实：压力积分载荷 C 必须由同一 D_I 的负转置生成，满足离散压力功恒等式。

## 目标

在真实 CCH 库中建立不可变 CchPressureLoad，保留 H02b 的权威散度 handle，通过代码关系 C=-D_I^T 准备只读稀疏矩阵；验证 XY/RZ/2:1 AMR 的 random u,p 伴随关系，供 H04 saddle solve 直接消费。

## 允许修改

- `include/Modules/Hydrodynamics/CchPressureLoad.h`、`src/Modules/Hydrodynamics/CchPressureLoad.cpp`：负转置生成、同源只读消费与 apply。
- root `CMakeLists.txt`：新实现加入已有 CCH 唯一 source list。
- `tests/contracts/hydrodynamics/cch_pressure_load_contract.cpp`、`tests/CMakeLists.txt`：新独立契约与标准注册。
- 本 TODO、API、报告及当前索引/路线图。结果归 `runs/baselines/H03_20261003`，临时审计归 `runs/checks/H03/<时间戳>`。

预算：3 个 production/build 文件 + 2 个测试文件；不修改 H01/H02b 算式。

## 禁止修改

Core、Solvers、旧 Hydro/Compatibility、Interface/Electricity/Plasma/App、case/字典、冻结 PDF、原测试/验收容限、CCH 布局和散度实现。不得实现 H04 saddle、H05 边界/轴消元、H06 黏性或后续质量/动量/界面力/时间推进；不得独立装配 pressure-gradient stencil、补 AMR correction、插值 virtual children 或引入 fallback。

## 数学约束

\[
C=-D_I^T,\qquad u^TCp+(D_Iu)^Tp=0.
\]

- C 为 Nu×Np，所有非零系数逐项精确为同一 D_I 的负转置；nnz(C)=nnz(D_I)，不额外查询面几何或保存 stencil。
- 内部 native face 的载荷为 A_f(p_positive-p_negative)，实际 fine subfaces 不被替代。selected 边界仅有一项负转置贡献；constant p 不因此成为完整未消元矩阵的零模。
- RZ 使用 D_I 已有 full 2πr 测度，轴的零散度列转为 C 的精确零行，速度 DOF 仍保留。
- C 是压力积分载荷，不乘 dt、不除 density/dual volume、不加入物理边界；压力加速度及 AMR 高阶精度属于后续 Hodge/saddle 任务。
- apply 必须验证 values 的同一 layout identity 和实际 mesh identity，错来源/非有限结果拒绝；不可变来源和值不污染。

## 实现要求

- [x] 保存开工源码/测试/PDF/程序和 case 身份；复核前置交付及冻结公式。
- [x] C 仅由 retained divergence.matrix() 的负转置生成，公开 retained divergence/layout 与只读 compressed C；复制共享，生命周期不借用 mesh。
- [x] apply 只使用 C 乘同源 pressure，velocity 改变不影响结果，有限输入造成输出溢出明确拒绝。
- [x] 独立实际边端点/有向图核对逐 face 载荷；检查内部常压抵消、selected 边界、RZ 轴零行、粗细实际子面及 exact negative transpose。
- [x] XY/RZ 轴上及离轴、均匀五档、非方/非二进制精确域、三档 AMR、半域/不连通/带孔域及网格重建，通过随机/零/常量与大/小幅值 u,p；拒绝 stale 来源。
- [x] 标准构建、T0/T1/T2、来源和交付审计通过；报告未运行物理项并将路线图下一任务更新为 H04。

## 验收

**T0：** local Release 构建，无新增 compiler warning；公开头作为测试首 include。architecture/self_test、repository_layout、lamp_help PASS，新增 .cpp 只有 CCH owner，依赖仍为 Core/Eigen，无旧 Hydro/Compatibility。

**T1：** 75 场景，每场景 32 组 u,p。硬门槛 η_adj=|u^TCp+(D_Iu)^Tp|/(||u||₂||p||₂)<1e-12；本验证同时要求 η_adj<1e-13，零尺度要求压力功和精确为零。独立逐 face 载荷误差 / Σ|A_fp_c|≤1e-13，零尺度精确零；内部常压严格零、轴行严格零、C 与 -D_I^T 逐系数精确一致。范数和压力功使用有限的独立缩放/累加，不能换分母以降低缺陷。

**T2：** uniform/refine/coarsen 重新准备后同一契约通过，旧 mesh/layout 绑定拒绝；这不是 velocity transfer。原 H01/H02b 契约、hydro_gradient_load_contract、interface_first_jet_app_contract 保持；旧 production、runtime、case/PDF 哈希保持。未通过 H04/H05，不声称方程精度、physical gauge、低黏运行或 Gate H。

## 失败规则

负转置、伴随门槛、selected 边界/零轴、只读来源或溢出拒绝失败即 FAIL，按本范围修复，禁止调整散度/容限/物理参数。若必须改前置/Core 或超过预算，拆出明确缺失 contract 后记 BLOCKED；不以旧 gradient 替代。

## 交付

唯一自动生成 C 及 apply、独立短契约、API/TODO/报告、构建/数值/回归/依赖证据和最终哈希清单。H03 通过仅允许进入 H04；整体流体与电喷雾目标继续进行。
