# BF-01 / Regularized color 的唯一只读派生报告

状态：**PASS（仅 BF-01）**。本轮只实现同一 Pure FirstJet publication 的只读 CSF color 与求积审计。BF-02 尚未启动；完整 BF-CSF 的静滴、无黏波、能量/谱和强变形能力仍待后续 gate。旧 H13 FAIL 和 BF-00 已登记的 momentum 契约 TIMEOUT 保持原判定。

任务格式和范围见 [BF-01 TODO](BF-01_TODO.md)，接口见 [BF-01 API](../api/BF-01_Cch_Regularized_Color.md)。用户关于旧数据的授权已落实为当前冻结源码的新 A/B 基线；本任务没有恢复或依赖已丢失 payload。

## 修改范围与唯一 owner

生产文件共 3 个：新 `include/Modules/Hydrodynamics/CchRegularizedColor.h`、新 `src/Modules/Hydrodynamics/CchRegularizedColor.cpp`、根 `CMakeLists.txt` 的单行 CCH 注册。测试文件共 2 个：新 `tests/contracts/hydrodynamics/cch_regularized_color_contract.cpp`、`tests/CMakeLists.txt` 的单行契约注册。满足默认 5/2 预算。

新对象保留不可变 `CchPhaseGeometry` 来源和统一的物理 `epsilon`，输出 native pressure 顺序的 `color_csf`、逐 cell 报告和全域诊断。身份校验复用原 mesh/layout/geometry/full-clock 契约，并检查同一 divergence matrix owner。没有第二份几何、独立输运状态或材料库存。

与 BF-00 的 1,397 份数值源码/输入逐项 SHA256 比较：既有文件只有两份 CMake 注册变化，其余逐字保持；仅新增上述类对和契约。Core/Interface、sharp 相测度/材料质量、压力/鞍点、唯一 F_rho、实际 R/uR、对流/黏性及旧 H12/H13 实现没有修改。详见 [来源差异](../../runs/checks/BF-01/20261009_143534/source_diff.json) 和 [注册 diff](../../runs/checks/BF-01/20261009_143534/registration.diff)。当前树身份为 `cfc713745564efb00b5a658c1baf0e32f196a874252e7bc8187b33e8cbc9415c`（按路径、NUL、文件 SHA256、换行排序拼接）。

## 数学实现与拒绝条件

逐点只读查询原 Q4 的 phi/gradient，带内使用有长度量纲的 `d_h=phi/|grad phi|` 和规范 `H_epsilon(-d_h)`。原标量未被重置为距离函数；正缩放在可表示条件下消去。H 端点采用等价稳定级数，delta 先计算 s/epsilon，避免巨大有限长度导致中间乘积溢出；不可表示的 delta 明确拒绝。

XY/RZ 子矩形采用真实物理容量和 GL8/GL16 正权平均；RZ 权重 2*pi*r 只乘一次。companion liquid/gas 权重按真实容量归一，精确再现常值并保持 boundedness，未 clip color 或修补 sharp 材料。整子矩形认证为带外时直接返回精确 0/1。

原 Q4 affine 包络和精确四次多项式导数恒等式给出整 panel 梯度上下界；平移仅查询同一多项式数学延拓，TwoSum 保留坐标舍入余项，没有物理 ghost/邻格数据。纯相和带判别向外取整；H 上下界由 16/17 项递减交错级数构造，范数采用缩放后的有向商、平方、和、sqrt 与乘积。

所有非 pure 接受 panel 必须证明 `gradientLower > up(relativeGradientResolution*gradientScale)`。未认证几何优先细分，与积分误差预算分开，不能被零 GL 差值放行。统一带宽若不满足 `epsilon>=max(hx,hy)/2`，cut cell 早拒绝；其他可能与带相交的子矩形必须细分，直至认证 pure 或明确拒绝。查询和深度耗尽均报告失败。查询上限同时计原 jet 和 affine 包络，默认每 cell 262144。

absoluteErrorEstimate 是 GL8/GL16 差值及全纯采样时 color 跨度的自适应**积分误差估计**，不是真误差的严格上界。带体积采样估计与区间上下界分别保存；全 panel 梯度下界比和实测点最小梯度比分开报告。保守认证可能拒绝尚不能证明安全的输入；不把这种拒绝冒称为已证明原几何退化。

## 实际构建与验收

工具链为 VS2022/MSBuild 17.11.9、MSVC 19.41.34120.0、CMake 3.29.5-msvc4，x64 Release；构建及全部输出在 D:。新增类由真实 lamp_hydrodynamics_cch 单一 owner 编译。最终 [Release 构建日志](../../runs/checks/BF-01/20261009_143534/build_color_delivery.log) 无新 warning/error。

| 验收 | 实测结果 | 门槛/证据 |
| --- | --- | --- |
| 新注册 color 契约 | PASS，123.714 秒 | 原 180 秒不变；[完整日志](../../runs/checks/BF-01/20261009_143534/color_delivery_ctest.log) / [JUnit](../../runs/checks/BF-01/20261009_143534/color_delivery_ctest.xml) |
| kernel | PASS | 401 点单调/有界/独立公式/导数、端点、4096 点 delta 积分、极大/极小物理长度 |
| 平界面位置扫掠 | 264 cell，最大平均误差 3.62514e-9 | 原 4e-9；XY/RZ、uniform/2:1 AMR、三位置、重复逐位稳定、正标量缩放 |
| epsilon 敏感性 | 176 额外 plane cell | 同 publication，epsilon/h=1,1.5,2；真实 color、sharp 差异和带体积保留 |
| 真四次多项式 | 44 cell，最大平均误差 2.24830e-10 | 原 4e-9；独立解析 phi/gradient 的 1024/2048 中点 Richardson，3.7 正缩放 |
| 圆/球 | 176 cell，最大独立误差 2.34841e-10 | 原 2e-7；同 publication 两 epsilon、1e-9/1e-10 两积分预算、独立 128/256 积分 |
| RZ near-axis | PASS | 轴上 proxy 与解析值及偶对称近轴极限；无 1/r 查询 |
| native face 低梯度反例 | 正确拒绝 | phi=(x-.5)^3+1e-5*(x-.5) 原几何合法；color 阈值 2e-4，Gauss 点未命中的 face 最小梯度触发未解析认证拒绝 |
| selected/来源/clock/退化 | PASS | 同 publication 子域 color 与全域相同；错误 layout/geometry/full clock、epsilon 和梯度配置拒绝；sharp 材料不变 |
| 常值纯相 | XY/RZ、正负四组合 PASS | 梯度零且带外，全部精确 0/1，无点采样 |
| T0 和既有相关回归 | 13/13 PASS，312.22 秒 | architecture/self-test/layout、native DOF/D/C/saddle/测度、相几何/黏性/材料、旧恒应力/变量牵引；[日志](../../runs/checks/BF-01/20261009_143534/regression_ctest.log) / [JUnit](../../runs/checks/BF-01/20261009_143534/regression_ctest.xml) |

完整原始 stdout 与 [660 行 cell CSV](../../runs/checks/BF-01/20261009_143534/color_cells.csv) 已保存；每行包含 epsilon、c、sharp alpha、带估计/上下界、输入标度、两种梯度比、误差估计及成本。最终记录中最大单 cell 总 source 查询为 30665、最大细分深度 12；均在原上限内。汇总见 [color_summary.json](../../runs/checks/BF-01/20261009_143534/color_summary.json)。

开发中的失败日志未删除：薄带漏采、曲面局部预算/查询上限和制造解先被原 publication 拒绝的版本均留在本时间戳目录；最终通过来自修正实际实现或合法反例初值，没有扩大 epsilon、放宽阈值、改 phi/消费者或延长 CTest。较早成功 JUnit 默认截断 stdout；最终用 CTest verbose 和更大的输出保存长度归档完整数据，物理/时间门槛保持。

## 未解决项与下一任务

BF-00 的完整 cch_momentum_geometry_contract 两次 180 秒 TIMEOUT 属于冻结的旧模态候选，未因本 color 改动消失，也没有以分项诊断替代 PASS。本轮数值代码未修改该候选，按原记录单列。

BF-01 只证明新 color owner 在上述范围的公式、平均、来源和拒绝行为。尚无 CSF 载荷、真实曲率转移、完整静滴/无黏周期/能量谱或生产时间推进；总体目标未完成，H14/H15 不启动。按已授权 BF 顺序，**下一次允许仅执行 BF-02**：复用同一 C 构造恒曲率实际 RHS；本轮未写 BF-02 代码或 TODO。
