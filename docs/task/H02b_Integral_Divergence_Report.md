# H02b / 唯一原生积分散度验收报告

**状态：PASS；父任务 H02 PASS。** 2026-10-03。按[任务书](H02b_Integral_Divergence_TODO.md)完成 D_I=B A_f 的唯一装配与计算，符合冻结 CCH §6.1 和 H-R1/H-R2 的散度部分。下一单任务 **H03：C=-D_I^T**；本次未实现压力载荷、流体求解或 Gate H。

## 实现与范围

新增不可变 prepared [CchIntegralDivergence API](../api/H02b_Cch_Integral_Divergence.md)，绑定 H01 的 native layout/mesh identity，唯一 owner 为真实 `lamp_hydrodynamics_cch`。每个参与原生面查询一次 Core geometry，negative/positive cell 分别装配 +A_f/-A_f；实际 fine subfaces 保持唯一 velocity DOF。内部列的两项系数精确抵消，RZ 轴上的零面积列保留但不存显式零项。apply 只调用同一 compressed matrix 乘法，不保存 B、额外几何或第二套散度。

返回 cell 积分净外流量，不作 Vc 归一化，也不加入密度、相分数、时间步或边界猜测。公开矩阵只读；算子复制共享，mesh/domain 析构不影响 retained 数据。错来源与有限输入产生的溢出明确拒绝，后续合法调用仍通过。

## 数学验证

独立参考由实际边端点和 FaceSegmentGraph 的 cell outward normal 求积分，不调用 production domain.face 或 CartesianMeasure。XY 使用端点边长，RZ 使用 2πr·dz / π(r1²-r0²)。所选域边界包含未选邻居；真实 coarse side 按相邻单元层级识别，由两条实际 fine edges 的通量核对矩阵贡献。

总计 **75 场景、2,400 组速度字段、1,308 个粗侧面聚合检查、32 个拒绝**：

- XY、RZ 轴上及离轴各五档均匀 n=2/4/8/16/32，加三组 4×8 非方网格，共 18 场景。
- 三种测度、AMR 基础 n=4/8/16、全域/半域/不连通域/带孔域，共 36 场景。
- 三种测度的非二进制精确坐标 8×6 / 两层 AMR，全域与带孔域，共 6 场景。
- 每种测度按 depth=0/1/2/0/1 重新建立网格与算子，共 15 场景；重新绑定通过，旧来源拒绝。这是算子重建，不是物理 velocity transfer。

每场景 32 组速度包含精确零、常量 1 和 1e-100/1/1e100 尺度的确定性随机字段。逐 cell 与全域边界闭合均使用运行前固定的 **1e-13** 通量相对尺度，零尺度要求精确为零。另以线性/二次向量场的准确面平均验证 cell 解析净流量，以 XY/RZ 多项式流函数验证局部零散度；这些准确测试输入不参与生产算法，也不声称实际速度重构已有相同精度。

| 指标 | 最大归一化缺陷 |
| --- | ---: |
| 独立面测度/矩阵系数 | 2.02761e-15 |
| 随机场局部积分闭合 | **1.87897e-15** |
| 全域净流与 selected 边界流量 | 2.72435e-16 |
| coarse side / 真实子面求和 | 1.89931e-15 |
| 解析非零散度积分 | 1.04492e-15 |
| 流函数局部零散度 | 4.76512e-16 |

所有矩阵 compressed、尺寸 Np×Nu，nnz≤2Nu、≤8Np，每行≤8；内部列系数和精确为零，轴列未丢失。32 拒绝覆盖 stale geometry/topology、同尺寸不同 layout、selected/RZ 错来源及输出溢出；独立等价 mesh、只读共享/寿命及溢出后的合法调用通过。MSVC 的 long double 与 double 精度相同，本报告不将独立参考声称为扩展精度结果。

## 构建、回归与审计

标准 local Release 构建通过，无新增 compiler warning。新算子契约最终 **1/1 PASS（0.78 s）**；此前独立回归窗口 **7/7 PASS（67.80 s）**，覆盖 lamp_help、architecture、architecture_self_test、repository_layout、interface_first_jet_app_contract、cch_dof_layout_contract 和 hydro_gradient_load_contract。它们是八项不同检查的两个执行窗口，不是一次完整 CTest。最后只改新测试夹具，生产与旧测试保持，未重复运行这些无变化的物理回归。

architecture 实扫 14 targets / 372 translation units；自测 133 项、31 个端到端拒绝通过。H01 保持 23 场景/29 拒绝；FirstJet App 保持 12 smokes/49 拒绝。旧 Hydro 检查的 OMP 线程上限提示仍保留，不属于新增编译 warning。

开工 1,007 份来源中 1,005 份不变，旧文件只改 root/tests 两份 CMakeLists，新增 production header/cpp 与一个测试 cpp。总预算 **3 production/build + 2 tests**。全部旧 production、H01/H02a、architecture checker、冻结 CCH PDF、LAMP.exe 及 243 份 case 输入身份保持。唯一新 .cpp owner 和 Core/Eigen 依赖闭包通过；未链接旧 Hydro/Compatibility。

首轮 partial selection 参考遗漏了未选邻居形成的边界，修正为实际 selected 边界后通过；补充非整齐坐标时，测试用严格浮点相等判断半边长，改为真实邻居层级识别及端点舍入核验后通过。两次修复只改独立测试参考/夹具，生产实现和 1e-13 数值门槛未修改。初次编译的 Vector2D 乘法调用修正及失败日志也保留。

本次未运行完整 72 项 CTest、Debug、pressure/velocity MMS、完整 benchmark 或低黏/无黏耦合。H02 仅证明拓扑积分散度；H03 伴随、H04 saddle、H05 轴正则及后续物理 gate 仍须依序完成。

## 交付与复现

修改源码：[公开头](../../include/Modules/Hydrodynamics/CchIntegralDivergence.h)、[实现](../../src/Modules/Hydrodynamics/CchIntegralDivergence.cpp)、[独立契约](../../tests/contracts/hydrodynamics/cch_integral_divergence_contract.cpp)、root/tests CMakeLists。当前 API、TODO、报告与文档索引/路线图同步；无 case 或工具入口变更。

```powershell
.\tools\build.cmd -Target quantum_cch_integral_divergence_contract
ctest --preset local -R '^cch_integral_divergence_contract$' --output-on-failure
```

在 CTest 可用的开发 shell 中执行；默认通用检查入口仍为 `tools/build.cmd -Test`。

证据：[最终数值汇总](../../runs/baselines/H02b_20261003/verification.json)、[算子完整输出](../../runs/baselines/H02b_20261003/operator_detail.log)、[7 项回归完整输出](../../runs/baselines/H02b_20261003/regression_detail.log)、[构建](../../runs/baselines/H02b_20261003/build.log)、[最终契约重编译](../../runs/baselines/H02b_20261003/operator_rebuild.log)、[来源审计](../../runs/baselines/H02b_20261003/source_audit.json)、[实际依赖图](../../runs/baselines/H02b_20261003/architecture_graph.txt)。

失败记录：[首轮 selected 边界参考](../../runs/baselines/H02b_20261003/initial_operator_detail.log)、[非整齐坐标夹具](../../runs/baselines/H02b_20261003/nondyadic_initial_operator_detail.log)、[初次编译](../../runs/baselines/H02b_20261003/initial_build.log)。交付与本地链接检查见[清单](../../runs/baselines/H02b_20261003/delivery_manifest.json)和[文档审计](../../runs/baselines/H02b_20261003/document_audit.json)。
