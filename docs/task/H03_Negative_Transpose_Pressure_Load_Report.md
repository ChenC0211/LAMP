# H03 / 负转置压力载荷验收报告

**状态：PASS，2026-10-03。** 按[任务书](H03_Negative_Transpose_Pressure_Load_TODO.md)实现冻结 CCH §6.2 式(6.4)–(6.6)及 H-R2：C 只由 H02 的权威 D_I 自动负转置生成，离散压力功伴随恒等式通过。下一单任务 **H04：单相 saddle solve**；H03 不代替方程精度、边界/Hodge 或 Gate H。

## 生产实现

新增不可变 [CchPressureLoad API](../api/H03_Cch_Pressure_Load.md)，构造只接收 CchIntegralDivergence，保留同一来源 handle，直接执行 `C=-divergence.matrix().transpose()`。只读矩阵为 Nu×Np，compressed 且 nnz(C)=nnz(D_I)；所有非零项精确为原矩阵的负转置，未添加 geometry query、pressure-gradient stencil、AMR correction 或虚拟子面。apply 验证 H01 layout 与 actual mesh，只用 C 乘 pressure；velocity 不影响结果。

复制共享 C 与 retained source，mesh/domain/原 divergence 对象析构后仍可使用。错来源与有限 pressure 造成的载荷溢出明确拒绝，不污染算子和值，后续合法调用通过。返回的是 pressure×area 积分载荷，未除 density/dual volume、未乘 dt、未接入物理边界或求解器。

完整未消元 C 保留所选域边界的一项贡献；常压在内部面严格抵消，边界载荷不被删除。RZ 的 full 2πr 测度直接继承 D_I，原零轴列转为精确零载荷行。这些是代数/测度语义，未实现 pressure gauge 或轴速度消元。

## 独立算子验证

共 **75 场景、2,400 组 u,p、2,616 个粗侧面真实子段、25 个零轴速度模式、32 个拒绝**。覆盖 XY、RZ 轴上及离轴，每种测度均含 n=2/4/8/16/32 五档均匀网格、4×8 非方网格、n=4/8/16 的两层 AMR 全域/半域/不连通/带孔域、非二进制精确坐标 8×6 AMR 全域与带孔域。另有每种测度 depth=0/1/2/0/1 的 15 个重新绑定场景；合计 18 uniform + 36 AMR selection + 6 non-dyadic + 15 rebuild。重建测试不声称物理 velocity transfer。

每个场景 32 对输入包含 u=0、p=0、常量 1 以及独立组合的 1e-100/1/1e100 随机幅值。实际 Cp 与 D_Iu 均走生产 apply，压力功由独立 dot 累加、范数由 hypot 累加，分母保持 H-R2 的 ||u||₂||p||₂；不使用面积或其他分母稀释误差。

独立逐面参考从实际 edge endpoints 的 XY/RZ 测度和 FaceSegmentGraph 的 outward normal 计算 `-Σs_cf A_f p_c`，不调用 production domain.face、CartesianMeasure 或 D_I 系数。逐面载荷误差除以 Σ|A_fp_c|，要求≤1e-13；零尺度精确零。exact negative transpose 检查尺寸、nnz、压缩、每项精确相反和全部 incidence；内部常压精确零，未消元 boundary 常压贡献保留。

| 指标 | 实际最大值 | 运行前门槛 |
| --- | ---: | ---: |
| η_adj：离散压力功伴随缺陷 | **7.38771e-16** | H-R2 <1e-12；本验证同时 <1e-13 |
| 独立实际面载荷归一化误差 | **2.25179e-15** | ≤1e-13 |
| C 对 -D_I^T 的系数误差 | 精确 0 | 精确 0 |
| 25 组零轴模式的 D_Iu / 压力功 | 精确 0 | 精确 0 |

32 拒绝含重建旧 mesh/layout 24 次，以及 fresh same-size layout、same-count stretched geometry、实际顶点变化、XY/RZ 来源、同尺寸左右 selection、非有限压力输入、refined mesh 和有限压力溢出共 8 次。公开头独立作为首 include，const API、不可默认构造、复制共享、来源恢复后的合法调用和独立等价 mesh 生命周期均通过。MSVC long double 与 double 精度相同，没有扩展精度参考声明。

## 构建、回归与范围审计

标准 local Release 构建通过，无新增 compiler warning。新 H03 契约 **1/1 PASS（0.97 s）**；另一个回归窗口 **8/8 PASS（68.69 s）**，涵盖 lamp_help、architecture、architecture_self_test、repository_layout、interface_first_jet_app_contract、H01 layout、H02b divergence 和 hydro_gradient_load_contract。共九项不同检查，未声称一次完整 CTest。

实扫架构 14 targets / 373 translation units，自测 133 项和 31 个端到端拒绝通过。H01 原数值记录仍为 23 场景/29 拒绝；H02b 的 75 场景/2,400 字段、六项缺陷及 32 拒绝汇总逐字保持。FirstJet App 12 smokes/49 拒绝通过。旧 Hydro 的 OMP 线程上限提示保留，未修改物理参数或线程条件。

开工 1,010 份来源中 1,008 份保持；只改 root/tests 两份 CMakeLists，新增 pressure-load header/cpp 和一个契约 cpp。预算为 **3 production/build + 2 tests**。H01/H02b 实现、全部旧 production、architecture checker、冻结 CCH PDF、LAMP.exe 及 243 份 case 输入身份不变。新 .cpp 唯一 owner 为 lamp_hydrodynamics_cch，依赖闭包仍为 Core/Eigen，没有旧 Hydro/Compatibility。

本次未运行完整 73 项 CTest、Debug、完整收敛 benchmark、pressure/velocity MMS、Hodge/边界求解、实际流体或低黏/EHD 耦合。数值与编译首轮即通过，未改容限、旧散度或物理参数。

## 交付与复现

生产：[头文件](../../include/Modules/Hydrodynamics/CchPressureLoad.h)、[实现](../../src/Modules/Hydrodynamics/CchPressureLoad.cpp)与 root CMake source list。测试：[独立契约](../../tests/contracts/hydrodynamics/cch_pressure_load_contract.cpp)及 tests 注册。API、TODO/报告、现行索引/路线图同步；旧阶段报告与任务结论保留历史。

```powershell
.\tools\build.cmd -Target quantum_cch_pressure_load_contract
ctest --preset local -R '^cch_pressure_load_contract$' --output-on-failure
```

在 CTest 可用的开发 shell 中执行；通用默认检查入口仍为 `tools/build.cmd -Test`。

证据：[数值汇总](../../runs/baselines/H03_20261003/verification.json)、[新契约完整输出](../../runs/baselines/H03_20261003/operator_detail.log)、[八项回归完整输出](../../runs/baselines/H03_20261003/regression_detail.log)、[标准构建](../../runs/baselines/H03_20261003/build.log)、[来源审计](../../runs/baselines/H03_20261003/source_audit.json)、[实际依赖图](../../runs/baselines/H03_20261003/architecture_graph.txt)、[开工身份](../../runs/baselines/H03_20261003/start.json)、[最终交付清单](../../runs/baselines/H03_20261003/delivery_manifest.json)、[文档审计](../../runs/baselines/H03_20261003/document_audit.json)。临时审计脚本随证据归档，未新增长期 tools 脚本。
