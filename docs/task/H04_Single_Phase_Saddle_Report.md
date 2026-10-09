# H04 / 单相封闭域 CCH 鞍点验收报告

**状态：PASS，2026-10-03。** 按[任务书](H04_Single_Phase_Saddle_TODO.md)完成冻结 CCH §6.4/§11.1 的封闭域单体 reference solve，均匀 XY 的 pressure/velocity 连续 MMS 通过。下一项 **H05：RZ true measure / axis regularity**。本项没有实现 AMR compatible Hodge、完整 production 边界、时间推进、低黏流体 App 或 Gate H。

## 方程和实现

[CchSinglePhaseSaddle](../api/H04_Cch_Single_Phase_Saddle.md)保留 H03 权威 pressure load 及同一 D_I。caller 提供的 layout-bound A_u 经复制、有限性、symmetry≤1e-12、LLT SPD 与有限 factor 检查；solver 不推断材料/时间步，也不产生第二套 pressure gradient。封闭域以 C 行 incidence 编号 free/prescribed velocity，同时消元 A_fb u_b 与 C_b^T u_b；RZ 零轴行固定零。

按每个内部面连通的 pressure component，以 retained 真 V_c 装配体积加权零均值增广约束。稀疏矩阵中的 C_f 和 C_f^T 直接同源，压力不采用 cell gauge pin；λ 仅为约束未知量，最终原生场仍为 Nu+Np。一次 SparseLU reference solve 返回速度与压力，随后核验完整增广残差、原 free momentum、全部 D_Iu 和每个域 pressure mean。没有物理 projection、显式 inverse、人工黏度或 fallback。

模块公开 typed velocity block/RHS、prepared solver、只读 native result/report。矩阵和 factor 可供重复 RHS，复制共享，不借用 mesh/domain/input 生命周期；只有 layout/actual mesh 同源的数据可求解。固定速度处的动量反力不作为 free momentum 误差；所有所选域边界必须显式给定 normal velocity，未实现 open pressure/traction anchor。

## 独立数学验证

**39 个 algebra 场：** XY、RZ 轴上/离轴；n=4/8、depth=0/2、全域/两块不连通/带孔域共36组，另三种测度各一组非对角 SPD 与非零兼容 prescribed velocity。参考由实际 edge endpoints 的测度、outward normal 和多项式 streamfunction 构造离散无散场，pressure 独立按真 cell volumes 逐域去均值；RHS 的压力载荷来自独立边图，未读取数值 C。A_u 是该验证显式选择的 SPD block，恢复结果仅证明 solver algebra，不声称它是生产 AMR Hodge。

| 指标 | 最大归一化缺陷 | gate |
| --- | ---: | ---: |
| algebra velocity recovery | 6.31841e-15 | ≤1e-9 |
| algebra pressure recovery | 2.33114e-14 | ≤1e-9 |
| 零边界速度能量/压力零功关系 | 3.60608e-15 | ≤1e-10 |
| 全 cell 原生连续性 | **6.66685e-14** | ≤1e-11 |
| 完整增广代数残差 | 5.06463e-15 | ≤1e-11 |
| 逐连通域体积加权压力均值 | 1.11961e-15 | ≤1e-11 |

以上 report 峰值覆盖实际 algebra 与 MMS 求解的相应检查；恢复误差和能量仅取对应 algebra 场。39组 algebra 的增广矩阵 symmetry 检查精确为零；连通/两块不连通压力 gauge 数正确。非零边界 + 非对角 A_u 检查了消元对 momentum RHS 的完整贡献，重复 RHS 逐字确定性保持。

**21 个拒绝：** 非 SPD/非对称/非有限/错尺寸 A_u、非法 tolerance、错 velocity/RHS layout、错 RHS 长度/NaN、free face prescribed payload、边界净流量不兼容、同尺寸拉伸/顶点变化/AMR stale mesh、有限极大 RHS、非零轴速度、过严代数与连续性门槛，以及两块域各自非零但全局精确抵消的边界流量。最后一项在 component compatibility 阶段拒绝，证明不能以全局零流量或增广 λ 掩盖局部不可压缩失败。失败后合法零输入可继续，已返回的 result 不变；独立等价 mesh 生命周期和 caller matrix 隔离通过。

## 连续 pressure/velocity MMS

独立连续方程为 a u+∇p=f、div u=0，a=2.5，单位 XY 域全边界 normal velocity=0；ψ=sin(πx)sin(πy)，u=(ψ_y,-ψ_x)，p=cos(πx)cos(2πy)。输入为准确 native face averages / cell-average pressure；A_u=a H_u，H_u 为均匀网格真实 dual volume。forcing 是解析 face-average f 经同一 uniform Hodge 投影，不由离散 C/D 或求解结果反造。

完整 native velocity 与全部 cell pressure 分别用 dual-volume/cell-volume 加权相对 L2 比较，边界不裁除。两种误差全部档位递减，末档独立重算阶与程序输出相符：

| n | velocity 相对 L2 | pressure 相对 L2 |
| --- | ---: | ---: |
| 4 | 3.94301e-2 | 1.92771e-1 |
| 8 | 9.77858e-3 | 4.47619e-2 |
| 16 | **2.43989e-3** | **1.09893e-2** |
| 8→16 观测阶 | **2.00281** | **2.02617** |

两项满足运行前 ≥1.8 门槛。n=16 的速度/pressure 误差约0.244%/1.10%，不能以阶数代替误差大小。短三档 operator MMS 不属于论文完整空间/时间研究，也没有验证 RZ/AMR PDE 或 H15 local Hodge。

## 构建、回归和身份

local Release 构建无新增 compiler warning，新契约最终 **1/1 PASS（0.32 s）**；另一个定向回归窗口 **9/9 PASS（69.82 s）**，含 lamp_help、architecture/self_test、repository_layout、FirstJet App、H01/H02b/H03、旧 hydro_gradient_load。共十项不同检查，未声称完整74项 CTest。

architecture 实扫14 targets / 374 translation units，自测133项/31端到端拒绝通过。H01的23场景/29拒绝及H02b/H03完整数值汇总逐字保持，FirstJet App仍12 smokes/49拒绝。旧 Hydro OMP 线程上限提示保留，未调物理/线程参数。

开工1,013份来源中1,011份不变；旧文件只改root/tests两份CMakeLists，新增header/cpp与一份短契约。总预算 **3 production/build + 2 tests**。前置H01–H03、全部旧production、architecture checker、LAMP.exe、冻结CCH PDF及243份case输入身份保持，新.cpp仅CCH一个owner，依赖仍为Core/Eigen，无旧Hydro/Compatibility。

数值/编译首轮通过；后续只补新测试的component不兼容与严格残差拒绝，澄清forcing是face-average投影。生产方程、default gate及MMS参数未改。未运行Debug、完整CTest/benchmark、时间收敛或完整流体/EHD；reference LU尚不代表大规模预条件器性能。

## 交付与复现

[公开头](../../include/Modules/Hydrodynamics/CchSinglePhaseSaddle.h)、[实现](../../src/Modules/Hydrodynamics/CchSinglePhaseSaddle.cpp)、[短契约](../../tests/contracts/hydrodynamics/cch_single_phase_saddle_contract.cpp)、root/tests CMakeLists、当前API/TODO/报告和索引/路线图已同步。短测试几何只在tests生成，未新增长期case分支或tools脚本。

```powershell
.\tools\build.cmd -Target quantum_cch_single_phase_saddle_contract
ctest --preset local -R '^cch_single_phase_saddle_contract$' --output-on-failure
```

在CTest可用的开发shell中执行；通用入口仍为tools/build.cmd -Test。

任务专属的[来源/数值审计脚本](../../runs/baselines/H04_20261003/engineering_scripts/source_audit.ps1)和[交付审计脚本](../../runs/baselines/H04_20261003/engineering_scripts/delivery_audit.ps1)随证据归档；只核对来源、已有日志、独立误差阶及文档链接，不进入生产或通用tools。最终文档更新后单独复跑repository_layout，原九项回归的完整日志保留。

证据：[数值与独立阶审计](../../runs/baselines/H04_20261003/verification.json)、[完整新契约](../../runs/baselines/H04_20261003/operator_detail.log)、[九项回归](../../runs/baselines/H04_20261003/regression_detail.log)、[构建](../../runs/baselines/H04_20261003/build.log)、[契约重编译](../../runs/baselines/H04_20261003/operator_rebuild.log)、[来源审计](../../runs/baselines/H04_20261003/source_audit.json)、[实际依赖图](../../runs/baselines/H04_20261003/architecture_graph.txt)、[最终交付清单](../../runs/baselines/H04_20261003/delivery_manifest.json)、[文档审计](../../runs/baselines/H04_20261003/document_audit.json)。
