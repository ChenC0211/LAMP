# I12a11r8 / 闭合与扫掠共享局部体积测量

## 目标及范围

只解决一个数值事实：I07/I10/I11 的 native Q4 单元体积与 I12 PhaseRegion 对同一完整单元的体积必须采用相同的局部 Green、小相直接积分及补量路径，消除 r7 已证实的两种测量差额。

允许修改 4 production：新增 `Detail/HermiteGalsRegionMeasure.h`、`HermiteGalsPhaseRegion.cpp`、`InterfaceCellMeasure.cpp/.h`；2 tests：既有 phase_region_contract 与 benchmark。另含 TODO/报告、相关 API、README 和 `runs/baselines/I12a11r8_20261001/`。不改 Core、Hydro、M、B/J、解析参考、预算、拓扑证书、scalar 根或求积/迭代上限；冻结文件不变。

## 数学约束与实现

从现有 PhaseRegion 提取只读局部积分 helper，保持原算术顺序、retained cuts、上下相积分、原点及 root refinement，原区域查询结果不得改变。完整 Q4 discovery 在同一 scalar 和完整 topology 上使用此 helper 计算 cut cell 的 authoritative 体积和 moment estimate；纯相继续采用既有证书语义。只对现行 discovery 入口启用，legacy evaluateVolumes、敏感度积分及原全局 fixed/free 账本路径保持。

全局 fixed/free moments 仍用原共同原点，供边界反对称/账本诊断；局部负相体积与总量单独计算，明确两条代数等价路径的浮点结果可不同，不能伪造账本字段。局部体积必须有限且严格有界，双 Green/容量一致性保持既有区域要求。新增计算成本需记录，不在本项改求积算法。

## 验收与失败规则

T0：Release 无 warning，30 项直接回归通过；逐项核对输出变化，仅允许本次数值路径造成的几何/求根输出及新增测试，原区域查询结果逐字一致；源码范围、81 冻结/参考文件、3 PDF、r7 archive 审计。

T1：XY/RZ × uniform/AMR、圆/环与液气互换的 cut cell 上，将 I07 结果与相同 scalar 的 whole-cell PhaseRegion 测量逐位比较；保留原独立解析区域精度、加性、owner/M 与拒绝测试。新的显式推进模式须完成原 r7 两个四步案例，不把诊断重现当推进 PASS。

T2：T1 通过后运行原 XY AMR n=32 四步、四组静止及四组 40 步长平移、12 组非多项式空间矩阵；保持逐级下降、梯度阶≥2、曲率阶≥1.5。推进或收敛仍失败则记录 FAIL，按具体失败继续分解，不调预算或跳过配置。本项 PASS 也不能代替旋转/vortex 等 I12 其余门槛。

T2 失败定位补充：本轮 XY 空间矩阵已确认误差随细化增长，额外运行既有 `--curvature` 初始静态非多项式矩阵（相同 12 配置、无源码改动）。对比 t=0 与 t=.04 的梯度/曲率，区分初始表示不收敛和推进过程放大；该对照不能替代运动矩阵或放宽门槛。若与主矩阵并行运行，记录运行时间受并发影响。

交付修改文件、原始日志、旧新输出比较、运行成本、未解决问题、源码/runtime 哈希和 archive。按完整证据记录 PASS/FAIL；只有无法在已有授权/依赖内继续时才报告阻塞原因。
