# I12a6 / 同源 native Q4 矩形液量报告

**PASS。** 已能从原完整 native Q4 snapshot 测量受支持矩形内的两相液量，不重建 scalar、不更改 M。它仍是空间测量；I12 整体 FAIL，后续继续扫掠积分和 driver。

## 实现

2 个新增 production 文件 `HermiteGalsPhaseRegion.h/.cpp`，新增 contract 与 CMake 注册共 2 个测试相关文件。范围见 [TODO](I12a6_Rectangular_Phase_Measure_TODO.md)，接口与限制见 [API](../api/I12a6_Rectangular_Phase_Measure_API.md)。

prepare 保留 state/view 和同源完整 graph；query 按原实际叶片切割矩形，用完整目标/补集区域维持每条 graph 的 owner。边界 cuts、free moments 与原点一致。逐 endpoint 检查覆盖，纯相使用原证书；切割区域直接积两相，稳定形成较小相的互补量。严格容量、双 Green ledger 与浮点误差估计门槛保持分离。根细化只替换请求的 retained query，重新计算其所有相关 moments；失败不发布任何修改。

## 验证

- 主程序及 11 个相关测试目标 Release 构建通过，12.73 秒，无 warning。
- **11/11 CTest 通过**，132.45 秒，其中新 region contract 42.67 秒。I07、I10/I11、I12a1–a5 和 CLI 的 **9 项完整输出与 I12a5 相同**；architecture 符合新增一个 TU/两个源文件：12 targets / 380 TUs / 899 files / 3084 includes / 5 legacy edges。
- 新合同 **212 组测量检查、11 类拒绝**。XY 圆/环形与 RZ 球/球壳，正反相，4×4 均匀与 2:1 AMR；整个域、跨叶片、局部矩形、1/1024 薄条、纯区域、单 ulp 宽纯相、按两方向分割可加性均与独立解析交集一致。两个 raw Green 公式及严格容量/互补一致性均检查。
- source identity/snapshot owner 保持，临时 state 生命期、复制/移动、精确 mesh clone、错误/陈旧 mesh、非法/域外矩形、无接触点的 partial active domain、正长度零迹拒绝及失败后的同源查询通过。全部源 M 保持原值。
- 81 冻结/参考文件、3 PDF、前序 I12a5 archive 哈希不变。旧 crossing probe 完整输出不变且仍退出 1；未重跑长程 benchmark，未更新 I12/Gate I 状态。

测试准备过程中修正了混合 auto 声明编译错误；RZ 制造解改为满足原轴偶对称约束的轴心球/球壳。pure partial-domain fixture 在 prepare 阶段即被原 Core 的显式接触规则拒绝，测试和 API 如实记录，没有为通过测试移除这条限制。原始失败与成功日志均在证据目录保留。

证据：`runs/baselines/I12a6_20261001/`；执行 `quantum_interface_phase_region_contract.exe` 复现完整矩阵，`--contracts` 单独复现来源/边界拒绝。目录包含完整日志、预算/冻结核验、前后输出比较、runtime 与 archive 逐项 SHA256。

## 后续限制

本接口测量调用者给定的 double 矩形，不证明其边界等于精确扫掠映射，也没有 BC 时间历史身份。浮点 moment estimate 不是严格误差界；其容限不能替代库存的 64 eps 门槛。下一任务必须明确这些时间积分语义，再将 curved Q4 接入 I12a3/a4 和完整半/全步推进。
