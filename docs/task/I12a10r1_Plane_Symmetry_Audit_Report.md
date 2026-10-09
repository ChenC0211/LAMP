# I12a10r1 / AMR 平面对称性丢失定位报告

**PASS：诊断解释全部 12 个拒绝；I12a10/I12 仍 FAIL。** 遵守 [TODO](I12a10r1_Plane_Symmetry_Audit_TODO.md)，仅改 benchmark 与 CMake，0 production。

`--plane-symmetry-audit` 重放 I12a10 全部 48 个配置。旧已发布 Q4 的外边界法向 DOF、预测 primary 的 vertex/edge 法向 first jets 都严格为零。真实 departure 上的旧 scalar 法向导数和 B/J 拉回差值也都为零；常速度 Jacobian 逐项等于单位阵。

全部 12 个 AMR 正向案例在 I06 mixed 后首次出现非零边界 dxy；对应 native Q4 边界 normal-derivative 的五个局部 DOF 不再全零，独立 native component discovery 重现相同拒绝。日志逐项列出实际边界、vertex、cell 与 slot，不再将 bubble 异常中的 cell=27 当作故障位置。例如 RZ AMR 正向、正相、q=.75 的 r=1 上 vertex 14 的 dxy 为 -8.4080890398278491e-15，对应 cell 5 slot 18 为 -5.2550556498924057e-16；局部张量 DOF 包含单元尺寸缩放，两者不是同一单位。

在独立 private Patch 中，仅将声明 symmetry 上的派生 mixed 置零，保持 primary、center 和所有其它 mixed 不变。重新生成 exact hanging 约束后，48 个案例的边界 Q4 导数全部严格为零，12 个原拒绝全部恢复 native discovery。该快照没有发布、没有传给 advance，也未更改实际 candidate 或 accepted identity/M，因此只证明故障层及约束效果，不代替完整推进验收。

代码检查显示 I06 已对 RZ 轴采用精确 mixed 约束，但 API 只接收 candidate，无法获知 I07/I10/I11 使用的外边界 symmetry 声明。下一独立修复应向现有 I06 兼容系统传入同一声明，对严格满足零法向 primary 的真实外边界施加派生 mixed 零约束；不能在求解后盲目抹零、容差接受非法输入或修改 first jets。

最终 Release 构建退出 0、3.46 秒、无 warning。11 项 CTest 中 10 PASS；唯一 FAIL 仍是原 crossing publication 的 36/48，退出 8，32.42 秒。全部 10 项旧测试完整输出逐字不变，旧 midpoint crossing probe 仍退出 1 且输出不变。其它 27 个可执行、所有 production、81 冻结/参考、3 PDF 和 I12a10 archive 哈希保持。完整日志、审计、源码及 runtime 位于 `runs/baselines/I12a10r1_20261001/`。
