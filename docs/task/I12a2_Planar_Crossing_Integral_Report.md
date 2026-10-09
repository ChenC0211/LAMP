# I12a2 / 平面跨面事件积分报告

**PASS。** 完成常速度、同源平面、单一 native face 的时间区间积分。父任务 **I12 仍 FAIL**：本任务没有接库存、reflux、bubble 或 benchmark 驱动。用户已授权继续推进至 I12 通过，归档后进入下一独立子任务。

## 实现与范围

修改 5 个 production 文件：`HermiteGalsTemporalPhaseFlux.h/.cpp`、`LevelSetHermiteGeometry.h/.cpp`、`Detail/LevelSetHermiteQuarticLeafData.h`；2 个测试相关文件：新增 `interface_planar_temporal_flux_contract.cpp` 与 CMake 注册。详见 [TODO](I12a2_Planar_Crossing_Integral_TODO.md) 和 [API](../api/I12a2_Planar_Temporal_Phase_Flux_API.md)。

常速度来源工厂同时持有值和生成实际 endpoint queries，复用现行 RK2。Core 在 strict TU 内从原 Q4 homogeneous DOFs 证明完整扫掠盒中所有叶子为同一物理平面；精确裁剪 face×time 上的负半平面，用有理多边形面积与一阶矩积分。RZ 使用真实 2*pi*r 权重。没有拟合近似平面、重建时间节点解析状态、M clipping、容差/CFL 修改或 midpoint fallback。

正反方向独立写入四个 gross integrated volumes，rate 只由区间派生；来源、旧 authority、实际面与区间组合保持 I12a1 语义。近纯液比例通过同一有理测度的气相补集求值；事件和多项式证明使用精确谓词，最终浮点输出仍有表示舍入。

## 验收证据

- Release 主程序与 20 个相关测试目标构建成功，261.11 秒，无 warning。
- 19/19 CTest 通过，107.52 秒。17 项前序测试完整输出与 I12a1 相同；architecture 仍为 12 targets / 379 translation units / 896 files，内部 include 从 3074 增为 3075，5 条 legacy edges 不变。
- 新算子共 741 组检查：730 组穿越矩阵（73 个相位、正反速度、XY 及四个 RZ 环带含轴）、原十进制反例、斜面一阶矩、静态截面、端点瞬时零迹、restart 子区间、轴面与静止零通量；包含同源非均匀时长区间组合。
- 原 q=1/4、3/4 返回 `0.003749999999999999`、`0.0012499999999999968`，与独立闭式参考在既有 double 舍入尺度内一致。
- 8 类拒绝验证：正时空测度恒零 trace、未认证 callback 速度、非 affine Q4、未知 coverage callback、域外扫掠、RZ 常径向速度、非有限常速度、陈旧 mesh。正常和失败调用均不改变旧 state、snapshot 或 M。初始 state 为结构 fixture；这些检查不代表库存推进或闭合轨迹。
- 69 个冻结文件、12 个集成参考与 3 份 PDF 身份一致，上一 I12a1 archive 未变。
- 原 `--crossing-probe` 仍退出 1，两项原拒绝及完整输出与 I12a1 一致。未重复完整 I12 长程 benchmark，因为本任务没有接驱动。

可复现命令：

```powershell
build/windows/tests/Release/quantum_interface_planar_temporal_flux_contract.exe
build/windows/tests/Release/quantum_interface_hermite_gals_benchmark.exe --crossing-probe
```

完整目标、构建/CTest 日志、source audit、runtime SHA256、旧 crossing 与逐输出比较在 `runs/baselines/I12a2_20261001/`。`evidence.zip` 和 `archive_verification.json` 给出逐条归档校验。审计脚本首次把新增系统 `<optional>` 也计入内部 include，预期 3076 不符；核对工具定义后修正为 3075，仅修审计预期，没有改测试或数值代码。

## 未解决与下一项

一般曲面、变化速度、未知 callback 的事件完备性、外部边界历史仍不在本算子能力内。单面积分不能代替严格库存、拓扑/曲率、动态 AMR 或长程验收。下一项单独实现 integrated-volume 共享库存装配，再接 I09 消费和完整步进；每项仍先写 TODO，并保持原有严格门槛。
