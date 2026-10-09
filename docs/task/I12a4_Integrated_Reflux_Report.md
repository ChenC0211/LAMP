# I12a4 / I09 区间积分消费报告

**PASS。** I09 已直接持有 I12a3 结果并汇总其最终 actual subface integrated volumes，库存与源 owner 不变。I12 整体仍 FAIL：一般曲面/变化速度时空积分及完整 driver 尚未接线。

## 实现

2 个 production 文件 `HermiteGalsReflux.h/.cpp`，2 个测试相关文件 `interface_integrated_reflux_contract.cpp` 与 CMake 注册。见 [TODO](I12a4_Integrated_Reflux_TODO.md) 和 [API](../api/I12a4_Integrated_Reflux_API.md)。

新增 integrated result/side/report 与重载入口；旧 I09 rate API 保留。两条入口共用同一个 native incidence、side ordering、2:1 配对与粗细诊断实现。新区间结果明确保留真实来源，四个 gross volume 与 signed volume 均直接累加最终实际子面，rate 仅由该源 dt 派生。没有重新积分、再次修复 M、复制 snapshot 或伪造 midpoint identity。

## 验证

- 主程序与 8 个相关测试目标 Release 构建通过，11.83 秒，无 warning。
- 按 TODO 运行 8 项直接相关 CTest，全部通过，12.65 秒。旧 I08/I09、I12a1/a2/a3 与 CLI 的 6 项完整输出与 I12a3r1 相同；architecture 保持 12 targets / 379 translation units / 897 files，内部 include 3078→3079，5 条 legacy edges 不变。未重复其它未修改算子在 I12a3r1 已通过的检查。
- 新合同 **33 组检查、8 类拒绝通过**：XY/RZ、0/1/2 层 AMR、四方向 coarse/fine、纯/空相与 counterflow、12 组实际平面 I12a2→I12a3→I09、修复后的 coarse/fine 实际子面、源生命期、精确 mesh clone 与陈旧 mesh、非法 side row、moved-from handles、不可表示平均 rate。
- 独立 oracle 通过扫描 actual face owner 与方向形成各 side，不使用生产 FaceSegmentGraph 分组；逐字段验证有序和完全相等。source 的 faces/rawFaces/M 与输入共享 owner；几何与 interval 保持。每个 native face 的双侧 incidence 和 side-based M ledger 均通过。
- 69 冻结文件、12 集成参考、3 PDF 以及上一 I12a3r1 archive 哈希不变。
- 原 `--crossing-probe` 仍退出 1，完整输出不变。旧长程 driver 尚未接入新通量，不宣称完整 I12 通过。

首次构建发现向旧 reflux 公共头引入 trace 定义会与旧测试的 Boundary 别名冲突；以新来源类型前置声明、仅在实现/新消费者引入 temporal header 解决。未改旧测试、Core tracer 或数值方法；原始失败与成功日志保留。

复现：`build/windows/tests/Release/quantum_interface_integrated_reflux_contract.exe`。构建目标、完整日志、逐输出比较、source/runtime 哈希和逐项 archive 校验在 `runs/baselines/I12a4_20261001/`。

## 后续

平面单面积分→共享库存→I09 汇总的算子链已通过。接下来必须覆盖实际 benchmark 的 curved native Q4 与运动时空源，再让半/全步 driver 消费；保持原 strict inventory、bubble/拓扑、精度与 Gate I 条件。
