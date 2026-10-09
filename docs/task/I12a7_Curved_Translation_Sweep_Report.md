# I12a7 / curved Q4 常速度扫掠积分报告

**PASS。** 对支持的法向常速度平移，curved native Q4 的穿面相体积已由同一旧 scalar 的完整空间扫掠测量获得。中点/四分点全为气相的六组步内事件仍得到正确非零液量。尚无域外 BC 历史或完整 driver 消费，I12 整体仍 FAIL。

## 实现

4 production + 2 test 文件，符合 [TODO](I12a7_Curved_Translation_Sweep_TODO.md)。Core strict TU 形成精确平移位置的内/外 double 包围，并区分真实 total volume 与 coordinate volume band。新的 Interface 入口核对原 interval 与 prepared region source 的双重身份，在同源曲线上积内外矩形，采用较小相稳定形成互补量。旧平面入口、库存、reflux 与 primary/characteristic 算式未改。详见 [API](../api/I12a7_Curved_Translation_Sweep_API.md)。

## 验证

- 主程序与 25 个相关测试目标 Release 构建通过，275.02 秒，无 warning。
- **24/24 CTest 通过**，156.60 秒；**22 项旧算子完整输出一致**。architecture 为 12 targets / 380 TUs / 899 files / 3085 includes / 5 legacy edges，只有 temporal 实现引入 region header 的一条内部 include 增量。
- 新合同 **66 组检查、12 类拒绝**：48 组 XY/RZ、两方向/正反相、均匀/AMR、二进制/十进制时间的独立圆/球扫掠解析量及区间切分；24 组十进制案例实际触发非零坐标包围带。两个独立已知 decimal endpoint 包围与带体积检查、8 组严格纯相 total/liquid 恒等，以及静止/纯切向零法向量通过。
- 6 组漏事件反例：XY 两方向正反速度的薄椭圆、RZ 轴向正反速度的小球，在 t/dt=1/4、1/2、3/4 上整面均为气相。完整扫掠仍得到约 `4.67095447993417e-4`（XY）和 `3.90111487804576e-9`（RZ），与独立解析交集一致。这里只初始化旧 scalar；没有任何中间时刻的解析状态重置。
- 异源/错误和陈旧 mesh、opaque velocity/coverage、非零切向运动、RZ 径向运动、出域、不可分辨内矩形、非法/过严坐标门槛及未解析闭合子单元环明确拒绝。accepted geometry identity 和原 M 不变。
- 81 冻结/参考、3 PDF、前序 I12a6 archive 不变。旧 crossing probe 仍退出 1 且完整输出不变；旧长程 driver 尚未接线，没有宣称完整 benchmark 通过。

初始漏事件 fixture 使用完全包含在单个旧叶片内的极小 XY 圆，原 Core 按既有 negative-fragment ownership 规则拒绝。该空间限制已保留为独立拒绝测试；漏时间事件测试改用跨越实际旧叶片边界的薄椭圆，使原空间拓扑已有解析能力。没有修改 Core 的空间拒绝规则，也没有把该原始失败记为成功。

日志、scope/冻结审计、前后输出逐字比较、runtime 哈希和逐项 archive 校验位于 `runs/baselines/I12a7_20261001/`；运行 `quantum_interface_curved_sweep_contract.exe` 复现完整合同，`--guards` 只跑身份/边界等检查。

## 后续

当前支持 XY 无切向常速度、RZ 轴向常速度且完整 active 原域内的扫掠。内/外液量仍是浮点值，不是严格数值下上界；两项误差门槛不改变库存 64 eps 预算。下一步补入流边界历史的显式来源身份，使整个面集合可在同一半/全区间内装配，然后实际验证完整推进和变化速度。
