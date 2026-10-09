# I12a9 / 常速度半步与整步来源消费

`HermiteGalsTemporalVelocitySource::constantArrivalView(actualCandidate, actualArrival)` 返回实际 candidate/time 的 velocity view，保留同一 immutable constant provider。必须是 owned constant，begin 等于原 previous 的真实时间、step 相同、声明 end 不超过原 end，且实际 arrival 位于原正向时间范围内。opaque 或 moved-from source 拒绝。派生 view 在原句柄析构后仍保持 provider；固定 RK2/B/J 没有变化。

benchmark 的 `--translation-sweep-smoke` 与 `--translation-sweep` 是明确的新驱动入口：整步 prediction 创建 source，半步从该 source 派生真实中点 view。所有 interior first jets 和 center 使用原旧 scalar 回溯，外边界 first jet 使用原 prescribed BC。完整常气相入流历史由制造解的整个入流平面与整段时间范围证明，不对 interior 重置解析值。

旧 region source 只准备一次。half/full 实际面集合分别经 I12a8 区间积分、I12a3 独立库存装配、I09 汇总；随后各自经过原 I10/I11。新驱动逐项检查共享旧 snapshot、速度/特征/BC owner，以及发布后 M 与最终共享通量形成的 M 完全相同。失败保持 accepted identity/M。

没有 quarter geometry 或 midpoint 通量求积。发布的 midpoint 仍是实际中间状态；完整非恒定速度驱动尚待实现。旧 driver 保留为独立历史对照，不作为新入口的 fallback。验收范围与状态见 [TODO](../task/I12a9_Translation_Driver_TODO.md) / [报告](../task/I12a9_Translation_Driver_Report.md)。
