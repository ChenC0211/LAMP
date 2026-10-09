# I12a9r1 / 横向根归属分辨报告

**算子检查 PASS；完整平移验收尚未完成，本项暂记 FAIL。** 原 XY uniform 的第五步根包围误拒绝已消除，实际接受全部 40 步至 t≈.4。随后误差参考计算触发独立解析圆面积求和错误，脚本退出 2，后续三组平移未运行，不能据此宣称完整 I12a9/I12 通过。

## 最小修复

符合 [TODO](I12a9r1_Root_Enclosure_Ownership_TODO.md) 的 1 production + 2 test 文件预算。唯一生产变更在 strict Core 的 MappedInterface 实现：未分辨的横向根代表点或包围端点归属触发最多 4 次已有精确 root-refinement，每次确认同一 simple-root 导数方向。原唯一 owner 检查仍为最终门槛，真实重叠、缺区域继续拒绝。不分配边界 owner、不丢弃微小弧；不改 scalar、外部共享 face-root、GL8 公式、B/J、bubble 或 M。

独立复现显示，仅允许包围端点接触共享边界不足以解决问题，因为粗根代表点本身也可能落在边界。本项在生产修改前更新 TODO，采用同根细化并保留原判据。

## 独立测试与记录

浅抛物线 `v=.625+(q-.5)^2/8-depth` 有独立可积的小帽。32 组 XY/RZ、两个 graph 方向、正反相及 depth=2^-24/2^-42/2^-48/2^-52 通过；最小帽宽约 8.42937e-8，Green 积分量级 1e-23，微小弧保留。真实重叠、缺区域、缺实际裁剪 cut、异源 snapshot 共 4 类拒绝通过。

原生产路径上的未分辨代表点误拒绝见 `reproducer_before_ownership.log`。最初纯相对微小积分断言比原 binary64 横向根分辨率更强，其失败日志保留；最终检查使用明确记录的浮点根尺度加积分误差尺度，另要求非零弧和非零积分。外部 clipping roots 则由调用者通过原 immutable refined API 预先细化，以独立核验微小弧宽。该测试没有宣称任意趋零相体积的固定相对精度或严格误差界。

一次新增 guard 的编译因误用不可复制的 LevelSetMappedSegment 失败；其后误运行旧可执行的初步结果已明确标为 invalid。修正为独立 fixture/move 后重新编译，新可执行的 32/4 检查通过；旧错误日志保留，不作为验收结果。

## 长程运行的新发现

原 `--translation-sweep` 命令和 driver 完全未改，用时 101.92 秒。XY uniform 已接受 40 步，随后 `errors()` 中的解析 cell-volume 求和检查失败。原 harness 未在后处理前更新 stage，错误行仍显示 `end_closure`；catch 再次调用同一参考函数又抛出，使进程中止。故不能把该行解释为 I10/I11 闭合拒绝，也没有剩余三组本轮长程证据。

只读端点诊断复现了参考原函数的消减误差：center=.51、radius=.19 时，`asin(d/r)` 形式在右端点相对四分之一圆 primitive 的误差约 -3.24223e-11；使用相同 d/s 的 atan2 形式误差约 3.47e-18。这个诊断只涉及独立误差参考，不向数值状态写入值。下一子项修复参考计算与失败标签，再执行所有平移和后续 I12 验收。

主程序与 27 个测试目标 Release 构建通过，33.15 秒，无 warning；**27/27 CTest 通过**，326.41 秒，**26 项旧测试输出逐字不变**。architecture 保持 12 targets / 380 TUs / 899 files / 3085 includes / 5 legacy edges。81 冻结/参考、3 PDF、I12a9 archive 哈希核对通过；旧 crossing 退出 1 且输出不变。源文件预算、runtime、逐项 archive 哈希审计与全部原始证据位于 `runs/baselines/I12a9r1_20261001/`。证据完整性 PASS 不等于尚未完成的长程验收 PASS；I12 与 Gate I 继续未通过。
