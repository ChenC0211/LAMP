# I12a4 / I09 区间积分消费

## 任务 ID / 标题

I12a4：I09 直接汇总同一区间的最终 actual subface integrated volumes。I12a3 经 I12a3r1 验收 PASS；依用户授权继续推进 I12。

## 目标

只解决一个事实：coarse side 的时间积分量来自同一最终子面量之和，不能重做空间/时间积分或再修正库存。

## 允许修改

最多 2 个 production 文件：HermiteGalsReflux.h/.cpp。最多 2 个测试相关文件：新 integrated reflux contract 与 tests/CMakeLists.txt。另含本 TODO、API/报告和 runs/baselines/I12a4_20261001 证据。

## 禁止修改

不改 I12a1 区间、I12a2 积分、I12a3 库存、舍入修复、geometry、driver、Hydro/Electricity、字典或容差。不将区间来源伪装成 midpoint identity，不加入 fallback。

## 数学约束

I09 持有同一 immutable integrated inventory result；同一 actual face 的最终四方向 gross volume 以及 signed volume 只汇总一次。coarse=sum(fine)，两侧 outward sign 相反。M 原样保留，结果明确关联真实 interval。side 上的 signed volume 按实际子面的缓存 signed volume 顺序求和，避免从平均 rate 再乘 dt；平均 rate 仅由 side volume/该源 dt 派生。

## 实现要求

共享原 I09 的 native incidence、2:1 pairing 和 side ordering 实现；旧 rate API 保留，避免复制两套拓扑规则。新增明确的 integrated result/side 类型，支持 XY/RZ 和已有修复后面量。过期 mesh、移走后的 handle 等明确拒绝。

## 验收

- T0：Release 主程序、目标测试、architecture、冻结文件/前序 archive 与 2+2 文件预算。
- T1：独立按 actual face owner 扫描验证 side membership、四方向 coarse=fine、同源 interval、真实 dt 平均量、source/M 零复制、internal incidence、XY/RZ、0/1/2 层 AMR、counterflow、修复后的面量。
- T1：来源生命期、错误 mesh/非法 row/moved-from 结果拒绝；coarse/fine volume residual 保持原精度。
- 直接回归：旧 I08 flux/I09 reflux、新旧 interval/planar/inventory contracts、architecture/CLI；本任务未改其它已通过算子，不重复无关测试。保留旧 crossing probe 的实际失败。

## 失败规则

配对、来源或体积汇总不一致直接失败；不重新求积、不第二次库存更新。超预算需求另拆任务。

## 交付

API/代码、独立测试、前后完整输出比较、原始日志与 SHA256 归档。通过后再推进一般时空源/完整步进等 I12 必要子任务，I12 整体通过前不进入 Hydro。

## 执行结果

**PASS。** 2 production + 2 test 文件；33 组检查、8 类拒绝与 8 项直接回归通过；6 项原完整输出一致。详见 [报告](I12a4_Integrated_Reflux_Report.md)。I12 仍 FAIL，继续 curved native Q4 与时空源/driver 的必要任务。
