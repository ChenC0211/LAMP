# I12a11r4 / cut 求根精度与纯相舍入预算分离报告

**FAIL：原 17 eps 纯相拒绝消除，但原四步配置在第二步 midpoint 闭合出现可达区间拒绝。** 完成 [TODO](I12a11r4_Pure_Phase_Roundoff_Budget_TODO.md) 的实现；未以子问题修复代替整项验收，I12 仍 FAIL。

## 实现和单位验证

2 production + 2 tests。cut-cell 的初始跳过、括区、求根与最终工作精度保持 r2 的 `min(用户上限,16 eps)`；完整证书确认的纯相改为 `min(用户上限,64 eps)`。Report 分别列出两个实际预算。原库存 64 eps、M、first jets/B/J、位移/拓扑/求积/迭代 gate 未变，未恢复旧 1e-10 级尾量。

首轮 Release 两目标构建 9.71 秒；bubble 单位验证 71.98 秒 PASS。原制造解、定向求根、纯相 4 eps、拓扑及 SAFEGUARD 共 23 条数值行与 r2 完全一致，Newton/bisection/endpoint 仍为 503/128/0。新 4 组 XY/RZ 正反相中，32 eps 尾量原样保留且不改系数，128 eps、原 1e-12、宏观尾量均拒绝；将用户上限设为 8 eps 时也拒绝 32 eps。候选 identity/M 不变。

## 原推进结果

同一 n=32 XY AMR、dt=.01、非多项式、四步请求，开启原库存越界审计与 r3 的只读纯相分解。191.03 秒后退出 2。

首步 cell 1513 的 17 eps 尾差与 r3 逐项相同，未被裁剪或覆盖；这次成功接受 t=.01。下一半步仍保留相同尾差，无净面交换、无误差增长。但在另一 cut cell 1454 发生：

```text
midpoint_closure
[Interface/HermiteGalsBubble primary=13 cell=1454 stage=volume bracket]
target outside certified attainable volume interval: directed search cannot extend bracket
```

失败前已接受 t=.01 的位置、梯度、曲率误差分别为 1.8390815e-10、3.3122612e-8、3.2320791e-3；初始曲率误差为 4.7399677e-6。严格体积求根没有保证曲率改善，不能将已接受一步视作精度通过。实际 root 区间/目标差额及全局与局部 Green 体积误差仍需下一项诊断，当前不提前归因。

没有库存越界诊断行；失败仍保持 accepted identity/M。由于原先决四步未完成，本轮 40 步多项式和 12 组非多项式矩阵记为 NOT_RUN。

## 完整回归与归档

完整 Release 重建 31.83 秒，无 warning；30/30 CTest 通过，203.97 秒。29 项原输出逐字一致；bubble 只增加上述纯相预算验证，原 23 条数值行保持。architecture 为 12 targets/380 translation units/900 files/3088 includes/5 legacy edges，与 r2 不变。旧 midpoint probe 退出 1、输出不变；81 冻结/参考文件、3 PDF、前序 archive 哈希保持。

源码、28 可执行及运行库、原始日志、逐项比较和哈希归档于 `runs/baselines/I12a11r4_20261001/`。数值 FAIL 与归档完整性 PASS 分开记录。下一项只读比较 cell 1454 的同一 scalar 在现有全局 Green 与已有局部 Green 体积测量下的目标/区间差额，并检查首步 bubble 的曲率影响，随后依据证据拆修复。
