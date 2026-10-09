# I12a11r4 / cut 求根精度与纯相舍入预算分离

## 目标

针对 [I12a11r3](I12a11r3_Pure_Phase_Gap_Audit_Report.md) 已确认的 17 eps 纯相尾差，修正 r2 将 16 eps 求根工作精度施加给全部单元的错误。cut scalar 仍解同一方程至 `min(用户上限,64 eps/4)`；纯相不具备改变液量而保持拓扑的自由度，其允许残差取 `min(用户上限,64 eps)`。不提高原库存预算，不恢复旧 1e-10 级尾差，不修改 M 或通量。

## 范围和实现约束

2 production：HermiteGalsBubble.h/.cpp。2 tests：bubble contract 与既有 benchmark。另含本 TODO/报告、bubble API/README 与 `runs/baselines/I12a11r4_20261001/`。

只有初始 pure guard 与最终 pure 核验使用纯相预算；所有 cut 的初始跳过、区间、求根终止、最终核验保持 r2 工作精度。Report 明确记录两个实际预算，其字段不参与算法分支。原 signed-partial/拓扑/位移/64 迭代/moment gates、B/J/first jets 保持。

## 验证

T0：Release 无 warning；完整 30 CTest。对比 r2/r3 各项输出，cut 求根迭代/数值应不变，只允许新纯相预算验证及完成此前被拒推进带来的差异。81 冻结/参考文件、3 PDF 和前序 archive 不变。

T1：原纯相 4 eps 小尾差仍通过，新增 32 eps（工作精度外、库存预算内）XY/RZ 正反相对照，原样保留 M/系数与残差；128 eps、原 1e-12、宏观 mismatch 均拒绝。用户更严格上限对纯相仍生效。cut 制造解的最终残差继续不超过其工作精度。

T2：r3 只读分解接线改为要求原 n=32 四步真实完成，保留首步 cell 1513 的 17 eps 证据，同时开启旧库存越界诊断并要求零越界/零拒绝。只有此项通过后再运行四组 40 步多项式平移/四组静止与原 12 组非多项式矩阵。曲率递减单独记录，不能用库存成功掩盖它。

## 失败与交付

原拒绝继续出现或产生新停滞/超预算差额则本项 FAIL，保留原日志并拆最小根因任务，不通过减小案例范围、放宽 64 eps 或重置几何实现通过。源码/runtime/日志/hash 归档。即使本项通过，旋转、vortex、时间/空间与其余 I12 gate 未完成时仍继续推进。
