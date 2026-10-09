# I12a10r2 / I06 声明外边界对称性闭合

## 目标

只解决一个事实：已严格满足 symmetry 法向 first-jet 条件的同一 primary，在 I06 mixed/AMR 闭合后仍保持整条声明外边界法向导数恒零。依据 I12a10r1 的 12 个可复现拒绝与 private 反事实证据。

## 允许修改

最多 2 production：HermiteGalsAmrClosure.h/.cpp；最多 2 test：既有 AMR contract 与 benchmark。不增加新 CMake 测试。另含本 TODO/报告、I06 API 补充、README 与 `runs/baselines/I12a10r2_20261001/` 证据。

## 禁止修改与数学约束

保持 B/J、first jets、M、I07–I11 标量/体积/拓扑/位移 gate、Core exact restriction、所有 tolerance 和冻结资料。沿既有 RZ 反射边界路线：显式声明的真实外边界上 normal first jets 必须严格零，派生 mixed 是边界约束未知量；先装配受约束兼容系统，再生成 exact hanging skeleton。不在解后按阈值抹零，不自动从小数值猜测 BC。

## 实现要求

- I06 接受与 I07/I10/I11 相同的 SymmetryPlane 列表，默认空，原调用结果保持。
- 检查 plane 的坐标/方向/真实外边界/重复/与自动 RZ 轴身份冲突/active trace。严格拒绝非零 vertex 或 master edge 法向值，包括极小非零。
- 对反射法向的紧凑导数使用原 axis 零约束；mixed 的系统行施加 q=0，精确约束后检查实际边界节点。无声明部分仍原五点、平均、线性系统及双残差检查。
- benchmark 初始化与完成半/整步均传同一声明。保留无声明原故障诊断，不能改为测试期望通过而掩盖原反例。

## 验收

- T0：Release 无 warning，architecture、冻结 81 文件/3 PDF/前序 archive 不变；说明新增 include 边。
- T1：XY 四外边界、均匀/AMR、非平面 Q4 的 mixed 再现与 C1；RZ 自动轴加外 symmetry，严格边界 Q4 法向 DOF；非法声明、非零 first jets、无 active trace 拒绝且 candidate identity/M 不变。
- T2：原 48 个平面 crossing publication 全通过真实 half/full 链；9r2 四组 40 步平移与四组静止不退化；原 CTest 全套及旧 midpoint crossing 结果按预期保持。数值差异逐项分析，不以总 exit code 代替。

## 失败规则与交付

任何未解释回归、精确边界非零或完整 crossing 失败均记录 FAIL 并保留证据，再拆下一项；不得提高 gate。交付源码/API/报告/原始日志/运行时及逐项哈希 archive。即使本项 PASS，非恒定速度与其它 I12 条件未通过时父任务仍 FAIL，继续推进，不接 Hydro。
