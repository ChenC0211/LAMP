# I12a5 / 同源 Q4 局部 Green 势

## 目标与范围

I12a4 已 PASS。为后续薄扫掠区域的液量测量，单独验证一个数学事实：平移 Green 势的原点不改变散度与液体区域体积，也不改变原 scalar、根身份和区域归属。

允许 3 个 production 文件：Core 的 LevelSetMappedInterface.h、Detail/LevelSetMappedInterface.h 与 LevelSetHermiteGeometry.cpp；最多 2 个测试相关文件：新增 shifted Green contract、tests/CMakeLists.txt。另含本 TODO、报告/API 与 runs/baselines/I12a5_20261001 证据。先写 TODO 再改代码。

## 数学与接口约束

新增显式带物理原点 `(x0,y0)` 的 quartic 积分重载。XY 势为 `(x-x0)`、`(y-y0)`；RZ 势为 `pi*(r-r0)*(r+r0)`、`2*pi*r*(z-z0)`。差值在原 strict Core 的 retained rational chart 内先形成，再转浮点，不能先舍入全局根坐标后相减。原点仅为积分势常数，可在域外，不是 RZ 物理位置。

旧入口和零原点保留原浮点算术顺序。原 scalar、graph、actual cuts、region owner、全部 sensitivity 列、根细化请求及有界 GL8 控制不变；误差仍只是浮点估计。nativeMoments 的分部积分公式不在本次修改范围。调用方必须给 fixed boundary 使用同一原点；独立 free moment 本身会随原点变化。

## 禁止

不实现 rectangle volume、时间通量、库存或 driver；不改特征线、容差、冻结文件、原测试预期、Hydro/Electricity。非有限原点明确拒绝。不得把原点平移当作新的几何来源或误差认证。

## 验收

- 独立多项式原函数 oracle：XY/RZ、两种 graph 坐标、两种相方向、非零/域外原点、受控及原 GL8 路径。
- 同一曲线分区加和、零原点与旧入口逐字段一致、sensitivity 不受常数势影响；薄区域/大平移坐标检查避免先舍入根再相减。
- 非有限原点、错误源 cut 等原合同拒绝；失败不改变来源。
- Release 主程序和直接相关 Core/I07–I11、I12a1–a4 检查，architecture/CLI；前序完整输出比较、81 冻结/参考与 3 PDF、文件预算、前序 archive 完整性。
- 报告保留当前 I12 FAIL；本项通过后继续子区域测量和时空积分。若新数学障碍超出本项，则归档失败并另写任务。

## 执行结果

**PASS。** 3 production + 2 test 文件；145 组检查、8 类拒绝及 22 项 CTest 通过；21 项旧输出一致。见 [报告](I12a5_Shifted_Green_Moments_Report.md)。
