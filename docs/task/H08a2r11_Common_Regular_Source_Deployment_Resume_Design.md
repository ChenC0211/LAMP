# H08a2r11恢复 / 同源实际来源部署设计

2026-10-05，实施前。仅恢复[r11 TODO](H08a2r11_Common_Regular_Source_Deployment_TODO.md)，上一轮r13为progress。r13冻结304项、14份当前交付与1383源码/输入哈希已核对。r12/r13先决均已通过；原r11失败报告保持历史。

预算仍production四文件：CchSpatialAdvector.h/.cpp、CchCharacteristic.h/.cpp；test两文件：momentum/characteristic。其余生产、phase/material/dual/KEP测试、研究driver/参考及所有原输入逐字保持。部署代码采用已由r12/r13隔离验证的原r11四份候选，再对当前r12测试仅执行必要迁移。所有候选源码仍需通过本轮标准production构建、实际生产测试和完整研究才能接受。

当前所有patch消费者重审：生产无读者，r13研究两观察器已退役，momentum仅一处共享存储身份断言，改为同一regular rows身份。实际owner删除旧patch结构、tangent SVD、faceB/cellLine/lift及其私有状态。prepare持有唯一regular准备，并由真实support knots构建pieces和primal/dual overlay；每次evaluate只apply一次。值/J、native中心查询及完整normal多项式来自同一系数，cellDivergence保留真实D_I平均量，不强制点散度常值。

运动证书采用原设计的物理坐标五点Q4完整值/缩放J迹及固定Bernstein换基，128epsilon工作量预算、无幅值floor；两时端均验。轨迹/J仍为同一Ralston3/dense/event过程。当前r12严格实际piece事件计数、常数数值J、tiny/zero/FD负控保持；仅将旧“必须出现切向跳跃”fixture改为同原输入的C1跨支持收敛/独立RK4/FD，保留完整三/四次不可见扰动负控。

顺序：标准六契约+benchmark构建；原9项短回归、原64布局/208basis/144结构/48算术矩阵及运动24时间族/288FD/2304分区；r13独立advector/closure/boundary参考。短门槛通过后，先原audit/audit32内部affine和rows/rows32最终native独立矩阵，然后标准run_case完整full 28窗口/84gate（旧4donor单列）。记录真实成本，检查每个原行、时间窗和门槛，不用observer exit0代替PDE验收。

临时证据位于 `runs/checks/H08a2r11_resume/20261005_190631_897`，正式full仍由通用启动器写入自己的时间戳目录。任何原适用门槛失败均按r11任务书保留证据并退役候选，不调容差/输入或新增隐式fallback。成功后返回H08a3真实FirstJet桥接，H11未启动。

测试迁移复核：原跳跃fixture的1e−11包含轨迹和J，不能仅换成候选中的1e−8点误差。恢复时为同一原sin输入增加独立RK4变分方程（8192/16384步、参考自差≤1e−12），实际2560步同时验point/J加参考自差≤原1e−11；原40/80/160收敛观察及独立FD保留。这是更换不再适用的分片线性精确公式，不放宽准确性。原24族时间门槛保持。
