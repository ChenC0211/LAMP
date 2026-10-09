# H08a2r11 / 实际同源部署前设计

2026-10-05；依据[任务书](H08a2r11_Common_Regular_Source_Deployment_TODO.md)。本轮仅此任务。材料仍仅四硬约束，不新增点态散度、dual库存配分或R=0要求。

改动预算：production 4（CchSpatialAdvector.h/.cpp、CchCharacteristic.h/.cpp），test 2（原momentum、characteristic）；不使用第五文件。

当前消费者：CchPhaseGeometry读取canonical atom、primal/dual视图及准备统计；材料/dual使用完整normal polynomial积分；最终native对流使用同一owner的native查询；Characteristic索引pieces、读取pieceVelocityJet及两端full clock。旧patches仅momentum测试读取，无生产消费者。删除旧patch结构/接口、SVD tangent拟合、faceB/cellLine/lift；准备统计改为真实regular权重与support pieces，tangentUnidentified遗留字段只报告低阶不可辨识行，不作为测试适用范围选择。

实施后纠正：上段消费者清单漏列了既有benchmark driver的两处旧patch观察器；它们没有生产数值用途，但删除API后导致研究构建失败，不能忽略。候选已退役，见[报告](H08a2r11_Common_Regular_Source_Deployment_Report.md)；先完成r12/r13再恢复部署。此注保留原设计和实际偏差，不追溯伪称实施前已经发现。

每次evaluate只做一次regular.apply，并持有同一系数、native值、边界、几何身份和完整时钟。cell值/J与normal polynomial全部调用此系数；cellDivergence直接保存原DI平均量。nativeVelocity定义为实际face中心的同一连续物理速度（不把面均值冒充点值）。piece闭包保留指定piece边界，C1保证支持边界两侧值/J相同。

准备按regular cells的实际supportKnots划分pieces；primal/dual线overlay同时加入其来源cell支持分割点，primal atom保留原来源face标识，由固定相邻cell取得完整法向四次式。源origin共享和独立width双chart保持。

运动证书改为物理坐标五节点完整Q4迹：两个速度与四个物理J分量（J乘物理cell宽度），不使用r²坐标或r加权以免提高次数。五节点到Bernstein系数是固定线性变换，舍入work累计变换绝对项；同一128 epsilon相对工作量门槛，不设绝对速度floor。现有Ralston3/dense/事件规则不改。Detail证书函数置于既有Characteristic头/源供独立负控；负控含旧三采样节点不可见的三/四次差、真实常数/导数差以及零/极小幅值。轴正则由原regular和消费者门槛验证。

原跳跃fixture改验实际C1、跨support轨迹与独立同场积分/FD，不能继续要求新场产生旧模型跳跃。原24时间族、288 FD点和2304分区比较保持，事件拒绝规则保持。先短契约，再完整native/affine及H10b-v2原28窗口/84gate；任何原适用精度失败不调整门槛。失败候选按任务书保全证据后退役；不开始H11。
