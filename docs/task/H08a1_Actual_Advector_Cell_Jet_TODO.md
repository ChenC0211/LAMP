# H08a1 / 同一实际速度重构的 cell 分支 first jet

状态：**PASS，2026-10-05**，见[报告](H08a1_Actual_Advector_Cell_Jet_Report.md)。[H08a](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)的第一个实现子任务；H10o冻结基线已复核。H08a整体尚未通过，下一[H08a2](H08a2_Piecewise_Characteristic_Source_TODO.md)，H11未启动。下列范围及门槛在代码前固定，实施前副本保存在本任务checks。

目标：为唯一 `CchSpatialAdvector` 增加其**实际 cell 分支**的物理速度 Jacobian 和分支范围，使后续时空运动来源能够消费同一重构；不再另拟合一个梯度或速度。这个任务只解决局部运动学来源，不宣称已得到跨跳跃特征线。

允许修改：`include/Modules/Hydrodynamics/CchSpatialAdvector.h`、`src/Modules/Hydrodynamics/CchSpatialAdvector.cpp`（2 production）；原 `tests/contracts/hydrodynamics/cch_momentum_geometry_contract.cpp`（1 test）；API、任务设计/报告和当前入口。临时执行、比对和冻结脚本只在 `runs/checks/H08a1/<时间戳>`。不新增生产文件、CMake target 或永久工具。

禁止修改：原空间重构/atom 通量数值、material F/R、dual/KEP、Core tracer、FirstJet publication/transport、App/字典、任何库存/bubble/几何体积拟合、压力/黏性/电场。不能平滑或平均跳跃后冒充原 actual B；不能包装为已满足连续迹要求的全局 `DifferentiableVelocity`。

数学约束：材料只有用户指定的四项硬约束；不恢复 quarter 或 pressure→dual 精确交换。新导数是原局部 B 的解析微分，XY 使用物理坐标，RZ 使用 X=r²/2 后正确变换回物理 (r,z)，r=0 使用真实单侧正则极限。cell/子面交界允许单侧 jet；分支选择与既有 cellVelocity 完全相同，不能把分支导数当成跳跃的分布导数或跨事件流映射 Jacobian。

实现要求：见[实施前设计](H08a1_Actual_Advector_Cell_Jet_Design.md)。原速度值计算次序保持；返回局部物理包围盒，内部 knot 选正方向、cell 上端选最后分支。复用现有 immutable owner 和来源/clock/regrid 校验，无新插值权重。完整 H08a 后续仍须完成跨分支运动来源、真实时空相积分及 primal/dual 同一积分接线；本子任务不能替代这些门槛。

验收（实现前固定）：

- **T0**：标准 Release 构建；architecture、architecture_self_test、repository_layout 与 phase/momentum/material/dual/KEP 共8项 CTest；2/1源码预算。旧空间原契约确定性输出保持（计时除外），H10o完整研究不因添加未接入查询而重复。
- **T1**：XY/RZ × uniform/2:1 × n4/8/16 × full/selected/孔洞/不连通 × affine/非线性/任意 native mean，共144布局字段组合；每个实际 cell/piece 检验。jet.value 与原 cellVelocity 逐位相同；返回分支范围与独立 native-side 枚举相同。以实际返回值在 volume 坐标的中心差商（分片二次精确）核验全部四个物理导数，scaled error ≤2e−11；可识别解析 affine 的导数 ≤1e−12；RZ Jacobian 与 u_r/r 的真实散度匹配原 D_I，scaled error ≤2e−11。
- **T1 边界**：内部 knot 取向、cell 两侧及固定轴单侧值/导数/正则极限；共有法向迹 ≤1e−12；保留并量出切向跳跃，至少有非零反例，不要求切向连续。独立大尺度分支差商检验轴极限，容差2e−11。非有限/越cell/错pressure/FP模式拒绝；来源、clock、副本/lifetime/真实regrid沿用owner原门槛并查询新接口。拒绝不得改变既有owner。
- **T2（本子任务）**：原空间/材料/KE短回归保持；本项无时间积分或密度新分支，不将其报告为H08a的moving/crossing、R收敛或完整耦合验证。后续H08a完整门槛保留。

失败规则：原B值改变、导数不一致、分支不明确、轴或来源校验失败即FAIL，不加独立平滑速度或历史Interface fallback。仅实际外部依赖可BLOCKED。下一任务必须解决实际跨分支运动契约，而不是忽略切向跳跃。

交付：生产局部jet接口、独立短契约、来源/运行证据、API及报告；明确尚未接入FirstJet/时空F。任务拆分依据为真实owner边界及默认5/2预算，不缩小H08a总体验收。
