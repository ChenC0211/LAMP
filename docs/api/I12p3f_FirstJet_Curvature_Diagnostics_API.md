# I12p3f / 接受态同源曲率统计

输入保持 [App](I12p3b_FirstJet_App_API.md)、[输出](I12p3c_FirstJet_Output_API.md)和[错误定位](I12p3e_FirstJet_Input_Diagnostics_API.md)不变。每个初始/接受态的 conservation.csv 在原 kappa_max 后追加四列：

| 字段 | 含义 |
| --- | --- |
| kappa_rms | sqrt(sum w_phys kappa² / sum w_phys)，物理界面加权 RMS |
| curvature_quadrature_samples | 当前 publication 中实际 quadrature 节点总数 |
| curvature_unavailable_samples | identity 不符、非有限曲率/权重或非正权重的节点数 |
| curvature_measure | 有效节点物理权重之和：XY 弧长、RZ 旋转曲面面积 |

曲率与 physicalWeight 直接读取 publication.quadrature。XY kappa 为 meridional，RZ 为 meridional+azimuthal，符合原 +div(n) 方向。归一化只除以当前界面 measure，不除以液相体积。使用 publication 原最大绝对曲率缩放平方以避免无必要的大数平方；所有统计是只读 IO 诊断，无反馈。

空界面没有 RMS，输出 NA、measure=0、counts=0；存在任何 unavailable 样本时 RMS 也为 NA，measure 只列有效样本权重。不能把缺失数据解释成零曲率或完整平均。当前 publication 在发布时已拒绝非法 sample/weight，因此正常接受态 unavailable 明确为 0；这是已验证数据的计数，非未实现字段的占位零。

原 kappa_max 包含 publication 已检查的 graph endpoints 和 quadrature nodes；RMS 只使用上述物理求积，不把端点重复计权。VTP 的九点显示 polyline 不参与统计。即使 writeFields=false 或空间帧稀疏，全部接受态都计算/记录；耗时归入原 IO stage。

实现及验收范围见 [TODO](../task/I12p3f_FirstJet_Curvature_Diagnostics_TODO.md)。当前六个 App 问题都有界面；空界面语义由上述 IO 分支定义，本项没有声称新增空界面 App 基准。
