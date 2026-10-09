# I12p2a / 纯路线收敛与跨面结果

**PASS（本任务）；完整 Gate I 仍未通过。** 仅修改 benchmark 和测试注册，没有修改 production 算式、容限或历史算例。

12 组静态非多项式 curvature、12 组非多项式平移空间细化、6 组旋转时间细化、48 组非对齐步内跨面全部通过。h 取原 n=8/16/32；平移使用原 dt=.01、4 步；旋转 t=.3，dt=.01/.005/.0025。独立审计从原始误差重新计算阶数，与程序输出一致。

| 配置 | 静态梯度阶 | 静态曲率阶 | 平移梯度阶 | 平移曲率阶 |
| --- | ---: | ---: | ---: | ---: |
| XY uniform | 3.958 | 3.304 | 3.520 | 2.817 |
| XY AMR | 4.095 | 3.413 | 4.283 | 3.191 |
| RZ uniform | 4.015 | 2.970 | 3.503 | 2.442 |
| RZ AMR | 4.082 | 2.970 | 3.503 | 2.460 |

表中为最后两档的观测阶，不据此声称普遍四阶或三阶。所有档位梯度/曲率误差递减，满足原末档 >=2 / >=1.5 的门槛。旋转位置的末档时间阶为 uniform **2.00417**、AMR **2.00468**，满足原 >=1.8 要求。

48 组 crossing 覆盖 XY x/y、RZ z、两种网格、正负速度/相位、1/4 和 3/4 步内穿越原 .5 面；midpoint/end 均使用实际旧 scalar 输运和纯 publication。最大逐 cell 几何体积误差 4.441e-16，first-jet/normal 误差 6.107e-16，绝对曲率 5.742e-14，满足原阈值。没有 M、守恒相通量或 ledger 声称。

原 28 CTest 全部通过（180.47 s），28 项完整输出逐字保持；旧 legacy crossing 拒绝保持。4 个新研究分别直接执行通过（4.00 / 19.65 / 100.94 / 5.34 s），并核验其 CTest 注册；没有重复运行新研究来扩大通过数量。仅 benchmark 可执行变化，其余 38 个原 runtime 哈希保持；无 warning，81 freeze、3 PDF 和 p1 archive 保持。

本项只证明上述固定网格/静态 AMR 及原输入下的收敛和跨面精度。动态 refinement/coarsening 的状态/几何传递尚未验证，应进入 I12p2b；不能以此宣告完整 I12 或低黏流体耦合通过。

证据：[汇总](../../runs/baselines/I12p2a_20261002/verification.json)、[完整误差/阶数及跨面数据](../../runs/baselines/I12p2a_20261002/convergence_verification.json)、[范围审计](../../runs/baselines/I12p2a_20261002/source_audit.json)、[旧结果比较](../../runs/baselines/I12p2a_20261002/regression_comparison.json)。
