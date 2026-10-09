# I12a11r9 / 平移曲率逐阶段审计报告

**PASS（只读诊断）；I12 数值验收仍 FAIL。** 对应 [TODO](I12a11r9_Translation_Curvature_Stage_Audit_TODO.md)。0 production + 1 test，仅扩展 benchmark 的显式诊断入口。

## 实现与不变性

`--translation-curvature-stage-audit` 使用原 XY uniform n=32、dt=.01、四步平移 driver，记录 4 个 accepted 起点、8 个 predictor、8 个实际 bubble-closed scalar，以及 8 份求根统计。独立二阶解析导数与 cardinal bubble value/gradient/Hessian 对照继续使用；最大 jet 缺陷为 1.1143e-13，低于原 1e-10 检查。private scalar 只读、不发布、不反馈推进；candidate identity、M 和未绑定状态均保持。旧工作精度入口仍默认只测首步并执行原双体积比较。

原四步最终 BENCHMARK 行与 r8 矩阵逐字一致：位置 7.8330e-8、梯度 2.6939e-5、曲率 1.0386e-2、体积 9.0911e-12、ledger 1.3878e-17、reflux 0、逐单元 M 误差 2.3410e-8。独立审计与 benchmark 的末步曲率误差因等价公式舍入相差约 8.9e-16，不将这两种公式的结果宣称为逐位相同。

## 阶段证据

下表均为各阶段实际界面求积节点上的最大曲率误差；不同阶段的最大值可能出现在不同位置。

| 步 / 阶段 | 起点 accepted | predictor | closed | closed / predictor | closed 最差 cell |
|---|---:|---:|---:|---:|---:|
| 1 / midpoint | 9.7002e-6 | 8.0572e-6 | 5.4950e-5 | 6.8200 | 370 |
| 1 / end | 9.7002e-6 | 8.7629e-6 | 5.2636e-5 | 6.0066 | 370 |
| 2 / midpoint | 5.2636e-5 | 1.2739e-4 | 4.8911e-4 | 3.8393 | 617 |
| 2 / end | 5.2636e-5 | 2.0692e-4 | 7.3548e-4 | 3.5544 | 617 |
| 3 / midpoint | 7.3548e-4 | 5.2019e-3 | 2.5705e-3 | 0.4941 | 682 |
| 3 / end | 7.3548e-4 | 1.4937e-2 | 4.9324e-3 | 0.3302 | 650 |
| 4 / midpoint | 4.9324e-3 | 8.3431e-3 | 7.8691e-3 | 0.9432 | 650 |
| 4 / end | 4.9324e-3 | 1.3911e-2 | 1.0386e-2 | 0.7466 | 650 |

前两步 closure 放大误差；第三步 predictor 已先出现较大误差，随后 closure 反而降低最大值。此证据要求追踪 accepted 场到下一步 predictor 的传播与重构，不能仅凭最终误差归因于 bubble 或 compact mixed closure。

## 求根证据与界限

8 次闭合累计 410 次 Newton、0 次二分、0 次端点跳转。所有 closed 最差单元均为实际激活的 cut cell，最终系数均不在证书端点；本案例不支持“端点选择导致误差”的假设。最大整体闭合相对残差为 3.4417e-15，原 active 16 eps / inactive 64 eps 预算不变。

第二步 midpoint 的 cell 617 容量 9.765625e-4，初始相对体积残差 4.6455e-12，系数修正 2.1122e-6，最终体积残差为零；最终 dV/dc 为 −2.1478341e-9，实际割线斜率为 −2.1478337e-9。该批 48 个激活单元、48 次 Newton、最大 1 次迭代。第二步 end 同一单元初始相对体积残差 1.8170e-10，系数修正 3.7404e-7，最终残差为零，dV/dc 为 −4.7439524e-7，割线斜率为 −4.7439489e-7；该批 49 次 Newton、最大 2 次迭代。

第二步 end 的全局最大系数修正为 6.9023e-6，显著大于当时最差曲率单元的修正。本轮尚未记录该全局最大修正的 cell 身份，也未区分旧场精确平移与新 Q4 重构之间的 Hessian 误差。系数 / h² 仅为量纲指标，不能代替曲率界。

## 验证与交付

Release 构建 3.68 秒，无 warning；四步诊断 146.59 秒、退出 0，4 个起点 / 20 个几何阶段 / 8 次闭合 / 8 份 root stats 完整，拒绝与诊断异常均为零。11/11 直接回归通过，47.33 秒，各项输出逐字等同 r8；旧 midpoint crossing probe 保持预期退出 1 和原输出。其它 27 个可执行、全部 production、81 项冻结/参考文件、3 PDF 和 r8 archive 哈希均不变。

原始日志、字段表、20 阶段及 8 份 root 数据与比值、回归对照、源码/runtime 哈希、交付清单及 archive 位于 `runs/baselines/I12a11r9_20261001/`。下一任务只读比较 actual accepted scalar 的精确常平移拉回与下一步 predictor，定位重构误差，并追踪大系数修正的单元身份。I12/Gate I 仍 FAIL，不进入 Hydro。
