此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_hysing_bubble/tc1/README.md。

运行入口：tools\run_case.cmd hydro_hysing_bubble/tc1 。结果：runs/hydro_hysing_bubble__tc1/<时间戳_标识>/。程序：LAMP。

# Hysing TC1：黏性主导的二维上升气泡

**状态（2026-09-29）：V9 使用质量流单元并行与分解复用继续原完整终时计算，未数值验收；V7/V8 因程序更新主动停止。** [V1](output/20260929_085129_849_V1_B4_TC1_59d62e/run.json)因圆与网格相切处的面积分类问题在初始化失败；统一几何分类后的 [V2](output/20260929_091247_567_V2_B4_TC1_init_fix_105d56/run.json)完成初始化，首个耦合步因 `volume-space refinement exhausted material resolution` 退出 3。后续诊断区分了约束求解误差与相切交点不一致，历史初态和失败信息保留。未缩短终时或改变圆形以避开问题。本目录用于集中验收气泡浮升速度、形变和相守恒；不是短程启动检查。统一问题来源、参考曲线及判据见[Hysing 算例说明](../../../../../cases/hydro_hysing_bubble/README.md)。

## 物理模型与工况

求解二维不可压缩两相 Navier–Stokes 方程，包含重力、黏性及恒定表面张力。无电场、无质量交换、无相变。计算域为 [0,1]×[0,2]；初始气泡圆心 (0.5,0.5)，半径 0.25，初速为零。正 y 方向向上，重力加速度为 (0,−0.98)。

| 参数 | 曲线内部气体（negative） | 外部液体（positive） |
| --- | ---: | ---: |
| 密度 ρ | 100 | 1000 |
| 动力黏度 μ | 1 | 10 |

表面张力 σ=24.5；目标终时 T=3。输入使用 `units code`，数值尺度与参考工况一致，不额外缩放计算域或物性。上下壁为无滑移壁，左右壁不可穿透且自由滑移；无压力出口，压力使用封闭域的体积平均参考。初态尚无求解压力，后续压力属于实际求解的中点时刻。

## 空间与时间离散

提交的基础输入采用 16×32 均匀单元，Δx=Δy=1/16；初始 Hermite 曲线 64 段，后续按统一几何投影调整材料节点。使用完整 C1 速度、相容压力、共享相质量通量、中心动量输运及同一中点的黏性、重力和毛细力。固定 Δt=0.002，共 1500 步；未启用 AMR、重距离、附加黏性或事件阈值。

这只是集中验证的基础分辨率，不预先承诺满足文献误差范围。冻结版本后与加密输入比较空间/时间敏感性；参考曲线和最终误差容差须在查看结果前固定，不按运行成败改物性。

## 运行及结果

在项目根目录运行：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/build.ps1
.\tools\run_case.cmd hydro_hysing_bubble/tc1
```

结果写入本目录的 `output/<时间戳_标识>/`，保留输入、程序副本、日志及 `output/fields`、`output/series`。场每 25 步输出，即 Δt_out=0.05；初态和末态均输出。每步历程记录时间、能量、方程残差及相库存。历史字段 `liquid_volume` 在本例指曲线内部负相的面积，即气泡面积，不能按字段旧名称误认成外部液体。

集中验收比较气泡质心高度、气相平均竖直速度、圆度、面积误差及能量收支。专用输出 `hydro_moment_phase_integrals.csv` 已接代码，选择 `negative` 行的 `centroid_y`、`mean_u_y`、`circularity`；平均速度按气泡内实际完整速度积分，不由质心差分或流场最大速度替代。`volume_fv` 与 `volume_geometry` 分别保留守恒库存和几何面积。三组[作者原始曲线](../../../tests/reference/hysing_bubble/README.md)已保存，MooNMD 的零面积列不可用于面积误差；当前输入及编译成功均不代表浮升轨迹已经通过。

## 本轮集中调试

[V7](output/20260929_112201_294_V7_shared_roots_dc8712/run.json)统一相切面区间与单元裁剪的多项式根后，首次接受首步至 t=0.002：5 次非线性迭代，原方程残差 1.42329e−11，气相面积漂移 2.77556e−17，几何—库存分数差 1.00098e−12。该版随后接受 23 步至 t=0.046，在 V8 首步接受后因程序更新主动停止，进程退出 −1，未完成原 T=3；原始日志、[观察记录](output/20260929_112201_294_V7_shared_roots_dc8712/observed_progress.json)和停止原因保留。

该运行程序为 `CB55BF5855C5E6F20A7E1B550AE52DD528F2FE5969731D04CE61F92630F2D002`，不含后来加入的可行自由最小值检查。没有把不同版本的结果合并为验收，也未完成气泡主量的参考比较。方法与后续见[流体 §71](../../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#v3-volume-and-crossing)。

[V8](output/20260929_122903_004_V8_polynomial_material_a7e5bc/run.json)使用相同字典、原网格/物性/dt/T，从初态运行。程序包含 Gram 材料装配、按多项式次数减少重复求积及非线性几何约束改进，SHA256 为 `47C5CEF5E4339EE5289F43D2227F18CE60C11EAB5BD4D085307E54E5C78838B6`；不含后来的 RZ 极点修改，本例为 XY。新旧首步动能相对差 4.39372e−16，材料装配计时及范围见[流体 §74](../../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#polynomial-material-and-rz-poles)。13:01 的观察为已接受 13 步至 t=0.026，仍运行，完整参考比较尚未完成。

[V9](output/20260929_134547_423_V9_parallel_current_9c5738/run.json)使用相同五份字典，程序加入当前完整方程的分解复用、动量积分批量计算及质量流单元并行，SHA256 为 `3BD3D8EF7E93F2668597750A112322A5FFF0097468E4B8F4EE28EAEACA1BE660`。首步与 V8 的动能相对差 1.75749e−15，两版均需 5 次非线性迭代。新首步接受后，V8 主动停止于 31 步、t=0.062、进程退出 −1；原日志与停止原因保留。13:52 的 V9 已接受 3 步至 t=0.006，继续原 T=3。不同版本轨迹不拼接为完整验收，方法见[流体 §76](../../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#parallel-mass-current)。
