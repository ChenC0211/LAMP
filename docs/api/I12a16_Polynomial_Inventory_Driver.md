# I12a16 / Polynomial 库存驱动说明

见 [TODO](../task/I12a16_Polynomial_Inventory_Driver_TODO.md) 与 [报告](../task/I12a16_Polynomial_Inventory_Driver_Report.md)。这是原 benchmark 的显式新路线，尚未通过完整驱动验收，不能作为 production 接口。

命令：

```text
quantum_interface_hermite_gals_benchmark.exe --polynomial-inventory-smoke
quantum_interface_hermite_gals_benchmark.exe --polynomial-inventory-full
```

smoke 对原 rotation/vortex、uniform/AMR 各推进一步；full 使用原 rotation 30 步、vortex 400 步。均 n=8，circle cx=.65、cy=.53、radius=.13；rotation speed=2、dt=.01，vortex speed=40、period=.4、dt=.001。失败从真实 producer/assembly/closure 抛出，保留旧 accepted 状态，并使新 ordinary CTest `interface_polynomial_inventory_driver_smoke` 失败，不设置 WILL_FAIL。

`Case::polynomialIntervalInventory` 默认 false；新命令同时启用原 boundedInventory 与 meshScaledGeometryBudget。原 translation/旧 rotation/vortex 命令保留原路线，没有失败 fallback。owned polynomial 使用原速度的真实两端系数：rotation degree1，vortex degree3；cosine 时间由原两端数据线性重构，其误差必须留给后续完整时间收敛测试，不能称为连续 cosine 真解。

end 与 midpoint prediction 都使用旧 accepted scalar；midpoint 的实际身份交给同一个 provider 的 `polynomialArrivalView`。vertices/edges/centers 仍走原 actual RK2/first-jet 路径，原显式边界 first-jet 条件保留。半步/整步 integrated flux 均保持同一个旧 snapshot、velocity owner、characteristic owner 与旧库存区域。

每个 native face 的 absolute volume tolerance 固定为关联最小 cell capacity 的 1e-6，max panels=4096，depth=40。这是积分误差预算，与几何 τ=max(1e-10,(hmin/L)^3)、Hessian h² 及 M 的64eps roundoff修复限制相互独立。代码完整保留 transfer→shared inventory assembly→reflux→I10/I11 两阶段发布流程。

输出 `POLYNOMIAL_PREDICTION` 记录实际两阶段预测完成；`POLYNOMIAL_INTERVAL` 仅在全部 face 和库存装配成功后记录；`POLYNOMIAL_FACE_FAILURE` 记录真实面、邻接cell、区间、误差预算及拒绝。没有边界相历史时不会按解析圆位于域内就给面量填零，亦不借用只适用于恒速的旧 inflow capability。

本项现有失败来自物理入流面所需的相历史未被新 polynomial integral 表达，下一任务须以显式来源与可验证边界覆盖补齐。未通过前不能声称完整步、时间阶、动态 AMR 或 I12 通过。
