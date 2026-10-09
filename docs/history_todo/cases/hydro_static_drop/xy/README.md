此页为 G02 整理时保留的历史算例说明。数值、限制和既有结果供追溯；其中旧脚本、旧程序名和输出位置不再是当前操作入口。当前说明见 cases/hydro_static_drop/xy/README.md。

运行入口：tools\run_case.cmd hydro_static_drop/xy 。结果：runs/hydro_static_drop__xy/<时间戳_标识>/。程序：LAMP。

# B2-XY：平面静滴

**状态（2026-09-29）：V8 已完成 3 步，因大网格内存调度主动停止，完整静滴未验收。** [V1](output/20260929_085129_833_V1_B2_XY_d4712e/run.json)在初始化时内存分配失败；保持全域矩阵稀疏、压缩边界支撑后的 [V2](output/20260929_091237_332_V2_B2_XY_init_fix_944840/run.json)已经保存初态，后续首步工作集增长到约 20.9 GB 时主动停止，详见[资源记录](output/20260929_091237_332_V2_B2_XY_init_fix_944840/resource_stop.json)。此前将 V2 称为“尚未完成初始化”不准确，已由实际 t=0 输出纠正；原失败证据保留。后续 V7/V8 状态见文末。

两相不可压缩 NS、黏性及数值曲率毛细力；圆心 (0,0)、R=0.4，域 [−1,1]²，四边自由滑移。两相 ρ=1、μ=0.0816496580927726，σ=1，La=120；初速零，压力由方程建立。平衡压差为 2.5。

64×64 均匀网格、64 个初始弧段，dt=0.004898979485566356，T=7.83836717690617，共 1600 步，每 40 步写场。完整模型、参考、误差定义与加密比较见[统一说明](../../../../../cases/hydro_static_drop/README.md)。

输入：[公式](../../../../../cases/hydro_static_drop/xy/system/geometry.dict)、[物性与边界](../../../../../cases/hydro_static_drop/xy/system/hydro.dict)、[界面](../../../../../cases/hydro_static_drop/xy/system/interface.dict)、[网格](../../../../../cases/hydro_static_drop/xy/system/mesh.dict)、[运行](../../../../../cases/hydro_static_drop/xy/system/run.dict)。

```powershell
.\tools\run_case.cmd hydro_static_drop/xy
```

现行初值为隐式 x²+y²−R²=0，种子 (0.4,0)、traceTolerance=1e−10；提取后仍使用同一 Hermite 曲线。结果在本目录 output/<时间戳>/，场与历史位于运行副本的 output/fields/ 和 output/series/。[隐式输入预检](output/20260928_233051_762_implicit_input_5d7b3c/run.log)和[原参数输入预检](output/20260928_230935_411_analytic_input_bdf638/run.log)保留追溯。

**集中调试状态（2026-09-29）：** [V6 原工况](output/20260929_094915_349_V6_material_profile_b09596/run.json)已完成初始化并保存 t=0 场与相库存。因首个耦合步成本高，在 512.69 秒主动结束诊断，接受步仍为 0，工作集峰值 16.043 GB；[停止说明](output/20260929_094915_349_V6_material_profile_b09596/diagnostic_stop.json)区分诊断中止与数值失败。没有完成静滴平衡、寄生流衰减或网格收敛验收。此前 V2/V4 同样已有初态，不能仅凭未及时刷新的控制台日志判为仍在初始化；详见[流体 §69](../../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#v1-local-coordinates)。

**V7 更新：** [原输入运行](output/20260929_100803_035_V7_block_solve_49ad23/run.json)使用同一联立系统的精确块消元，在 8 次非线性迭代后接受首步至 t=0.004898979485566356。方程残差 4.2354e−11、散度残差 9.3974e−16、库存漂移零、几何—库存分数差 1.0032e−12。第 2 步进行中主动停止，退出 −1，见[接受态及停止记录](output/20260929_100803_035_V7_block_solve_49ad23/diagnostic_stop.json)。总诊断时间 1572.52 秒包含初始化、首步和部分次步，峰值工作集 14.291 GB；重复分解和装配仍昂贵。大场输出间隔为 40 步，本次只有初态大场和首步历程，未达到原 1600 步终时，不作为静滴数值验收通过。方法与后续见[流体 §70](../../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#v2-exact-block-solve)。

**V8 更新：** [新程序原工况](output/20260929_131632_829_V8_reused_block_76ae37/run.json)包含分片多项式材料积分及当前鞍点方程的分解复用，五份字典与 V7 哈希一致，原终时保持。程序 SHA256 为 `B8D956D3EAB7756206CBA92B4A951FA97266809560B801A9E190B5F6297266C4`。首次分解约 38.09 秒；首步内后续轮次已复用该分解，仍按当前完整方程与真实非线性残差判断接受。13:26 的观察为第 5 轮求解进行中，尚未接受首步，不将线性求解通过视为物理步完成。方法和动态观察见[流体 §75](../../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#retained-block-factor)及[本包证据](../../module_history/output/20260929_130951_reused_block_factor/README.md)。

**13:37 后续观察：** V8 已接受首步至 t=0.004898979485566356，8 次非线性迭代仅做一次速度/压力分解；方程残差 4.05386e−11，液量漂移 1.11022e−16，动能 6.15992e−19。随后继续第二步。该冻结程序不含[下一包单元并行](../../../docs/QUANTUM_CCH_Hydrodynamics_TODO.md#parallel-mass-current)，不移记新程序检查结果，也未完成原 1600 步验收。

**14:22 调度停止：** V8 最后接受 3 步至 t=0.014696938456699067，最后一步方程残差 4.12773e−11，液量漂移 −1.11022e−16。为让原 64×64 B3 无黏振荡获得大网格内存名额而主动停止，程序退出 −1、启动脚本退出 1；[停止原因](output/20260929_131632_829_V8_reused_block_76ae37/stop_request.json)及[运行元数据](output/20260929_131632_829_V8_reused_block_76ae37/run.json)明确记录 `stopped_for_resource_scheduling`。这不是数值失败，也不是完整静滴通过；原终时、网格和判据保持，后续需收齐完整轨迹。
