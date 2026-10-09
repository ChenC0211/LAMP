# I12p3c / FirstJet 正式输出

**G04 现行输入/来源布局（2026-10-03）：** 使用 [五文件 schema](G04_Interface_Dictionary_Schema.md)，五个源文件均记录哈希；effective_config 按 owner 汇总，并在 meta/effective 导出五份可回读字典。下述 G03 单文件输入说明已被替代，其 NA/null 诊断语义保留。

**G03 现行差异（2026-10-03）：** 输入改用 [通用数据配置](G03_Interface_Transport_Input.md)。生产入口不持有解析参考，analytic_volume/absolute_volume_error 为 NA，accuracy.json 的参考与解析误差为 null，真实体积漂移保留；完成标识为 COMPLETED。其余空间帧、series、来源、失败及接受态序列规则沿用本文。

独立规定速度 App 沿用 [p3b 输入](I12p3b_FirstJet_App_API.md)，本任务接通 `writeFields=true`。`writeFields=false` 现在只关闭空间帧，仍保存 series/meta；旧 p3b 中“无输出目录”的限制属于此前阶段。

```text
Interface
{
    Benchmark
    {
        method hermiteGals;
        problem polynomialVortexXY;
        resolution 8;
        dt .001;
        steps 400;
        initialMesh amr;
        nonpolynomial false;
        Regrid { atSteps (100 200 300); }
        writeFields true;
        writeEverySteps 100;
    }
}
```

`system/run.dict` 仍为 `Run { app Interface; }`。`writeEverySteps` 是 [1,1000000] 整数，启用字段时默认 1，由原 ConfigurationTrace 记录；关闭字段时不允许显式提供无效的 schedule。初始接受态 step=0 和最后接受步始终输出，其余在 step 能被该频率整除时输出。帧编号从 000000 连续递增，不能用帧号推断物理时间。

执行 `bin/QUANTUM.exe --case <case> --check-case` 只预检。正式执行要求 case 内没有 `output`，包括符号链接；已有结果绝不覆盖。为重跑创建新 case 副本，可沿用项目的运行副本脚本。正式输出目录固定：

```text
output/fields/solution.pvd
output/fields/volume/volume_000000.vtu
output/fields/interface/interface_000000.vtp
output/fields/rejected/
output/series/conservation.csv
output/series/solver.csv
output/series/performance.csv
output/series/events.csv
output/meta/frames.csv
output/meta/effective_config.dict
output/meta/configuration_trace.csv
output/meta/manifest.json
output/meta/build_info.txt
output/meta/accuracy.json
output/meta/failure.json                 # 仅失败时
```

ParaView 打开 `solution.pvd`，同一时间有 volume/interface 两组。字段关闭时不生成 PVD/VTU/VTP，frames.csv 只有表头。未提供拒绝候选场快照，因此 rejected 目录为空；数值失败不会生成假接受帧。尚未求解 Hydro，暂不生成 energy.csv。

VTU 用实际 quadtree 叶单元四边形连接关系，包含 `cell_id`、`level`、`phase_volume_fraction`、`geometric_liquid_volume`、`cell_capacity` 和 `phi_center`。体积分数直接由同源几何测量除以 cell capacity，保留原始值，不裁剪、不回写 M。phi_center 从当前 snapshot 的单元中心 chart 查询。XY 的 volume 表示面积，RZ 为完整 2πr 权重体积，具体测度见 manifest。

VTP 每个已有 native graph 用 9 个同源 sample 连成显示 polyline；不另外拟合几何或曲率。点场为 `normal`、`curvature`、`meridional_curvature`、`azimuthal_curvature`、`level_set_residual`、`native_cell_id`，线场为 `component_id`。法向由负相指向正相，组合曲率为 +div(n)，RZ 轴端使用 publication 原有极限。显示折线不是求解器内部的精确曲线，不得用于计算验收面积或替代场的曲率。

两个 VTK 文件都保存 `TimeValue`、`geometry_identity`、`accepted_step`、`mesh_revision`。frames.csv 和 PVD 登记的时间来自实际 accepted state；网格 revision 仅在真实叶网格内容改变时递增，网格内容 SHA-256 另保存在 manifest，不用内存地址冒充持久 identity。

conservation.csv 保存每个接受态（含初始）的同源体积、V0、(V−V0)/V0、解析总体积与绝对误差、cut cells、graphs 和 kappaMax。库存、相通量收支以及尚无 Hydro 的质量/动量/散度量为 **NA**。这里的逐步解析体积误差用 raw geometric total；accuracy.json 的最终 volume error 保留旧基准逐 cell 显示量求和的运算顺序，二者可能有舍入末位差别。最终位置、梯度、曲率及峰值漂移也保存在 accuracy.json。

solver.csv 的 accepted/rejected 是 App 状态，线性/非线性迭代及线性残差为 NA；`physical_residual` 记录 publication 的最大界面残差，按本地 min(cell width) 归一化。events.csv 保存真实重网格 kind、细化/粗化/不变单元数及按**全域容量**归一化的体积变化；这与 conservation.csv 按 **V0** 归一化的时间漂移不同。

performance.csv 每个接受态记录 `interface_initialization`（step=0）或 `interface_step`、`interface_mesh_transfer`、`io`、`step_total` 的实际 steady-clock 秒数，以及 leaf/native-face/cut/graph/primary-vertex/primary-edge/center 数量；未求解的 Hydro/Electricity/nnz 为 NA。未发生网格切换时该阶段耗时为 0。step_total 不含最后的 performance CSV 和 manifest 序列化；当前是 App 阶段统计，核内 trace/closure 细分与规模研究属于 I12p3d，不能据此声称效率 gate 已完成。

effective_config.dict 按 Run 与 Interface 块展示真实配置，原始输入仍在各自 system 文件；configuration_trace.csv 保存原始/生效值、来源位置、默认标记、单位及 consumer。manifest 包含独立 run_id、UTC 开始时间、case 路径、编译器/build type、明确的任务版本、输入文件 SHA-256、几何/网格标识、接受步/时间/帧数与退出状态；Git 信息不可得时为 unavailable。输入哈希已用系统 SHA-256 独立验证。

数值准备失败保留前一 accepted state。若数值已接受而 IO 失败，failure.json 明确报告该实际接受时钟，正常帧仅包含已完成写出的接受数据；程序仍非零退出。场文件完成后才登记索引；PVD、frames 和 JSON 分别以完整临时文件替换，避免写到一半破坏旧索引，它们不构成跨多个文件的原子事务。Windows 对原子替换的短暂 permission-denied 最多重试 6 次，总退避 630 ms，记录 IO retry；持续失败仍报错，绝不删除旧 manifest 或切换数值路线。

本输出不含 restart、Hydro 或电场。只读接口见 [InterfaceFirstJetOutput.h](../../include/Apps/Interface/IO/InterfaceFirstJetOutput.h)，验收见 [报告](../task/I12p3c_FirstJet_Output_Report.md)。
