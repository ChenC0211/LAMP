# I12p3b / 纯 Hermite–GALS 独立 App

**历史 API，2026-10-03 当前为 [G04 五文件输入](G04_Interface_Dictionary_Schema.md)。** 下文记录当时的命名 benchmark 入口和验证证据，不作为当前可运行 schema。

`bin/QUANTUM.exe` 已支持显式 `Interface.Benchmark.method hermiteGals`。这是六类规定速度的独立 Interface 验证入口，复用 [p3a 时间步](I12p3a_FirstJet_Step_API.md)和 [p2b 动态网格](I12p2b_FirstJet_AMR_Transfer_API.md)，不求解 Hydro/Electricity。它不携带独立 M，不拟合守恒体积；库存及相通量收支为 NA。

在一个新的 case 目录创建 `system/run.dict`：

```text
Run { app Interface; }
```

`system/interface.dict`：

```text
Interface
{
    Benchmark
    {
        method hermiteGals;
        problem translationCircleXY;
        resolution 8;
        dt .01;
        steps 40;
        initialMesh uniform;
        nonpolynomial false;
        Regrid { atSteps (10 20 30); }
        writeFields false;
    }
}
```

从项目根目录执行，替换 `path/to/case`：

```powershell
bin/QUANTUM.exe --case path/to/case --check-case
bin/QUANTUM.exe --case path/to/case
```

预检只解析、验证并报告配置，不建网格、不生成输出。它不能代替实际步的 CFL、几何或边界事件检查。正常退出为 0，输入或运行失败沿用 InterfaceApp 非零退出；错误包含 attempted step、最后 accepted step/time 和失败阶段。

| 字段 | 要求与 consumer |
| --- | --- |
| `method` | 必须显式 `hermiteGals`，选择 FirstJet；旧 `conservativeGals` 保持旧含义 |
| `problem` | 下表六个值之一，选择解析初值、外边界数据、规定速度及误差参考 |
| `resolution` | 整数 [4,256]，单位域初始 Nx=Ny |
| `dt` | 有限正值，`dt*steps` 有限；原 tracer 的半单元 CFL 仍实际检查 |
| `steps` | 整数 [1,1000000]，实际 End 层步数 |
| `initialMesh` | 可省略，默认 `uniform`；另可选 `amr` |
| `nonpolynomial` | 可省略，默认 `false`；`true` 将初值/边界/reference 的 Q2 标量 q 改为 exp(q)−1 |
| `Regrid.atSteps` | 可省略表示无网格切换；整数严格递增、正且小于 steps，在对应物理步完成后切换 |
| `writeFields` | 本阶段必填 `false`；`true` 明确拒绝，正式输出属于 I12p3c |

每个字段和默认值均使用已有 `ConfigurationTrace` 记录来源、单位、consumer。stdout 报告所有生效值和规定物理参数。`system` 只允许上述两个文件；未知字段/块、重复字段、错误布尔/数值、未消费的 legacy 选项均拒绝。原 `cfl/endTime/Adaptation` 不能用于本入口，不能把目标 schema 的平层 `Interface.method` 当作已接通输入。

计算域为 [0,1]²，code units；XY 用面积 dx dy，RZ 用轴对称体积 2πr dr dz。AMR 为原基准的一级中心框 [0.25,0.75]²；每次 regrid 交替 uniform/AMR。uniform→AMR 用 exact refinement，AMR→uniform 显式使用 FirstJetProjection，不能称粗化为精确恒等变换。

| problem | 初始圆/球参数 | 规定速度 |
| --- | --- | --- |
| `staticCircleXY` | 中心 (0.43,0.51)，半径 0.19 | 0 |
| `staticSphereRZ` | 轴上中心 (0,0.51)，半径 0.27 | 0 |
| `translationCircleXY` | 同 staticCircleXY | (0.2,0) |
| `translationSphereRZ` | 同 staticSphereRZ | (0,0.2) |
| `solidRotationXY` | 中心 (0.65,0.53)，半径 0.13 | 绕 (0.5,0.5)，角速度 2 |
| `polynomialVortexXY` | 同 solidRotationXY | 40 cos(πt/0.4) ρ²(−y′,x′)，x′=x−0.5、y′=y−0.5、ρ²=x′²+y′² |

输出 `vortexPeriod=0` 表示其它问题不使用此参数，不是零周期运动。解析场只初始化一次，并为模块允许的真实外边界事件提供 first jet；内部输运从旧 publication 回溯，不能重置为解析解。初始 master/hanging 归属由 Core skeleton 查询，mixed closure、tracer/J、pullback、publication 和 mesh transfer 均调用现有模块。

App 唯一持有 accepted mesh/publication/clock。若该步安排重网格，物理步与零时间 transfer 一并准备成功后才提交；失败保留原接受态。接受时间沿用旧基准的逐步浮点加法。stdout 的 `[Step]` 只在提交后出现，`FIRST_JET_APP_REGRID` 记录实际转移；`FIRST_JET_APP_RESULT` 给出位置、梯度、曲率、解析体积误差和最大相对初值体积漂移。重网格体积变化按**整个计算域容量**归一化，时间漂移按**初始液相几何体积 V0**归一化，两者不能混用。

本任务的 `writeFields=false` 不生成 `output`/`cache`。尚无新路线 VTK/series/meta 文件或 restart。验证结论及范围见 [报告](../task/I12p3b_FirstJet_App_Entry_Report.md)；完整 Gate I 仍需后续输出与效率验证。
