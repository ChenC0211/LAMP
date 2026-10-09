# Interface 通用规定输运输入

**历史布局，已由 [G04 五文件 schema](G04_Interface_Dictionary_Schema.md)替代。** G03 的通用表达式和 provider 思想保留，下述单文件示例不再作为当前可运行输入。

当前生产 schema，2026-10-03；替代 I12p3b 的 `Interface.Benchmark`。运行命令仍为 `tools\run_case.cmd <算例名称>`。`system/run.dict` 为 `Run { app Interface; }`，`system/interface.dict` 配置通用数据，不读取问题名。

```text
Interface
{
    method hermiteGals;
    coordinates xy;
    resolution 8;
    dt .01;
    steps 30;
    Domain { min (0 0); max (1 1); }
    Constants { cx .65; cy .53; radius .13; omega 2; }
    Initial { phi "(x-cx)*(x-cx)+(y-cy)*(y-cy)-radius*radius"; }
    Velocity
    {
        type polynomial;
        origin (.5 .5);
        scale (1 1);
        Coefficients
        {
            c01 { x "-omega"; y "0"; }
            c10 { x "0"; y "omega"; }
        }
    }
    Boundary { type prescribed; phi "((x-.5)-(cx-.5)*cos(omega*t)+(cy-.5)*sin(omega*t))^2+((y-.5)-(cx-.5)*sin(omega*t)-(cy-.5)*cos(omega*t))^2-radius*radius"; }
    writeFields true;
    writeEverySteps 7;
}
```

此例表示 `u=-omega*(y-.5), v=omega*(x-.5)`。要增加圆、椭圆或其他支持的隐式初始界面，复制配置并修改表达式；不添加问题枚举。选择 `Boundary none` 时没有规定外边界数据，若实际反向轨迹请求边界数据，程序明确拒绝；它不是自动外推/周期边界。

| 字段 | 当前含义 |
| --- | --- |
| `method` | 仅 `hermiteGals`，纯 FirstJet，无守恒库存与 flux ledger |
| `coordinates` | `xy` 使用 `x,y` 与面积；`rz` 使用 `r,z` 与完整轴对称体积，域的径向下界必须为 0 |
| `Domain.min/max` | 必填，两分量有限坐标，逐分量 max > min，支持非单位矩形域 |
| `resolution` | 每方向基础单元数，4..256 |
| `dt/steps` | 有限正 dt，1..1000000 步，有限总时长；实际半单元 CFL 在推进中检查 |
| `Constants` | 可选有限数字常量；必须被表达式引用；变量、pi 和函数名不能作为常量名 |
| `Initial.phi` | 必填，仅空间表达式；负值为液相，梯度从同一表达式链式求导，不能引用 t |
| `Velocity.type constant` | 必填 `value (u v)`；其余 polynomial 字段不允许 |
| `Velocity.type polynomial` | 必填 origin、正 scale、Coefficients；已有不可变 endpoint provider 提供值和梯度 |
| `Boundary.type prescribed` | 必填 phi 空间/时间表达式；只在算子请求外边界到达点/时间数据时使用 |
| `Boundary.type none` | 不允许 phi；未提供外边界数据，不代表零梯度/自动扩展 |
| `initialMesh` | 可选 uniform（记录默认）/amr |
| `Refinement` | initial amr 或有 Regrid 时必填 `min/max/level`；域内矩形，level 1..2；未使用时禁止提供 |
| `Regrid.atSteps` | 可选非空列表，严格递增且在末步前；按计划在基础网格与指定加密网格间切换，细化 exact、粗化 FirstJetProjection |
| `writeFields` | 必填布尔值；false 仍输出 series/meta |
| `writeEverySteps` | true 时可选 1..1000000，默认 1；false 时禁止提供；始终输出初态/末态 |

Q3 polynomial 定义：`xi=(first-origin.x)/scale.x`、`eta=(second-origin.y)/scale.y`；`cij` 乘 `xi^i eta^j`，i,j 各为 0..3，x/y 分量对应第一/第二坐标速度。系数是仅含 `t` 和常量的表达式，在前一/到达时间层求值；时间 slab 的插值和 RK2 沿用现有数值路线。未列出的系数为 0，并记录来源默认。这不是任意空间函数的隐式多项式拟合。

例如规定时变剪切：

```text
Velocity
{
    type polynomial;
    origin (0 0);
    scale (1 1);
    Coefficients { c01 { x ".1*cos(t)"; y "0"; } }
}
```

表达式复用 Core ScalarExpression：`+ - * / ^`、括号、有限数字、pi、sin/cos/tan/exp/log/sqrt；幂指数须为常数，空间偏导在固定时间求值。公式必须整体加引号，最大 4096 字符、嵌套深度 64；不执行脚本或文件操作。解析、未知变量、类型、范围错误携带文件、行号和字段。预检不建立网格/输出；预检不证明全域/全时段表达式有定义，推进期间非有限值或无定义数据仍会拒绝并保存接受前缀。

生产入口不持有解析解，也不执行 benchmark 误差 gate。`conservation.csv` 的 analytic_volume/absolute_volume_error 为 NA，`accuracy.json` 的参考和解析误差为 null；真实几何体积漂移、曲率与性能保留。`effective_config.dict` 导出全部实际数据和默认系数；取其中 Interface 块可作为新 interface.dict 使用。独立 benchmark 仍由 `tools\build.cmd -Preset local-benchmarks` 构建，再用通用启动器运行 `benchmarks/hermite_gals`。

当前范围是规定速度下的几何输运。Hydro/EHD 原生速度 view 耦合、任意外部场 reader、周期边界及一般动态 AMR 策略要通过后续对应 contract 接入，不能用专用 `problem` 分支代替。

既有 tracer 对同时跨两壁的角点反向轨迹仍明确拒绝；提供表达式不会绕过此接受条件。例如两个方向均非零的常速会在某个外角点触发该限制。这项行为保持原数值路线，未在输入迁移中扩充。
