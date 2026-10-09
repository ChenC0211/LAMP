# Interface 五文件字典 schema

状态：当前冻结输入，2026-10-03。保留 G03 通用表达式与纯 Hermite–GALS/FirstJet，替代其单文件布局。完整模板在 [旋转算例](../../cases/hermite_gals_rotation_xy/system)。以下五份文件均放在 system，运行命令保持 `.\tools\run_case.cmd hermite_gals_rotation_xy`。

## 文件与可运行示例

run.dict 管唯一主时钟和输出，Run 块只选 App：

```text
Run { app Interface; }
Time { units code; dt .01; endTime .3; }
Output { Fields { enabled true; writeEverySteps 7; } }
```

geometry.dict 管计算域、坐标与初始几何：

```text
Geometry
{
    coordinates xy;
    Domain { min (0 0); max (1 1); }
    Initial { phi "(x-cx)*(x-cx)+(y-cy)*(y-cy)-radius*radius"; }
}
```

mesh.dict 管网格，与界面方法无关：

```text
Mesh { resolution 8; initialMesh uniform; }
```

fields.dict 管独立给定速度及输入表达式共享参数：

```text
Fields
{
    Constants { cx .65; cy .53; radius .13; omega 2; }
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
}
```

interface.dict 只管界面方法与边界数据：

```text
Interface
{
    method hermiteGals;
    Boundary
    {
        type prescribed;
        phi "((x-.5)-(cx-.5)*cos(omega*t)+(cy-.5)*sin(omega*t))^2+((y-.5)-(cx-.5)*sin(omega*t)-(cy-.5)*cos(omega*t))^2-radius*radius";
    }
}
```

此例为偏心圆在 `u=-omega*(y-.5), v=omega*(x-.5)` 中旋转。Boundary 表达式只提供实际请求的外边界到达点/时间数据，不覆盖内部推进。

## 字段、范围与默认

| 文件/字段 | 当前契约 |
| --- | --- |
| run: Run.app | 必填 Interface；不在 Run 块内放时间或输出 |
| run: Time.units/dt/endTime | 必填 units code；有限正时间；endTime/dt 须对应整数步数 1..1000000，使用公共 fixedStepCount，steps 为派生量 |
| run: Output.Fields.enabled | 必填 true/false；false 仍保留 series/meta |
| run: Output.Fields.writeEverySteps | enabled=true 时可选 1..1000000，默认 1；false 时禁止提供；初态/末态始终输出，格式固定 ASCII VTK XML |
| geometry: coordinates | 必填 xy/rz；变量 x,y 或 r,z，测度为面积或完整轴对称体积；rz 的径向下界须为 0 |
| geometry: Domain.min/max | 必填有限两分量向量，逐分量有序；支持非单位矩形域 |
| geometry: Initial.phi | 必填空间表达式，不含 t；负值为液相，梯度从同一表达式链式求导 |
| mesh: resolution | 必填 4..256，每方向基础单元数；非正方形域的单元物理尺寸可以不同 |
| mesh: initialMesh | 可选 uniform（默认）/amr |
| mesh: Refinement | initialMesh=amr 或有 Regrid 时必填 min/max/level；域内矩形，level 1..2；否则禁止提供 |
| mesh: Regrid.atSteps | 可选非空正整数列表，严格递增且小于派生总步数；在基础/指定加密网格之间切换，细化 exact、粗化 FirstJetProjection |
| fields: Constants | 可选有限表达式参数，整个工作流仅此一处声明；必须被 Initial、Boundary 或 Velocity 使用，变量/pi/函数名为保留名 |
| fields: Velocity.type constant | 必填 value (u v)，不允许其他 polynomial 字段 |
| fields: Velocity.type polynomial | 必填 origin、正 scale、Coefficients；未写的系数为 0 并记录默认；不可用任意空间表达式隐式拟合 |
| interface: method | 必填 hermiteGals；无 problem、debug、守恒库存、失败 fallback 或占位 Numerics |
| interface: Boundary.type prescribed | 必填 phi 空间/时间表达式；只提供外边界数据 |
| interface: Boundary.type none | 不允许 phi；实际需要外边界数据时拒绝，不表示零梯度、周期或自动扩展 |

Q3 速度使用 `xi=(first-origin.x)/scale.x`、`eta=(second-origin.y)/scale.y`，cij 乘 xi^i eta^j，i,j 各为 0..3。x/y 分量是第一/第二坐标速度，包括 rz；系数表达式只含 t、pi 和 Constants，在前一/到达时间层求值，时间插值及 RK2 沿用原路线。

表达式为 Core ScalarExpression 的有限算术：+ - * / ^、括号、pi、sin/cos/tan/exp/log/sqrt；幂指数为常量，空间偏导固定时间求值。整个公式加引号，最多 4096 字符、嵌套深度 64；不执行脚本、include 或文件操作。Constants 是表达式参数 owner，不保存材料数据、solver 容限、网格或时钟控制。

## 验证与输出

五文件分别校验根块和字段，未知/重复/错放/未使用输入报错；不需要的 materials、hydro、solvers、units 等不能作为占位文件附加到此工作流。预检不生成网格或 output；实际 CFL、表达式在全时段的定义域与数值接受条件在推进时检查。Regrid 与 endTime、Refinement 与 Domain 的跨文件错误定位到 mesh/geometry/run 的真实源字段。

输出记录五份输入的 SHA-256。effective_config.dict 按所属文件汇总真实生效值；output/meta/effective 下的五份同名字典可复制到新运行的 system 后回读。configuration_trace.csv 保留真实源文件/行号、默认值和 consumer，不把所有字段记在 Interface 名下。

生产解析误差为 NA/null，真实体积漂移、曲率和性能保留；COMPLETED 表示配置的推进窗口完成。解析精度 gate 仍在独立 cases/benchmarks。算法、接受条件和误差预算不因文件拆分改变。

## 兼容性与后续扩展

G03 单文件输入及早期 Benchmark/problem 路线均拒绝。24 个模板及最小 fixture 同步迁移；没有重复 key 兼容别名，不按目录或 case 名推断参数。新增配置遵守 [字典约束](../03_LAMP_Dictionary_Rules.md)，同次更新实现、模板、测试、输出和 API。

本入口仍是独立规定速度输运；Hydro/EHD 耦合、任意外部速度 reader、周期边界和一般动态 AMR 策略尚未接入。既有 tracer 仍拒绝角点同时跨两壁的反向轨迹。未实现功能不能通过把字段放到 interface.dict 或提供空字典启用。
