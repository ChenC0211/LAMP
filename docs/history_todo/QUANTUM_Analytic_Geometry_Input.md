> **状态：historical/reference（历史资料，仅供追溯）。**
> 本文原始正文保留为历史记录，其中“当前”“下一步”、TODO、旧规则与通过记录均指写作当时；不再自动产生新任务，也不能证明 V2 gate 通过。
> 现行任务唯一入口：[分模块任务清单](../roadmap/05_LAMP_Module_Task_Roadmap.md)。文档权威见[开发总纲](../00_LAMP_Development_Charter.md)，完整入口见[文档索引](../README.md)。
> 旧字典/职责说明可供查询既有实现；V2 新增或修改仍遵守现行六份规范。历史正文中的链接保留原样，目标缺失时不能作为本机已验证证据。
> 本目录中的旧设计、TeX/PDF 和索引同属历史参考；总纲指定的三份冻结方法 PDF 仍以 `docs/` 下的权威文件为准，不因本目录的历史标记而降级。

# 解析几何字典输入

2026-09-28：参数方程与隐式方程均已接入普通两相入口，共用表达式求值和最终 Hermite 构造。Release 编译及隐式 XY/RZ、参数式振荡三个字典预检通过；尚未进行本版物理运行，不能据此判断曲率收敛或长期输运精度。

## 用法与适用范围

在 geometry.dict 中使用 `ParametricInterface.type analytic`（参数方程）或 `type implicit`（隐式方程）。当前正式消费者为 `cchMomentMidpoint` 的两相 Hydro 入口，配合 `Interface.transport parametricGals`；已有形状输入保留。它定义初始自由界面，不定义实体壁面、运行期强制形状或解析流场。

二维平面采用参数 s∈[0,1]，分别填写 x、y。例如二阶形状扰动液滴：

```text
Domain { min (-0.5, -0.5); max (0.5, 0.5); dimension 2; }
ParametricInterface
{
    name liquid;
    type analytic;
    constants { R 0.1; epsilon 0.05; }
    x "R*(1+epsilon*cos(4*pi*s))*cos(2*pi*s)";
    y "R*(1+epsilon*cos(4*pi*s))*sin(2*pi*s)";
    segments 64;
}
```

name 与 hydro.dict 的 interfaceLevelSet 相同。曲线须光滑、无自交、逆时针闭合；s=0 与 s=1 的位置及一阶导数一致。内部为 negative 相，外部为 positive 相；长度使用代码单位。

轴对称采用 r、z。下例为从下轴极点走向上轴极点的球滴子午线：

```text
Domain { min (0, -1); max (1, 1); dimension 2; }
ParametricInterface
{
    name liquid;
    type analytic;
    constants { R 0.4; }
    r "R*sin(pi*s)";
    z "-R*cos(pi*s)";
    segments 32;
}
```

Hydro 使用 axisymmetricRZ。两端须满足 r=0、dz/ds=0，起点 dr/ds>0、终点 dr/ds<0，其余自由面位于域内；RZ 使用真实旋转测度。轴线闭合段只用于相区域，不成为有表面张力的自由面。附壁及固定供液口初态继续使用已有专用输入。

多分量沿用 type components：每个子块可用 type analytic 或 type implicit 和自身的公式、constants、segments，也可与已有形状并列。常数只在所属分量内生效，不从外层隐式继承。

## 隐式方程与分量选择

二维圆滴可写为：

```text
ParametricInterface
{
    name liquid;
    type implicit;
    constants { R 0.4; }
    F "x^2+y^2-R^2";
    seed (0.4, 0);
    traceTolerance 1e-10;
    segments 64;
}
```

这里的计算域应包含整个圆，例如前述静滴的 [−1,1]²。`F` 填写等式左侧，不写 `=0`；F 在所选边界内侧须为负、外侧为正。`seed` 是目标边界上或非常接近边界的一点，用来选择这一条连通边界。程序只追踪这个分支，不自动发现同一公式中的其他零集；多个分量分别放入 components 子块并各给一个种子点。

轴对称球滴使用以下形状块，计算域 r∈[0,1]、z∈[−1,1]：

```text
ParametricInterface
{
    name liquid;
    type implicit;
    constants { R 0.4; }
    F "r^2+z^2-R^2";
    seed (0, -0.4);
    traceTolerance 1e-10;
    segments 32;
}
```

RZ 的 seed 必须位于下轴极点，r=0；曲线经过 r>0 区域，到上轴极点结束。两极要求 F_r=0、F_z≠0 和正确方向。公式须在零集附近光滑，包括轴附近用于求值的延拓区域。当前支持规则闭合 XY 边界或带两个轴极点的 RZ 子午线；不支持带孔液域、相交零集、梯度消失的拓扑临界初态、附壁隐式边界或自动枚举所有分量。

`traceTolerance` 是初始零集追踪的长度容差，省略时为较短域边长的 1e−9。它控制步进估计和 |F|/|∇F| 的法向修正残差，不是 F 的绝对数值容差；单纯缩放 F 不应改变几何定义。输入过小而低于坐标舍入尺度、过大而不能分辨局部曲线，均会拒绝。

种子点的初始法向修正限制在较短域边长的 1/32 内。追踪步长自适应缩小，最多 20000 个接受步；局部梯度或分支不能分辨时明确报错，不改种子或任选另一条曲线。适用性仍以用户声明的规则、充分分辨的分量为前提，这不是任意解析函数全部零集的拓扑认证。

## 公式规则

| 项目 | 规则 |
| --- | --- |
| 变量及常数 | 参数式用 s；隐式 XY 用 x、y，隐式 RZ 用 r、z；均有 pi，区分大小写 |
| 用户常数 | 可选 constants 块；值为有限数字，名称以字母开头，随后可含字母、数字和下划线，不能覆盖内置名称 |
| 运算 | +、−、*、/、^、括号及正负号；显式写乘号 |
| 幂 | 指数必须为常量表达式；-s^2 表示 -(s^2)，连续幂右结合 |
| 函数 | sin、cos、tan、exp、log、sqrt；三角函数使用弧度 |
| 定义域 | 不允许除零、非有限值或导数；log、sqrt 参数须严格为正；负底数非整数幂不合法 |
| 字符串 | 坐标公式用引号包住；无赋值、条件分支、文件访问或脚本执行 |
| 初始分辨率 | segments 为 8–512 个参数弧段，与流体网格分别设置 |

未知名字、重复常数、缺少坐标、坐标名与维度不符以及混用其他形状参数会报错。闭合点或轴极点只容许浮点舍入量级的误差，不把一般未闭合曲线强行封口。

## 与数值方法的关系

公式的值和一阶导数由同一个表达式按链式法则计算，无需用户填写导数，也不用差分估计。参数式直接给出位置及切向；隐式式沿 T=(−F_y,F_x)/|∇F| 追踪，使用四阶 Runge–Kutta 的步长加倍估计和法向 Newton 修正，并检测回到起始截面或抵达上轴极点的事件。这个几何初值计算没有物理时间，不是另一套界面输运方法。

两种入口最终共用端点位置、端点切向和中点位置的四次 Hermite 构造。隐式追踪的暂存点随后释放；后续输运不再调用初始公式。

有限段曲线是公式的近似。初始 FV 库存由这条实际曲线建立，不用解析体积覆盖；它与解析体积的偏差属于初始化误差，之后的库存漂移属于输运误差。segments 加密是几何收敛的一部分；traceTolerance 不控制最后有限段曲线的全部表示误差，二者应分别判断。预检检查表达式求值、端点条件和生成多项式的域内范围，不证明任意公式在整个参数区间光滑、无自交或已充分分辨。

## 完整示例与运行

示例见[平面静滴](../cases/hydro_static_drop/xy/README.md)、[轴对称静滴](../cases/hydro_static_drop/rz/README.md)和[无黏振荡](../cases/hydro_drop_oscillation/README.md)。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/build.ps1 -Target QUANTUM
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/run_case.ps1 -Case hydro_drop_oscillation -CheckOnly
```

去掉 -CheckOnly 才会推进物理时间；本轮仍按主计划集中开跑。结果保存在 cases/<case>/output/<时间戳>/，包括原始公式字典和程序身份。初始速度为零，压力由首步联立方程建立；初始压力占位列不能作为已有物理解。

当前版本编译记录见[implicit_geometry_compile.log](../build/implicit_geometry_compile.log)，退出 0；程序 SHA256 为 F34C6C5EB46F26B60275B0E95E98C69A2C90F5379B68694616C2FAA6AD634733。输入证据：[隐式 XY](../cases/hydro_static_drop/xy/output/20260928_233051_762_implicit_input_5d7b3c/run.log)、[隐式 RZ](../cases/hydro_static_drop/rz/output/20260928_233052_458_implicit_input_fe6204/run.log)、[参数式振荡](../cases/hydro_drop_oscillation/output/20260928_233053_091_implicit_input_88608e/run.log)。这些仅是初始化入口的预检，未运行 CTest 或物理时窗。
