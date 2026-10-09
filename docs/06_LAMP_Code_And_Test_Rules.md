# LAMP 代码、测试与退役约束

状态：现行工程约束，2026-10-02；补充 [架构规范](02_LAMP_Architecture_Rules.md)。

生产依赖方向固定为 Core → Solvers → Modules → Apps。内核不读字典、不写文件、不拥有 App 主时钟；App 不复制数值算法。每个物理量、几何和接受状态只有一个 owner。

当前 Interface 入口只有纯 Hermite–GALS / FirstJet：输运 `phi + grad(phi)`，几何体积为诊断，不执行守恒库存闭合或失败后切换其他路线。旧独立 levelSet / conservative / conservativeGals / parametricGals 入口已退役。

生产 App 不得用 `problem` 枚举、case 名或目录名选择专用初值、速度、边界、解析解和精度断言。通用数据 provider 与算法选择可以进入生产；具体问题的数据来自配置，解析参考和 benchmark 验收归 `tests` / `cases/benchmarks`。当前独立 Interface 使用 [五文件 schema](api/G04_Interface_Dictionary_Schema.md)，增加支持范围内的新工况只改字典。引入新物理模型/provider 才需要修改源码并定义通用 contract；不能把一个新算例当成新模型。

`Modules/Interface/Compatibility` 和 `lamp_interface_compatibility` 仅维持现有 Hydro/EHD 调用；其入口边逐条记录在架构检查中。新增消费者、GALS 主线依赖兼容层、静默 fallback 均禁止。新 GALS–Hydro 耦合达到对应 gate 后，按消费者逐项移除兼容接口。仍可复用的通用几何测量不是第二条输运路线。

每个长期源码/工具必须说明所属模块、职责及调用入口。按 owner 和数据流拆文件；仅在职责独立或编译边界必要时拆分，禁止每次调试都新增一个变体文件。删除无调用者的失败实现及其专属检查；失败证据留在归档。项目/程序名为 **LAMP — Levelsetbased Adaptive Multiphysics Platform**；现存 C++ `Quantum` 名字空间及 `quantum_*` 库符号暂保留源码兼容，不表示另一程序。

| 检查类型 | 固定归属与运行规则 |
| --- | --- |
| 短契约/单元 | `tests/contracts`，显式 CMake 注册；默认 CTest；有独立断言，通常秒级 |
| 最小 App/CLI/输出 | `tests/contracts` + `tests/fixtures`；复制到 `runs/checks` 后检查 |
| 完整时间窗/空间、时间收敛/论文研究 | `cases/benchmarks`，默认不加入 CTest；按需构建研究驱动 |
| 完整物理算例 | `cases`，通用启动器运行；README 说明目标、参考、验收、实际能力状态 |
| 临时定位 | `runs/checks/<任务>/<时间戳>`；关闭任务时删除或归档，不能变成生产开关 |

研究驱动只调用同一生产算子；参考值不得参与接受状态、误差修正或容差调节。回归检查的“能运行”、精度、收敛与物理验收分别记录。缺失 fixture 标注缺失并禁用，不计作通过；禁止永久堆积无配置的幽灵测试。

提交/交付前：配置和构建、目录与依赖检查、受影响的短测试必须通过；只有算法/误差预算改变时追加对应完整 benchmark。任务报告列明通过项、未运行项、缺失项与退役范围。任何未通过项都不能靠放宽数值断言消除。
