# I12a12 / 同源多项式速度区间源

目标：只解决一般运动的速度数据来源及连续区间包围事实，使实际 RK2 endpoint/substage views 与后续事件覆盖使用同一份不可变数据；不以黑盒采样一致性冒充事件认证。

允许修改：新增 `include/Core/Geometry/PolynomialVelocitySlab.h`、`src/Core/Geometry/PolynomialVelocitySlab.cpp`；现有 `HermiteGalsTemporalPhaseFlux.h/.cpp`；根 `CMakeLists.txt`（严格浮点 TU 接入）。共 5 production/configuration 文件。新增 `tests/interface_polynomial_velocity_contract.cpp` 与 `tests/CMakeLists.txt`，共 2 test-related。本 TODO、设计、API/报告、README 及本项证据。

禁止修改：现有 Core tracer 的 RK2、B/J、边界/CFL 规则；geometry、M、bubble、inventory/reflux；旧 constant/opaque provider 的语义；冻结文件、原 gate、既有 benchmark 参数。无新全局速度未知量、无 Hydro 接线、无失败 fallback。

数学约束：速度为固定坐标基底中的 tensor cubic 空间多项式，时间上在 owned old/end 两份系数之间线性插值；值及空间一阶导数来自同一 Horner 运算。使用既有 outward interval arithmetic 对闭空间矩形×闭时间段包围同一速度/导数。结果只包围速度，不提前认证 characteristic、相事件或守恒通量。所有数据深拷贝，实际 substage 保持 old/source owner 和真实 candidate/time identity。

实现要求：先写设计。Core 提供中立的 owned polynomial slab 求值/包围；加入既有 strict/no-LTO geometry 编译对象，复用现有区间算术，不在 Interface/App 复制舍入内核。Interface 显式工厂构造两个实际 endpoint views，公开只读 polynomial capability 和同源实际 substage view。opaque callback 即使数值恰好是多项式，也没有此 capability。特殊 t=old/end/midpoint 分支与范围包围都必须处理；不引入经验容限掩盖包围失败。

验收：T0 Release 无 warning；architecture 与已有 interval、trace、vertex、edge、平移 smoke、跨面、有界库存等相关回归保持旧数值输出。T1 覆盖 affine rotation、空间 cubic differential rotation、时间变号、RZ 轴向/线性应变、非单位坐标基底；独立解析值/导数、实际 RK2 B/J 与原 callback A/B、三档 dt 的旋转特征二阶收敛、空间有限差分 J 检验。包围覆盖解析内部极值及密集独立参考样本（采样仅验证实现，证明来自区间运算）；包含中点/端点、缩小 box、overflow/非有限/非法 time/scale/rounding-mode、moved-from、别的 source/candidate、input mutation/销毁后的 lifetime 与 accepted M/identity 不变。

失败规则：任何原回归、包围或 coherence 检验失败记 FAIL 并保留数据；不得扩大原容限/改 tracer。确实缺少接口时记录具体依赖，但不得把本算子通过当作 general phase flux 或 I12 通过。

交付：源实现、可执行测试、逐项数据、源码/冻结/runtime/旧输出审计、API 与报告及 archive；下一项才实现同源 RK2 空间—时间包围/事件覆盖，不并行扩张本项。

完成记录（2026-10-02）：**PASS（本项）**。8 类来源、5832 个解析样本、32 次实际 B/J A/B 与 45 项拒绝/契约检查通过；20 项回归通过，18 项旧数值输出不变。见[报告](I12a12_Polynomial_Velocity_Interval_Source_Report.md)。一般相通量与 I12 / Gate I 仍未通过。
