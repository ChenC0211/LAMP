# I12a11r21 / 几何精度闭合与有界通量的平移驱动接线

目标：只解决半步/全步完整推进的同源接线事实，将 r17/r20 闭合和 r18/r19 共享面通量共同用于原 I12 恒定平移驱动，验证真实 first-jet 输运、守恒更新及发布。

允许修改：`tests/interface_hermite_gals_benchmark.cpp`、`tests/CMakeLists.txt`（0 production + 2 test-related）；本文、Report、README、驱动配置说明及证据。复用现有编排函数，不复制 tracer/mixed/bubble/库存算法。

禁止修改：production 算子、Core、守恒/舍入容限、冻结文件、原 I12 收敛 gate、rotation/vortex 的旧路径及其它物理模块。不更换已完成 r20 的方法；失败不退回几何通量或放宽预算。

数学约束：半步和全步的面量均来自同一旧 accepted state、同一 prepared 库存区域、owned velocity、实际 map 和入流历史；master first jets/center 按原回溯从旧几何输运，M 从唯一共享区间面量更新。原 I10/I11 全部完成后才替换 accepted handle；old state 和失败隔离保持。没有解析内部重置。

实现要求：显式选择 bounded 库存契约，并复用原驱动供旧 A/B。此项为无量纲单位方域 benchmark，预先固定 `maximumHessianComponentChange=min(hx,hy)^2`，系数为 1；几何容限保持默认 1e-10。该 O(h^2) 新增 Hessian 扰动预算随加密缩小，不使用解析误差或失败结果调节；它是验证配置，尚非生产字典默认值。记录实际逐阶段最大请求预算、几何差/pure 差、实际 Hessian 扰动、修正计数及半步/全步库存/面修正诊断。零速度与入流仍明确处理。

验收：T0 原平移/平面跨面/相关 owner 契约输出不变；无 warning。T1 新驱动 XY/RZ×uniform/AMR 的零速 4 步和 Q2 平移 40 步到 t=.4 全部完成；原静止、ledger/reflux 门槛保持，Q2 position/gradient/curvature 均 <=1e-10（原 Q2 控制精度），逐 cell 解析库存相对差 <=5e-12（原平面控制门槛）。48 组任意步内平面跨面用原解析验收。新增 CTest smoke 用同一配置。可记录非多项式探针，但不以局部接线代替原完整 h/dt/曲率矩阵，后者下一项执行。

失败规则：任何数值/身份/门槛失败即 FAIL 并保存阶段与已接受时刻，缺必要 API 才 BLOCKED。不得用未执行矩阵宣告 I12 通过。

交付：驱动接线、实际配置/完整过程诊断、旧/新对照及回归、源码/冻结/hash 审计、报告和 archive。
