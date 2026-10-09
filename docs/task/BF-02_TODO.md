# BF-02 / C 共用的恒曲率 CSF 载荷

状态：**PASS（仅BF-02）**。详见 [报告](BF-02_Report.md)：新两组契约86.6855/9.58667秒，200次实际saddle及10/10回归通过；9,528条cell/22,002条face数据与3/2范围hash归档。本轮只执行BF-02；下一次允许仅执行BF-03，完整BF-CSF及旧H13后续gate仍未通过，BF-03/H14/H15本轮未启动。

任务 ID / 标题：BF-02 / 同一压力C上的恒曲率实际CSF动量RHS。

目标：从已验收的同源color生成 `g=C*c` 与 `S_sigma=+sigma*kappa0*g`，以实际闭域鞍点验证符号、静态自由行平衡、真体积gauge和壁反力。

允许修改：新 `include/Modules/Hydrodynamics/CchBalancedCsfLoad.h`、新 `src/Modules/Hydrodynamics/CchBalancedCsfLoad.cpp`、根 `CMakeLists.txt` 单一CCH注册，共3 production。新 `tests/contracts/hydrodynamics/cch_balanced_csf_load_contract.cpp` 与 `tests/CMakeLists.txt` 注册，共2测试。API/TODO/报告和现行状态索引；输出位于 `runs/checks/BF-02/20261009_160120`，全在D:。

禁止修改：BF-01 color、Core/Interface、旧H12/H13、C/D、H_rho/质量、F_rho/R/uR、KEP、黏性、压力/鞍点及原物理/数值/180秒门槛。不得启动真实变量曲率映射或时间推进。

数学约束：唯一C.apply(color)，实际RHS为正号；sigma有限非负、kappa0有限，sigma=0或kappa0=0精确零。自由面 `C(sigma*kappa0*c)-S=0` 到原1e-12尺度；真实AMR子面/RZ测度只由C承担。闭域每压力分量按真体积去均值，fixed rows保留reaction，不能假定壁C1=0。所有mesh/layout/pressure owner/color owner/geometry/full clock同源，错配拒绝。

实现要求：新不可变prepared owner保留原CchPressureLoad和CchRegularizedColor，只调用已有typed CchNativeValues/apply。暴露原color梯度与实际momentumRhsLoad、sigma/kappa及来源检查；不增梯度stencil、不复制几何、不输出新状态、不读写文件。

验收：T0 Release无新warning、architecture/self-test/repository_layout及精确文件预算；T1独立真实面incidence核验全部rows和虚功，随机平面offset/正标量输入标度、圆/球、XY/RZ、uniform/2:1 AMR、selected/annular/disconnected组件，sigma/kappa零、正/负kappa及sign负控；T2实际旧鞍点搭配原sharp两相mass/黏性，继承H12静态速度<=1e-11、压力<=1e-10、自由行balance<=1e-11/代数load<=1e-12，rho比1/1e3/1e6含mu=0，无floor；错误owner/clock/regrid拒绝。保存逐cell压力/color与完整逐面gradient/RHS/free residual/wall reaction。相关旧H12和BF-01回归不变。

失败规则：错误符号、重复面积/密度因子、错误边界行、任意新stencil、来源错配放行或原门槛未过均FAIL；不可通过解析压力写入solver、加人工黏性、改color/旧消费者、放宽门槛或隐藏例外换PASS。真正超本任务范围的依赖记录BLOCKED，正常修正仍在同BF-ID进行。

交付：3/2精确文件和基线diff、唯一owner/API、实际构建/原始测试、真实全row/pressure数据及hash归档、未解决项、PASS/FAIL/BLOCKED及下一次是否允许仅执行BF-03的明确判断。
