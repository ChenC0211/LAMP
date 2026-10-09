# H08a2r13 / benchmark来源观察器退役旧patch依赖

状态：**PASS（观察器与独立参考）**，见[报告](H08a2r13_Benchmark_Source_Observer_Report.md)。原37表/733364行保持，完整双端/五点物理查询及14负控通过，9/9回归；production 0。下一唯一任务恢复[r11](H08a2r11_Common_Regular_Source_Deployment_TODO.md)。来源为[r11构建审计](H08a2r11_Common_Regular_Source_Deployment_Report.md)，r12消费者预算已PASS。

目标：研究驱动观察当前公共物理source，解除对旧tangent patch/rank/slope实现的依赖，使同一公共owner换为正则场后完整研究入口仍可构建和独立验收。

预算：production 0；research≤1（既有cases/benchmarks/cch_convection/driver.cpp），相关参考/短契约≤2。实际依赖超出先再拆任务。无新的driver/target、临时生产兼容壳或重复速度拟合。

先审计spatial/closure观察器及所有CSV读者；旧rank/support/slope观察与实际物理值/完整atom来源严格区别。退役没有数值用途的旧导出，保存历史基线；当前来源的观察应使用实际完整多项式和现有身份/准备报告，不编造正则场的旧patch权重。参考端以独立物理积分验证原208basis、零gross、反号、selected、tiny和AMR，原完整研究输入与PDE门槛保持。

T0标准研究构建/目录/依赖；T1所有改动schema的读回及原行覆盖核验，未改表逐字保持（计时单列）；T2对应现有研究模式和完整次数负控。交付接口/CSV退役清单、独立参考、原范围对照，不把编译通过称为PDE通过。

通过后恢复[r11](H08a2r11_Common_Regular_Source_Deployment_TODO.md)同步实际来源与运动证书，再跑原最终native/affine及H10b-v2 28窗口/84gate。H11未启动。
