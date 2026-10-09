# H08a2r13 / benchmark公共来源观察器报告

结论：**PASS（观察器及独立参考）**，2026-10-05。关闭 r11 的研究驱动旧 patch 构建依赖。production 0、短契约 0；仅修改既有研究 driver，新增同目录独立参考并更新说明。实际生产来源仍为旧 trace/lift；本轮不把隔离候选或观察器通过记作部署/PDE通过。

依据：[TODO](H08a2r13_Benchmark_Source_Observer_TODO.md)、[实施前设计及五点补充](H08a2r13_Benchmark_Source_Observer_Design.md)。材料四硬约束保持；未增加 R=0、quarter/pressure→dual交换要求，Pure FirstJet独立，无库存/bubble。

## 改动与退役

| 观察入口 | 本轮处理 |
| --- | --- |
| advector_faces | 删除 rank/deficit/support；物理字段、来源和full clock保持 |
| closure_faces | 删除 mean_x/mean_y/scale/slope；原208非轴单位面通量基、坐标、native邻接和面积保持 |
| closure_atoms | 公共实际atom完整c0…c4/e0…e4、物理长度/面积、双端积分 |
| closure_views | 实际primal/dual atom索引和原邻接，不生成假patch权重 |
| closure_sources | 实际geometry/advector identity、full clock、规模和view数量 |
| closure_samples | 每atom五个公共cellVelocity法向值，GL3节点加两个独立点 |

当前活动 cases/tests/tools 无其他旧CSV读者。旧H10j0/H10l1参考的三系数、旧rank库存及trace/lift slope重算保留在冻结历史中。新 [verify_source.py](../../cases/benchmarks/cch_convection/verify_source.py) 从CSV按完整物理次数独立积分/读回，不读取production weights，不复制拟合。H10i独立几何诊断 `spatial_source_patches.csv` 不依赖production patches，保持；momentum契约的patch指针断言仍由r11部署迁移。

GL3对四次速度乘RZ线性测度精确，但三个值不能唯一检查四次物理迹。因此最终导出使用五个不同值节点；三点过渡运行仅留作中间证据。负控明确构造在GL3三个节点及物理积分均为零的三/四次扰动，最终五点查询会拒绝。未引入新的生产分支、owner、驱动或CMake target。

## 验收

T0：标准 `build.cmd -Preset local-benchmarks -Target LAMP_cch_convection_benchmark` 构建通过。架构/自检/目录及phase、material、momentum、characteristic、dual、KEP共 **9/9** 通过，199.98秒。最终五点仅改变研究导出；生产/测试源码和六个契约可执行文件哈希保持。六个契约原 **10,644** 行逐字保持（原计时字段单列）。包括原4个selective-law反号见证、tiny/multiroot/轴/selected/AMR/F/R/uR/KEP与来源拒绝；闭流函数研究本身不冒称覆盖反号见证。

T1：同一当前生产库的 before/final 对照运行 `advector/closure/boundary/trace-contract`。原 **37表/733,364行** 全保留覆盖：34表逐字相同，2表仅退役上述字段，1表仅准备/更新时间变化。新增4表283,728行；最终共41表1,017,092行。C++差异严格限于两处观察函数，其外全部字节保持，原full输入、RK4、native对流、donor及PDE gates未改。

T2：当前来源与冻结r11候选各运行三模式并用同一个最终参考验收。隔离候选四份来源/运动源码与r11冻结字节一致，当前driver可直接编译；未加入patch兼容壳。原advector 64布局/128更新、56适用解析布局，boundary 32布局/64更新，closure 208基全部通过。

| 物理覆盖/误差 | 当前生产来源 | 隔离正则候选 |
| --- | ---: | ---: |
| closure实际atoms / 五点查询 | 35,968 / 179,840 | 43,680 / 218,400 |
| closure实际view entries | 67,712 | 75,424 |
| 原native faces / cells | 12,016 / 4,912 | 12,016 / 4,912 |
| 严格零gross的原D_I cells | 4,552 | 4,552 |
| basis面均值预算误差 | 1.108e−16 | 7.065e−17 |
| basis cell边界积分预算误差 | 4.423e−17 | 3.511e−17 |
| basis物理速度查询预算误差 | 2.298e−16 | 2.640e−16 |
| advector原解析point最大绝对误差 | 2.221e−16 | 2.332e−15 |
| advector平均散度最大绝对误差 | 1.188e−14 | 1.188e−14 |

预算为完整物理chart与面积/输入工作量，阈值仍1e−12且无常数地板。零gross原D_I仍严格为零，浮点积分余项不被强行置零。平均散度由native incidence独立重算，点散度与旧体积坐标割线不作为新场常值约束。真实view端点、平面、无缝子面、面积、邻接朝向及覆盖逐项检查；完整材料内部反号/唯一F/KEP仍由原生产契约验收。

14项负控全部拒绝：双端/单端截断高项、缺c4列、积分号、宽度、full clock、缺采样、错误物理查询、basis输入号、非法共享atom、邻接号以及上述两个GL3不可见扰动。完整次数有理多项式参考的截断误差为0.05424；两项GL3不可见扰动在新增点产生明显非零误差，见原始JSON。

## 证据与下一步

检查目录：[H08a2r13/20261005_184026_443](../../runs/checks/H08a2r13/20261005_184026_443)。关键文件：`verification.json`、`tables_verification.json`、`reference_binding.json`、`final*_reference.json`、`negative_controls_final.json`、`LastTest_final.log`。最终运行使用 `final_*` / `final_candidate_*`；`after_*` / `candidate_*` 是三点中间快照。所有参考结果绑定最终脚本、CSV及实际runtime哈希。

冻结目录：[H08a2r13_20261005](../../runs/baselines/H08a2r13_20261005)。保存before/final源码和程序、隔离候选、编译选项、原CSV、脚本/日志及manifest；无可执行临时文件进入生产树。

本轮未运行完整full/28窗口84gate，也未实现同源时空积分、FirstJet实际运动桥接或完整流体App。下一唯一任务恢复 [r11同步实际来源及运动证书部署](H08a2r11_Common_Regular_Source_Deployment_TODO.md)，先真实production全部契约，再重跑原内部affine/native及H10b-dissipation-v2 28窗口/84gate，随后返回H08a3。旧4项donor失败不改判；H11未启动。
