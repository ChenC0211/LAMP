# I12a11r11 / 弱灵敏度 cut-cell 闭合残差来源审计

## 目标与范围

只确定 r10 已定位的三个实际 end 闭合（step 1 / cell 648、step 2 / cell 650、step 3 / cell 363）为何产生大系数修正：分解其闭合前 M − predictor volume。0 production + 1 test：`tests/interface_hermite_gals_benchmark.cpp`，另含本 TODO/报告、README 与 `runs/baselines/I12a11r11_20261001/`。禁止改变冻结模块、primary、M、实际 B/J、flux、I06/I10/I11、全部容差或 benchmark gate。

## 数学约束及实现

复用原 XY uniform n=32、dt=.01、四步 driver，显式只读入口。用实际 old PhaseRegion / raw-final integrated faces / assembled M，将差额分为：old M − old geometry、old geometry + raw exchange − shifted geometry、flux repair、assembly rounding、shifted geometry − predicted geometry。逐项有符号相加应还原实际残差，原 64 eps capacity 只作算术分解缺陷检查，不作为放宽闭合条件。

记录旧/平移/预测的两相体积、双 Green 分量、moment estimate、初始 dV/dc，以及随后真实 root 的初始/最终体积、系数和灵敏度，核验诊断预测体积等于实际 root 初始体积。private predictor fixture 仅作为原区域 API 的只读测量适配，不发布、不反馈实际推进。记录 predictor 的 native Q4 DOFs，便于复算。candidate/old identity 与 M 不变，异常单独计数后继续原推进。

## 验收与失败规则

T0：Release 无 warning；11 项直接回归完整输出与 r10 一致；其它 27 可执行、全部 production、81 冻结/参考文件、3 PDF 和 r10 archive 不变。

T1：3 个原指定 cut-cell 残差完整分解、与 3 次实际 root 对齐，所有诊断异常为零；四步最终 BENCHMARK 行逐字等同 r10。报告每个差额的量级及测量误差估计，不把 error estimate 当严格上界，不把相关单元身份直接当因果证明。

证据不完整或诊断扰动 actual solution 则 FAIL；诊断 PASS 不改变 I12 数值 FAIL。只在具体外部依赖不可用时记录 BLOCKED，不因数值问题扩大容差。交付源码、日志/耗时、原始分解、哈希、报告和 archive；再据证据确定最小修复或必要的具体方法修订。
