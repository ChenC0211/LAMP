# H10d / 全部native算子一致性与失败候选复核

2026-10-04，**H10d FAIL / H10未关闭**。已独立定位全部native误差，保持原内部affine反例修正；周期闭合/全中心P2候选不相容，硬C/P1加P2惩罚候选仍有加密系数增长并在最细准备失败。候选实现已归档退役，生产恢复H10c精确源码；只有被动研究观察器保留。下一项[H10e](H10e_Bounded_Operator_Consistency_TODO.md)，不进入H11/H08a。

## 现行算子的独立定位

[任务书](H10d_Final_Operator_Consistency_Closure_TODO.md)要求先定位再实现。[rows/rows32](../../cases/benchmarks/cch_convection/driver.cpp)独立组装原K、共享J、实际研究周期端口及uR，与真实apply核对≤1e−12，实际最大1.17e−15；内部、外部、固定轴的平方误差贡献必须加回整体。原四组合n=8/16/32的15,344个native行全部保留，固定轴约束0不删除其行/质量。试验前与最终恢复后的全部逐行CSV字节相同。详见[逐行统计](../../runs/baselines/H10d_20261004/row_summary.json)与对应CSV。

uniform内部L2截断贡献近二阶；外部贡献约1.5阶，normal velocity两个周期重复采样拥有半dual却各自单边导数，峰值误差约O(h)。AMR内部三次余项也随加密放大，不能把全部失败归给边界：

| 原算子内部指标 | n=8 | n=16 | n=32 |
| --- | ---: | ---: | ---: |
| XY 2:1 L2贡献 | 0.0789495 | 0.0257895 | 0.0102056 |
| XY 2:1 无量纲行尺度 | 2.00677 | 2.72131 | 7.33570 |
| XY 2:1 三次导数矩/h² | 6.28235 | 11.43410 | 45.35263 |
| RZ 2:1 L2贡献 | 0.0836741 | 0.0279489 | 0.0116139 |
| RZ 2:1 无量纲行尺度 | 2.07981 | 2.82104 | 8.82787 |
| RZ 2:1 三次导数矩/h² | 6.42860 | 12.27720 | 50.62169 |

上述L2是分类平方误差的平方根贡献，不能相加为完整L2。列中的三次量与行尺度是观察器实际wrapped矩，不是可接受状态修正；精确低阶矩不约束其加密增长。

## 已退役候选及可识别性

[实现前设计](H10d_Final_Operator_Consistency_Design.md)将实际matching exterior segment配对、平移image和native重复采样等值关系交给生产算子；one-ring通过真实周期pressure邻接闭合。两个原half-dual及真实质量不删除，等值速度class用其总质量形成共同导数，同时逐行保持`load=-N+uR`和原material-KE/动量总和。无AMR特例、问题名称或研究载荷修正，未改变H08/H09质量/F/R。

精确全中心P2候选在uniform n=8/16完成组装/实际观察器核验，但首个XY 2:1 n=8联合解失败。独立COD检验该实际组装矩阵：1,440约束/3,864pair、rank=1,422，最小二乘最大残差约7.01e−4，残差平方和约1.8e−5；仅720个硬C/P1约束rank=714、残差为0.002441×256eps（约1.39e−16）。记录在[实际rank日志](../../runs/baselines/H10d_20261004/checks/paired_hard_rank/run.log)。这证明本图和目标的全P2硬约束数值不相容，**不证明所有相容空间算子不可能**；C/P1与skew本身仍可共同满足。

随后实现前补充了硬C/P1、P2余项惩罚与同一one-ring最小改动。统一固定`lambda=(h/L)^2`，原硬矩不放松；没有根据AMR或参考误差选择参数。但实际n=8/16内部余项仍变差：

| 惩罚候选内部指标 | n=8 | n=16 |
| --- | ---: | ---: |
| XY 2:1 L2贡献 | 0.256167 | 0.125026 |
| XY 2:1 行尺度 | 7.99544 | 19.62470 |
| XY 2:1 三次矩/h² | 34.35331 | 87.51564 |
| RZ 2:1 L2贡献 | 0.250721 | 0.127797 |
| RZ 2:1 行尺度 | 9.03034 | 23.31795 |
| RZ 2:1 三次矩/h² | 45.21745 | 107.14934 |

[最细试验](../../runs/baselines/H10d_20261004/checks/penalty_rows32/run.log)在首个uniform XY n=32准备达到49,920迭代，求解器报告残差指标为50,505.84×256eps（仅此混合指标不能单独确定哪类矩或系统残差失败）并exit 2。没有形成四组n=32数据或完整轨迹，不能把部分输出当作完整通过。未继续提高预算、调惩罚参数或增加黏度。

## 恢复后的核验与未运行项

按退役规范，两个生产文件恢复为H10c冻结字节，原短契约未改变；候选header/CPP/driver/设计与实际binary保存在[rejected_paired_candidate](../../runs/baselines/H10d_20261004/rejected_paired_candidate)。该配对API不是现行能力，没有作为fallback留在生产。

最终研究driver仅插入被动rowAudit和两个观察模式，原Study、rhs、RK4、full、donor与84个gate文本精确保持；独立[scope核验](../../runs/baselines/H10d_20261004/research_scope.json)可逆移除观察器后等同H10c。README解释输出含义，退出0仅表示观察器核验通过。

原n=8/16八组合内部affine再次通过，dof911为−1.0000000000000742，误差7.41629e−14，全部内部最大1.58492e−13≤1e−10。最细内部affine先前H10c证据及同一恢复生产源码保持，不冒称本轮另跑audit32。[旧短契约](../../runs/baselines/H10d_20261004/LastTest_final.log)124场景/1,356 trace/6,144 load/24实际跨界/19拒绝/4 regrid及1,516独立矩/3,504pair保持，数值SUMMARY与H10c逐字相同。载荷/动量/KE原1e−12及常速度1e−10门槛不改。架构、自测、目录检查[3/3通过](../../runs/baselines/H10d_20261004/ctest_layout.log)，加短契约本轮受影响检查共4/4。

新闭合候选与有界余项前置失败，**本轮未重跑原完整H10b，也未运行新增配对全场景短契约**。最新完整研究仍是[H10c](H10c_Native_Convection_Moments_Report.md)原28窗口/15,116步、84gate中10空间+4donor失败。前次19/19回归是历史证据，本轮没有重新声称19项运行。没有H11、H08a或一般Hodge实现，也没有完整CCH流体App。

## 冻结与下一入口

开始重核H10c的100个payload。结束审计1,042受保护源码中1,041不变，仅研究driver被动观察器变化；2 production和原test字节不变。346 case/fixture资产中344不变，仅driver/README变化。工具、CMake、输入、旧App及冻结CCH PDF保持，证据和当前交付复制于[runs/baselines/H10d_20261004](../../runs/baselines/H10d_20261004)，独立verification记录工程PASS与数学FAIL分开。

H10e从最终native Taylor展开和显式端口出发，先查精确C/P1、受控二次/三次余项与有界one-ring系数能否共同成立。失败候选的精确P2或惩罚不能当成已验证路线；也不能仅凭此次失败自行改变真实质量定义。空间离散通过之后才按用户顺序启动H08a；Pure FirstJet几何独立、不恢复库存/bubble，原材料R/uR与KEP继续保留。
