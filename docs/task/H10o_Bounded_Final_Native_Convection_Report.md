# H10o / 有界最终 native 对流算子报告

2026-10-05，**H10o PASS：显式闭合端口的实际算子与 H10b-dissipation-v2 完整研究通过**。原不透明 trace 入口及其内部 affine 回归保持。H10 的本次空间/KE gate 已收口；下一单任务为用户指定的 [H08a 同源时空材料通量](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)。尚无完整 CCH Euler/NS App，不把本次规定速度研究当作任意切界面耦合或一般生产边界的 PDE 验收；H11 未启动。

## 实际改动与四材料硬约束

按[任务书](H10o_Bounded_Final_Native_Convection_TODO.md)和[实施前设计](H10o_Bounded_Final_Native_Convection_Design_Change.md)，仅修改 CchKepConvection.h/.cpp、一个既有 KEP 短契约及原研究 driver。原 H08/H09 的 F、真实三层质量、R、E_sp/G、uR 公式及 FirstJet 全部保持。材料仍只有唯一共享 F、内部反号/AMR 子面和、来源一致、KEP/skew 四硬约束，不恢复 quarter、交换恒等或 R=0，更没有库存、bubble 或几何体积拟合。

闭合入口接收显式段配对和平移关系，由实际 pressure one-ring 生成最终共享 skew 对。一次矩使用真实 M 与同一 actual B；不再把 flux-centroid 矩当点 PDE 目标。常数/可识别一次矩、二次余项界和行 l1 界作用于同一个最终 K。以中心共享系数为最小改动起点，已有 Eigen 稀疏 LDLT/ADMM 与活动约束精化求解；全局共享保证反号，每对支撑仍局限于事前 geometry one-ring。没有 AMR/问题名称特例、耗散或失败 fallback。

二次移位矩界为 s*h/L，行 sum|K| 界为 2s，s 为原实际面流量绝对值之和。独立 OSQP 只在 runs/checks 中验证可行性，不成为项目依赖，也不把求解系数喂给生产。L1 约束通过支撑半空间验证/加入同一凸问题；实际四网格族没有触及 l1 界。局部 SVD 检查可识别性，最终复算原尺度残差，不能以迭代停止条件代替验收。

周期 normal 采样按声明关系等值；两份原 half-dual 质量保留。按其总质量形成同一加速度后逐原 native 分配载荷，代数上仍为 -N+uR。所有边界/轴行继续计入误差；不猜测周期性，不丢弃原生自由度。错误配对、非零未配对流量、不同 alias 值与原来源/clock/mesh 错误明确拒绝。

## 原反例、独立矩阵和完整研究

原 n8/16/32 内部 affine 先于完整研究完成，15,344 原行全部保存，内部最大误差 **9.335e-13 <= 1e-10**。非周期 affine 用原显式解析 trace 入口，未错误套用周期关系。新闭合入口另检查可允许的 transverse affine，以及独立最终一次/二次矩。

新算子的独立矩阵审计保留 **15,344 native 行**，包括重复 half-dual、外边界和固定轴；矩阵与真实 apply 最大尺度误差 **1.341e-15**。分类误差平方贡献加回整体。n32 AMR 原尺度常数/一次残差均小于 5e-16；相对于 2s 的最大 l1 比为 XY 0.606861、RZ 0.635787。最大支撑距离/h 为 7.017835，来自原 one-ring 允许对，未扩大图。

原 **28 窗口 / 15,116 RK4 步 / 928 保存时序行**完成，进程 exit0。`H10b-dissipation-v2` **84/84 gate PASS**；原四项累计 donor 损失阶仍在 legacy_gates.csv 中 FAIL。其余80项、网格 n8/16/32、T=1、32n步及独立20/40/80对640步时间研究的源码定义保持，未调整精度门槛。

| 网格族 | 截断误差阶 n16→32 | 终值阶 n16→32 | n32 终值 L2 | n32 终值 Linf |
| --- | ---: | ---: | ---: | ---: |
| XY uniform | 1.991655 | 1.990209 | 0.04029228 | 0.04028410 |
| XY AMR | 1.825562 | 2.003767 | 0.03103141 | 0.04944562 |
| RZ uniform | 1.991655 | 1.990270 | 0.03657952 | 0.04028410 |
| RZ AMR | 1.769790 | 1.958653 | 0.03154876 | 0.05239670 |

donor 仍走原同 F/质量的独立路径；从24份原始 edge/node 快照重新组装实际 q、R、rhs、D 与材料/边界功。最大独立尺度误差 **1.486e-14**，小于原 1e-12。完整 CSV、实际退出码、源码和运行二进制在冻结证据中绑定。

## 工程、成本和能力边界

标准 Release 六目标构建通过，必要 **8/8 CTest** 通过；原124场景和来源拒绝保持，新增12组 XY/RZ × uniform/AMR × rho1/1e3/1e6 的实际闭合端口、独立矩/载荷、常数、transverse affine 与错误端口拒绝通过。原 selected/tiny/counterflow/三层移动几何验证继续覆盖原一般 trace 入口，不能据此声称这些场景的新闭合点 PDE 都已验收。

准备成本仍偏高：OMP=1 时 n32四组合的准备和瞬时逐行审计总计 **447.69秒**；XY/RZ AMR 分别7,351/8,051次ADMM、32/21次活动约束精化。完整研究总计 **928.29秒**，包含重新准备和全部时间步。这是联合稀疏求解，不是每个 patch 独立常数成本的小求解；尚未提供步间 warm start/几何因子复用或动态材料更新时间基准，不能声称已达到大规模生产效率。

不透明外部 trace 的原 H10c 构造保持，**不具有本次闭合入口的全 native 有界性保证**；一般生产边界归 H14/H15。新闭合点 PDE 强制一次矩目前验证光滑纯材料族，不把它推广成任意 sharp-cut 精确 affine 条件，更不会反向修改材料或 FirstJet 以满足它。H08a 必须继续在真实非均匀材料/时空来源下验收 R 余项和相容性。旧 H10b/c/m 报告保持原 FAIL，不追溯改判。

一次中途重链接因原 audit32 可执行文件仍运行而出现 LNK1104，失败日志保留；该数值进程正常结束，没有超时重启。随后最终标准构建和上述验证完成。临时独立优化器安装及中间输出留在 runs/checks，不增加项目依赖、target、永久工具或生产变体。

## 交付

[冻结基线](../../runs/baselines/H10o_20261005/README.md)包含实际 source/runtime、原始结果、设计前快照、独立重组脚本、范围核验及 manifest；当前检查目录为 runs/checks/H10o/20261005_084402_075。1,044源码成员中仅2 production、1 test、1 research变化，原 App/PDF/输入/CMake及材料 owner 保持。API见[最终 native 对流](../api/H10c_Cch_Native_Convection_Moments.md)。总体项目 goal 仍 active；本轮没有执行 H08a 或 H11。
