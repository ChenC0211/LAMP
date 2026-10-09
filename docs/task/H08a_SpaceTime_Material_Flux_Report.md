# H08a / 同源时空相通量积分报告

2026-10-06，**PASS：本任务规定的积分、材料一致性与收敛范围**。r11 基础设施保持封口；没有另立 source/Jacobian/低秩任务。完整 CCH 流体 App、一般 sharp-cut PDE、全耦合时间二阶及效率不在本次通过范围内。

任务书：[H08a](H08a_SpaceTime_Material_Flux_Upgrade_TODO.md)；实现前范围及原门槛：[设计](H08a_SpaceTime_Material_Flux_Implementation_Design.md)。本报告不改判历史失败。

## 唯一积分与独立几何

`CchMaterialTransport::prepareSpaceTime` 从同一原始 Pure FirstJet、同一完整 endpoint characteristic source 生成每个采样时刻的真实几何。两端 B 的原系数直接定义步内场；每个时间节点调用同源空间相率，再沿真实 canonical atom 做 Gauss 3/5 自适应时间积分。不是三个独立 snapshot 的插值，也不是中点相率乘 dt。

液气积分仅密度加权一次为 F_rho；primal/dual 复用同一 atom 积分及原 incidence。平均通量率为 F/dt，中点即时率只作诊断。几何质量来自各自真实 publication；步积分缺陷 r_int=质量差+有向材料积分和，缺陷率 R=r_int/dt，继续原 uR、材料动能与 KEP/skew。Pure FirstJet 不读 R、不接受材料修复，无 M_K/bubble，无 R=0 或 quarter commuting 强制。

Interface 只增加中立已解析运动消费，保留原 RK2 入口。给定整条外边界 FirstJet 是显式数据规则；内部越界且缺少可用历史仍拒绝。预算不足、错误来源/完整 clock 或旧网格来源明确拒绝，无中点/旧守恒路线 fallback。

## 必要修正及原失败

冻结 r11 程序在内部主顶点的零场、轴向和对角常平移上仍拒绝，虽实际连接均有完整 C1 证书。原输入/旧程序/失败记录保留；当前仅通过真实 C1 面连通选择合法 chart，保持 point/J、三阶公式、原根/工作预算。此修正直接服务 FirstJet 消费，不扩大 source 精修范围。

端点穿面还复现了 RZ Green 测量的舍入容量越界。只把已有 Green 势的自由坐标原点选在靠近界面的既有 box 角点，标量、根、图、区域和求积规则不变；体积不裁剪，原严格容量/误差预算保留。采用连续内部原点的中间候选曾触发原 phase-region 失败，故未接受。重网格测试最初误建非 2:1 网格，修正为原 FirstJet 平衡嵌套 fixture 后通过；这些失败均保留在检查目录。

生产改动共 10 文件，原 characteristic/material 两份测试扩展，完整研究归 `cases/benchmarks/cch_spacetime_material`；临时脚本只归 `runs/checks/H08a/20261005_232800`。未改 regular 准备权重、最终 native 对流算子、压力或黏性。

## 当前程序验证

| 项目 | 结果 |
| --- | --- |
| 标准 Release、架构/目录及受影响契约 | 24 项 PASS；最新 11 项架构/材料 + 保持不变的 11 项数值契约 + 当前 vertex/edge 两项 |
| 移动平面，四布局 × 密度比 1/1e3/1e6 | 原 12 组及实际粗细边界切换 6 组，解析通量最大相对误差 3.227821e-08 ≤ 1e-7 |
| 实际粗细边界非零相通量切换 | XY/RZ × 三密度，24 个共享子面检查；两侧 level 不同，面相从全气变全液，沿用原通量预算 |
| 起止端点穿面 | 8 组，最大误差 1.766975e-15 ≤ 1e-7 |
| 同源 dense B/J | 8448 比较，最大 6.192315e-16 ≤ 1e-12 |
| 实际内部 FirstJet 骨架轨迹 | 3268 条，point/J 最大 4.574751e-15 ≤ 1e-12 |
| 原 characteristic 验收 | 24 时间族及原 FD/分区/可信参考/负控保留；最小 point/J 阶 2.998745/2.998115 |
| atom/primal/dual、真实 R/uR 与材料 KE/skew 独立复算 | 22 组，最大 1.071199e-16 ≤ 1e-12 |
| 实际 regrid 后 transfer/rebuild/resume | XY/RZ × uniform/AMR 四组通过；旧来源拒绝，无库存 |
| 原 H10b 完整回归 | 28 窗口、15116 步、84 PASS；四旧 donor FAIL 保留，8 份结果表与 r11 逐字相同 |

H10b 原 24 份 donor 功快照独立重组通过，同一 F、R/uR、耗散及材料/边界功保持原 1e-12 门槛。结果见[标准运行](../../runs/benchmarks__cch_convection/20261006_012544_870_H08a_corner_ec79cee7/run.json)。

## 固定窗口空间研究

非多项式标量 exp(r²−a²)−1 的平移圆/轴上球，T=0.04、rho 比 1e6、n=4/8/16。每步消费上一真实 end publication，不用解析解重置内部状态。参考液相通量为独立 swept-box 解析体积；误差参考不参与生产接受。

下表为整个固定窗口内 L1 缺陷率除以固定物理输运尺度 rho_l U A_projected；归一尺度不随面数量增长。

| 布局 | n4 的 R | n8 的 R | n16 的 R | 最后空间阶 | n16 最大 dual R |
| --- | --- | --- | --- | --- | --- |
| XY / uniform | 1.873755e-04 | 6.852458e-06 | 2.934472e-07 | 4.545 | 3.700592e-07 |
| XY / 2:1 AMR | 4.921339e-05 | 4.870454e-07 | 4.642264e-08 | 3.391 | 2.481005e-08 |
| RZ / uniform | 1.093568e-04 | 8.392932e-06 | 4.119707e-07 | 4.349 | 4.218857e-07 |
| RZ / 2:1 AMR | 1.091624e-04 | 8.324407e-06 | 4.119866e-07 | 4.337 | 4.219018e-07 |

四布局的通量、几何测量误差以及 primal/dual R 均随空间加密下降，最细网格满足原 1% 物理输运尺度门槛。R 保持非零。RZ 球在 AMR 布局中没有穿过细化 patch；平面短契约覆盖真实粗细穿面，不能将本曲面结果扩大为任意 RZ sharp-cut AMR 验收。

## 时间自收敛及误差底限

固定 n8 和 T=0.04，取 1/2/4/8 个真实连续步，8 步作同网格时间参考。分别复算窗口平均缺陷率、对齐子区间的积分差、逐子步常值缺陷率历史和终态几何体积。下表历史误差为时间积分的 L1(R_dt−R_ref)/T（有量纲质量/时间）；所有实际步都参与，不靠整窗正负抵消证明收敛。

| 布局 | 1/2/4 步的缺陷率历史误差 | 历史误差末档阶 | 窗口平均误差末档阶 |
| --- | --- | --- | --- |
| XY / uniform | 1.305989e-01 / 5.704740e-02 / 2.051721e-02 | 1.475 | 1.627 |
| XY / 2:1 AMR | 5.164453e-02 / 2.454268e-02 / 1.056386e-02 | 1.216 | 1.654 |
| RZ / uniform | 4.992303e-02 / 2.199559e-02 / 8.599179e-03 | 1.355 | 1.637 |
| RZ / 2:1 AMR | 4.962233e-02 / 2.185437e-02 / 8.537537e-03 | 1.356 | 1.635 |

这些自误差均下降。固定 h 的绝对 R 随时间加密趋向非零空间投影余项，并不要求降至零。历史范数使用逐步常值时间表示，其阶数不能直接当作完整流体时间推进阶数；8 步参考也不是连续精确解。本项未证明完整耦合二阶，不能据此新增 source/Jacobian 精修先决任务。

## 成本、复现与来源

最新 28 窗口研究耗时 4167.68 秒，单线程、非独占计时；最大单区间几何采样 1444，最大最终时间面板 105。自适应误差量是数值指标，不是严格区间包络。完整 publication 在多个时间节点的成本仍明显，本次不宣称效率通过。

```powershell
.\tools\build.cmd -Preset local-benchmarks -Target LAMP_cch_spacetime_material_benchmark
.\tools\run_case.cmd benchmarks/cch_spacetime_material -Threads 1
```

[最新标准运行](../../runs/benchmarks__cch_spacetime_material/20261006_012544_866_H08a_corner_547e7f08/run.json)保留程序/输入及完整输出；程序 SHA256 `DE983B74BB5C7579CFAA857C23BD36F3D12955EA16E74A497B5308AB8285F0EF`。[独立复核](../../runs/checks/H08a/20261005_232800/current_complete_verification.json)绑定 1386 源码/输入和当前受影响测试，[冻结证据](../../runs/baselines/H08a_20261006/README.md)保留源码、程序及原失败。原 Green 算术上的首轮 28 窗口仅作前版比较，没有替代本次运行。

H08a 的规定积分目标完成；r11 继续作为最后基础设施关口。后续按路线图单任务推进 Hydro，材料四项硬约束保持，完整流体时间步/边界及 EHD 的未完成目标仍保留。
