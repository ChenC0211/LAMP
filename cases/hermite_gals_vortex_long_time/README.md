# 可逆差速旋涡长时间验证

独立规定速度 Interface 基准，使用当前五文件 schema；用于检验非线性变形、两次形状回归和重复 AMR 转移的误差积累。数学与生产实现不修改，解析参考只读检查保存结果。该窗口可提供论文验证图与原始数据，但单次长运行不代替完整网格/时间收敛或真实颈缩、破碎、Hydro/EHD 验收。

## 配置与命令

| profile | 配置 | 目的 |
| --- | --- | --- |
| dynamic_amr | 基础 16×16，中央 level-1 加密；每 125 步交替细化/粗化，共 15 次 | 长时间重复网格转移 |
| fixed_amr | 同样的基础网格和初始加密，始终保持 448 个叶单元 | 固定 AMR 对照 |

两组均为 dt=0.001、endTime=2、2,000 步，每 25 步输出场（含初末态，共 81 帧）。dynamic 在 256/448 个叶单元之间切换，因此对照检验的是“网格策略与转移的综合影响”，不能把差异全部归因于 transfer。

```powershell
.\tools\run_case.cmd hermite_gals_vortex_long_time/dynamic_amr
.\tools\run_case.cmd hermite_gals_vortex_long_time/fixed_amr
```

结果分别在 `runs/hermite_gals_vortex_long_time__dynamic_amr/<时间戳>` 和 `runs/hermite_gals_vortex_long_time__fixed_amr/<时间戳>`。ParaView 打开对应 output/fields/solution.pvd。

## 独立参考与预先验收

设 a=x−0.5、b=y−0.5、ρ²=a²+b²：

\[
(u,v)=40\cos(\pi t)\rho^2(-b,a),\qquad
\beta(t)=\frac{40}{\pi}\sin(\pi t).
\]

粒子的 ρ 恒定，角度增加 βρ²；此映射面积 Jacobian 为 1，液相准确面积始终为 π·0.13²。初态圆心 (0.65,0.53)，半径 0.13；在 t=0.5/1.5 变形最大，在 t=1/2 准确回到初态。一个速度周期含两次方向相反的变形—回归过程。

参考验证器使用逆旋转解析 level set 的独立二阶链式求导，并以 forward 参数曲线法向/曲率和 finite difference 交叉检查；不导入 Core 表达式或数值求解器。每个 profile 的 reference/verification.json 保存研究参数、五输入哈希和运行前固定的阈值：面积相对误差峰值≤0.001，法向位置缺陷/R≤0.02，法向 L2≤0.05，曲率相对 L2≤0.1，末态位置缺陷/R≤0.01。参考契约独立于生产字典；修改研究输入必须显式制定新的参考与哈希，不能自动放宽阈值。

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File tools/postprocess/analyze_interface_reference.ps1 -Run runs/<算例>/<时间戳>
```

验证器自动使用本机 ParaView 的 pvpython（NumPy/Matplotlib），或用 -Python 指定满足依赖的 Python。它在该运行 output/analysis 下建立带时间戳的新分析目录，保存 accuracy_history.csv、verification.json、REPORT.md 和 PNG/SVG 图。退出码 0=PASS、1=精度/完整性 gate FAIL、2=来源/结构/工具错误；失败运行须显式 -AcceptedPrefix，不能标为完整窗口 PASS。

面积比较覆盖全部接受态；位置/法向/曲率比较覆盖保存的 81 帧，采用 VTP polyline 弧长权重。位置缺陷为 |phi_exact|/|grad(phi_exact)|，不是严格欧氏距离；另报告双向采样点到折线的 Hausdorff 近似。法向/曲率参考在数值采样点求值；不把这些指标称为完整 authoritative geometry 的连续误差范数。形状图直接读保存的界面，仅参考曲线由解析公式生成。

当前数据可用作后续原生 CCH 接线前的 Interface 基线；流体和电场尚不在本算例中。

2026-10-03 两组 2,000 步完整窗口和独立参考验证均 PASS；实际曲率、面积、耗时及原始数据见 [验证报告](../../docs/task/Interface_Long_Vortex_Verification_Report.md)。曲率相对 L2 峰值为 0.174%（dynamic）/0.103%（fixed），末态为 0.0483%/0.0366%。
