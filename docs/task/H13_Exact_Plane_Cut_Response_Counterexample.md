# H13 / 精确平面 cut 与高密度法向响应反例

2026-10-09 / 035647：**当前联合候选仍FAIL，不具备采用条件，H14未启动**。同一H13任务的顺序诊断；0 production/permanent test修改。原静滴六组回升和n8完整周期结果保留，本反例补充当前方法在一般cut、高密度、无黏目标下的缺口。

## 输入与连续真值

单位方形闭域，液体在下、气体在上；精确平面FirstJet `phi=y-H`，`H=.5+theta/n`，theta取0/.25/.5。几何曲率为零，制造法向应力 `J=-cos(k*x)`、k=2*pi。这是瞬时无黏法向应力响应，不是将制造应力作为真实曲率，也没有接受的界面更新。应力和响应都线性，可缩放到任意小幅值。

令rho_g=1、rho_l=rho，法向加速度真值为

`a_Gamma(x)=A*cos(k*x)`，

`A=-k/[rho*coth(k*H)+coth(k*(1-H))]`。

液、气压力分别为

`p_L=-rho*A*cos(k*x)*cosh(k*y)/[k*sinh(k*H)]`，

`p_G=A*cos(k*x)*cosh(k*(1-y))/[k*sinh(k*(1-H))]`。

二者各自满足Laplace方程；`a=-grad(p)/rho_phase`，上下及两侧外边界normal加速度为零；界面normal连续，`p_L-p_G=cos(k*x)`。界面法向运动真值唯一。bulk切向加速度存在两相跳跃，符合此无黏制造解；审计只用连续的界面normal作为运动真值，并分别记录原生面平均和common source的查询，避免把bulk切向跳跃当作source必须逐点再现的条件。

C++从上述真值的解析反导数生成native face averages。独立Python用80位压力/PDE/跳跃核验和16点Gauss积分重算face mean；不调用生产source/tracer生成真值。独立cell/dual相面积按精确平面计算，原实际alpha、mass逐位吻合；没有第二份接受几何或独立库存。

## 实际算子与完整参数积

调用原CchPhaseGeometry、CchMomentumMass、CchPressureLoad、共同source及`assembleCchInterfaceStress`，用原public saddle解 `M*a+C*p+F=0`；strict mu=0，没有黏性。n8/16/32/64 × 三个theta × 密度比1/1e3/1e6，共36组自然完成，246.18秒，原saddle所有门槛通过。

每组分别记录：实际native加速度，真native face averages经过同一source的界面query，实际求解结果经过同一source的界面query，原cell/face/q/pressure/mass/force。前者检验native负载响应；第二项单独检验插值，排除了压力、曲率和求解器误差。

## 反例

| n | 对齐，theta=0，rho=1e6的模态比 | quarter cut，theta=.25的模态比 | middle cut，theta=.5的模态比 |
| --- | --- | --- | --- |
| 8 | 0.956818 | -122.136 | -771.629 |
| 16 | 0.987525 | -26.3790 | -203.318 |
| 32 | 0.996756 | -5.22680 | -46.4026 |
| 64 | 0.999181 | -0.555891 | -10.4379 |

模态比为 `int a_source*cos(k*x) ds / (A*int cos^2(k*x) ds)`。负号表示实际normal模态与解析回复方向相反，而非频率定义差异。

代表性细网格证据为case36：n64、theta=.5、H=.5078125、rho=1e6。解析A=-6.26194239918e-6，实际source模态+6.53616810766e-5；actual normal点误差/|A|为681.172，native face误差/|A|为1536.30。最差native face6337为水平面，y=.53125、x区间[.015625,.03125]，真实质量.000244140625，纯气相；实际a=-.009625576902，真值=-5.340060021e-6。这不是仅靠液相质量变大就能压低的模式。

即使输入真native face averages，当前source仍有独立精度缺口：theta=.25/.5的normal点误差末档阶为0.99556/0.99522，内部x∈[.25,.75]的末档阶0.99258/0.99798；normal Fourier模态也约一阶。n64/middle cut的真值插值最大相对误差5.16228%、内部3.69322%、模态比.962521。对齐平面的真值插值约3.93阶、actual模态约1.985阶。良好的对齐端口结果不能覆盖一般cut。

## 应力配对与压力商空间

记 `g_u=sum_q W_q*k_q*b_qu`，`V_u=sum_q W_q*b_qu`，`a=C alpha`。当前实际力可以写成

`F=g-diag(kbar)*(a+V)`。

因此局部锚点引入的体积余项不是自动的pressure gradient。case36中

`g^T a_native = +3.26808405e-5`，

`-F^T a_native = 3.25406049e-6 = a_native^T M a_native`（数值余项保留），

`a_native^T C p = -4.38599e-24`。

这里是力—加速度配对，不是一次实际时间步的能量。surface gradient与实际加速度的正配对与回复方向相反；压力不做功，不能用压力解误差解释。局部体积余项配对约-3.59349e-5，是原算子中的真实项，不隐去或令R=0。单纯通过同一公式的装配/虚功短契约，不证明这一项在原真实M/C下具有正确的物理回复作用。

没有用此瞬时诊断宣称实际时间推进增长率或H14稳定性；实际capillary的FirstJet曲率变化与全state谱仍需原范围验证。但当前候选的准确性和一般高密度目标已被直接反例否定。

## 独立核验

36组的解析PDE/跳跃80位残差最大4.22e-81；C++/独立16GL face mean差归一化5.83e-15；真native离散连续性1.48e-15；alpha/mass差0；plane坐标/长度差0、normal4.44e-16、曲率2.27e-13。独立原C*p差0，free momentum相对残差1.99e-14；固定normal壁面的反力显式保留，不把壁面反力行误算为自由动量失败。故意污染真值幅值被拒绝。

数据/程序：[flat_cases](../../runs/checks/H13/20261009_035647/flat_cases.csv)、[flat_audit](../../runs/checks/H13/20261009_035647/flat_audit.json)、[实际原算子程序](../../runs/checks/H13/20261009_035647/flat_cut_response.cpp)、[独立审计](../../runs/checks/H13/20261009_035647/audit_flat.py)。完整fields/quadrature/pressure及原运行日志在同一证据目录。

最小复现：n8、theta=.5、rho=1e6，原36组case9。仅64个cells；单独命令运行1.03秒，输出case编号改为1，全部cell/face/q/summary数据与case9逐字一致，独立审计通过。[单算例数据](../../runs/checks/H13/20261009_035647/flat_minimal/flat_audit.json)与[逐字核对](../../runs/checks/H13/20261009_035647/flat_minimal/full_case9_replay.json)保留。在新的D:运行目录执行：

```powershell
& 'D:/Chen/LAMP/runs/checks/H13/20261009_035647/build/Release/flat_cut_response.exe' 8 0.5 1000000
& 'D:/Chen/LAMP/.venv/Scripts/python.exe' 'D:/Chen/LAMP/runs/checks/H13/20261009_035647/audit_flat.py' .
```

## 下一步的实质约束

当前授权的normal值/切向导数面矩改造不足以闭合一般cut；需要先推导共同source的内部normal trace与原真实M/C/应力测试的相容关系，分别证明真native输入的normal精度和实际应力响应的密度/位置一致性。不能只改局部锚点、靠对齐平面的好频率、删弱项、增加黏度或将制造解数据写进production。

原[联合授权](H13_Native_Flux_Port_Method_Change.md)限定master normal值/切向导数规则，保持原内部/横向/mixed和其他冻结算子。涉及这些规则的实质扩展，须先形成具体方法、可比性、预算与原全范围验收，再依总纲§2.1作方法决策；自动续跑不构成授权。本证据不预先选择phase绑定、Hodge变更或另一条production source。两个已提出的KEP/reference预算补充仍独立待批，其通过不能迁为该联合候选采用条件。
