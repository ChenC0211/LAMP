# H08a2 / 实施前设计变更：实际分片B的事件化midpoint来源

2026-10-05，代码前固定；**本候选最终FAIL并退役**，见[报告](H08a2_Piecewise_Characteristic_Source_Report.md)。原Core连续速度RK2和Interface生产入口不变。[任务](H08a2_Piecewise_Characteristic_Source_TODO.md)。下述原设计及门槛保留，不作为当前生产API。

## 来源及文件预算

允许5 production/build文件：原 `CchSpatialAdvector.h/.cpp` 补充实际分支闭包及单侧查询；新 `CchCharacteristic.h/.cpp` 拥有prepared几何查询与不可变时间来源；根CMake显式注册唯一库owner。允许2 test/build文件：新 `cch_characteristic_contract.cpp` 和原tests/CMakeLists。没有App、字典、F/R或Interface改动。

几何准备保留同一advector reconstruction，索引其实际cell/knot矩形，用包围盒树查找闭包；不制造新的mesh/field。时间来源保留两个同prepared重构、同边界law的真实advector实例，严格核验previous Initial/End、arrival End、step、begin/end和geometry identity；采用线性物理时间插值，value、空间J、时间导数来自这两个相同端点。endpoint值相同即支持冻结场，不把它声称为一般时间精确。没有任意callback或重新拟合速度。

## 离散运动与事件

使用现有显式midpoint的二阶连续扩展；改变的是跨不连续分支的事件语义，不是另一个平滑速度。以反向elapsed e、反向速度 w(e,x)=−B(t_arr−e,x) 记：

```
w0 = w(e,x)
wm = w(e+h/2, x+h*w0/2)
P(s) = x+s*w0+(s*s/h)*(wm−w0), 0<=s<=h
P(h) = x+h*wm
```

中间采样必须落在所选实际分支闭包；必要时二分h，不能外延采样。maxStep由调用者显式给定、与起点坐标无关；remaining与二分后的h依赖elapsed，其导数必须传播。对四条矩形边界分别解这个二次dense曲线的全部实根，选择第一正向外穿根；不会只看最终点而遗漏dense曲线先出后入。无事件时就是原midpoint离散步。事件点的法向坐标取真实边界坐标，根残差按浮点尺度核验，不用坐标扰动改变拓扑。

**微分同一离散映射**：传播位置导数D、elapsed导数E、h导数H；端点插值的w_e也参与。令Q=∂P(s)/∂x_arr（固定s，包括h的依赖），v=∂P/∂s。边界n·P=c给 `S=−nᵀQ/(n·v)`，随后 `D_new=Q+v*S`、`E_new=E+S`。最终固定时刻的remaining导数是−E。不能只乘连续ODE saltation矩阵替代本离散事件导数；也不能把事件时间当常数。

外部或selected-domain边界命中返回真实法向/位置/时间，无内部Jacobian。内部横穿使用同一几何树找后继分支，需唯一并指向其内部。相邻法向一致性由原Hdiv来源保证，不以修正B满足。零时间不可解析、掠切根、同时多面事件、后继不唯一明确拒绝。零轴法向允许沿轴的单侧轨迹，不允许负r或反射。

在跳跃面上的起点，流映射通常只存在单侧导数。结果标明单侧范围；可显式指定实际初始分支，否则只接受唯一指向内的分支。不能把单侧结果伪装成全局可微回调。终点恰落内部事件也有单侧约定，并在固定终时正确补偿elapsed导数。这些限制后续FirstJet耦合必须处理，不能在本项自动消失。

## 代码前验收门槛

T0：标准Release；新characteristic及原architecture/self-test/layout、phase/momentum/material/dual/KEP共9项短CTest；原5数值契约确定性输出保持（计时除外），5/2范围核验。

T1：首先真实原生场 `u_x=1, u_y=分片线性x`（由原native均值/同一重构产生）跨切向跳跃，独立分段积分的point/J误差≤1e−11；至少一个真实cell事件和一个AMR内部knot事件。一般静态/时间变化源的同一完整离散映射由独立中心差分核验，scaled Jacobian误差≤5e−7。XY/RZ、uniform/2:1、n4/8/16、同cell/单/多事件、RZ轴、selected/孔洞/不连通、外部命中、单侧起点和拒绝/复制/lifetime/regrid明确覆盖。无事件连续场与原Core单步结果/J≤1e−12。

T2：上述四网格族、n4/8/16、静态及线性变化的同一可识别affine场，固定T=0.25，maxStep=T/20、T/40、T/80。32个内部初值，值/J对解析流误差RMS下降，末两档阶均≥1.8，最细scaled RMS≤1e−4；报告事件数量和误差，不能删掉穿面点来通过。该项是运动学时间验证，不冒充父H08a的相积分/R收敛或完整流体时间推进。

门槛失败即本项FAIL。若发现新数学缺口，保留实际反例并据其选择下一动作，不回退到三个snapshot插值、平滑B、守恒库存或冻结速度冒充一般来源。
