# H10l1 / 实施前实际速度场 closure 与边界设计

2026-10-04；前一轮progress。启动复核H10l的328 payload/115 current canonical及k/j/j0/i2/f/e/d/c归档通过。按[任务书](H10l1_Spatial_Advector_Closure_And_Boundary_TODO.md)只改现actual advector h/cpp、现phase/material两tests及既有research observer/文档，0新源码/targets/tools。H10l失败候选不恢复；H09/KEP与原H10b gates保持。

## 先定位原反例

先为原evaluate拒绝增加实际cell、原D_I/V、实际field.div及native gross上下文；只改诊断，不改计算。既有research driver只读输出n4四布局原单面basis的实际输入、patch系数、原DI/physical measures和拒绝。由原片段算术重算受影响cell，与真实拒绝比较；独立解析完整face均值是输入均值，零均值线性变化的全face积分必须为零。区分其分段相消舍入与真实系数/几何错误，不能只因旧gate使用1e−300就放宽标度。

## 零均值基函数与积分

若确认反例来自零均值线性迹的分段相消，采用同一真实face的局部体积切向坐标：端点t0/t1、L=t1−t0、s=t−t0，B=b+k(s−L/2)。b由原native物理均值及RZ因子生成；k仍由原one-ring weights生成，不换stencil或目标。

对从face起点到长度w的prefix，积分为 `w*(b+k*(w-L)/2)`。完整face w=L时切向贡献精确为0；不先积出两份大数再相消。每个cell两侧分别按真实face计算prefix并取反号差，除原normal separation。两侧各自分段即精确积分同一piecewise trace，原atom overlay仍按两侧knots分割。原lift继续使用这些prefix与实际平均导数；没有把field.div赋成expected D_I，没有第二散度矩阵。

原D_I/V仍独立核验，1e−12 native gross阈值保持；全零gross要求精确零，去掉经验1e−300 floor。独立atom积分和cell边界积分用真实polynomial work归一化，避免以应为0的signed sum当误差尺度；不要求一般零均值face的子段通量为零。tiny/large动态范围明确验收，有限值无法表示时拒绝，不填零。

## 显式不可穿透来源

在同一header/cpp定义只读`CchAdvectorBoundary`：绑定原momentumGeometry instance和调用者明确列出的native velocity DOFs。仅允许该selected domain的真实boundary（包含明确选定的selected边界及axis），重复、越界、内部面拒绝。law是静态齐次normal条件，没有独立clock或速度值；每次实际evaluate验证相同geometry/mesh和指定面native mean严格0。

现四参数evaluate保持无law语义；新typed重载接受law，原Data保留只读绑定副本。每个actual B revision仍取得唯一advector stamp并携带原geometry/full clock，下游source不需要额外伪造标识。law副本共享身份，失效mesh/regrid/wrong instance拒绝。

只在明确指定面把normal slope固定0，配合严格零输入均值，使整个normal polynomial为0；其余face，包括mean==0的内部/一般边界，仍保留原重构。faceB、cell lift、native/cell velocity、canonical atoms及phase/material consumers全部消费同一实际field；没有material/H09后置压F。不可穿透绑定不供外部相数据，也不绕过原zero-trace外部流入拒绝。

## 固定验收与边界

现phase test增加全部非轴basis/线性组合的D_I、独立physical mean/分段共享、cell实际速度散度检查；跨XY/RZ uniform/2:1、零净/近相消、全零、tiny/large输入。现material test增加原四closed场景的明确law与未绑定对照：所有boundary atoms/实际subsegments严格零，原D_I与B散度一致，原Flux仍唯一，内部counterflow保留，rho1/1e3/1e6真实phase/材料率及incidence核验。完整来源/clock/instance/copy/lifetime/regrid/nonfinite及selected拒绝保留。

既有64布局/128advector更新和64布局材料研究继续检验原56 tangent缺口、同源相积分/真实AMR子面和；内部affine重跑，原full/RK4/donor/gates不改。最终必要8项CTest与3项文档检查、源码预算/compiled runtime/独立CSV/flat manifest。若必要身份或支持范围不满足即FAIL并保留证据，不借后续任务或放宽gate关闭本项。通过后另起H10l重试，不在本任务修改dual/KEP或进入H11/H08a。
